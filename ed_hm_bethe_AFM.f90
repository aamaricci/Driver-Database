program lancED
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                       :: iloop,Le,Nso,iorb,ispin
  logical                                       :: converged
  real(8),dimension(5)                          :: Wbethe,Dbethe
  integer                                       :: Nineq,Nlat
  !Bath:
  integer                                       :: Nb
  real(8),allocatable                           :: Bath(:),Bath_prev(:)
  !  
  complex(8),allocatable                        :: Hloc(:,:,:,:)
  real(8),dimension(:,:),allocatable            :: Dbands
  real(8),dimension(:,:),allocatable            :: Ebands
  real(8),dimension(:),allocatable              :: H0
  real(8),dimension(:),allocatable              :: de,dens
  real(8),dimension(:),allocatable              :: Wband
  character(len=16)                             :: finput
  real(8)                                       :: wmixing
  !
  !The local hybridization function: [Nlat/Nineq,Norb,Norb,Nspin,Nspin,:]
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss,Weiss_prev
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal
  !
  integer                                       :: comm,rank
  logical                                       :: master
  logical                                       :: betheSC,wGimp,mixG0,symOrbs,reGF

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  Nlat =2                       !two sub-lattices
  Nineq=1                       !only one is inequivalent (B=-A)


  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(Wbethe,"WBETHE",finput,default=[1d0,1d0,1d0,1d0,1d0])
  call parse_input_variable(Dbethe,"DBETHE",finput,default=[0d0,0d0,0d0,0d0,0d0])
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(betheSC,"BETHESC",finput,default=.false.)
  call parse_input_variable(wGimp,"wGimp",finput,default=.false.)
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  call parse_input_variable(symOrbs,"symOrbs",finput,default=.false.)
  call parse_input_variable(reGF,"REGF",finput,default=.false.)
  !
  call ed_read_input(trim(finput),comm)

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Nlat,"nlat")
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(Nineq,"nineq")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(Nspin/=2.OR.Norb>5)stop "Wrong setup from input file: Nspin/=2 OR Norb>5"
  Nso=Nspin*Norb

  allocate(Ebands(Nso,Le))
  allocate(Dbands(Nso,Le))
  allocate(Wband(Nso))
  allocate(H0(Nso))
  allocate(de(Nso))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))

  !Because we solve only one sublattice and retrieve the other from the solved one,
  !we only need lattice information (H(e),etc.) for this one sublattice, say A
  Wband = Wbethe(:Norb)         !band width
  H0    = Dbethe(:Norb)         !Local energy
  do iorb=1,Norb
     Ebands(iorb,:) = linspace(-Wband(iorb),Wband(iorb),Le,mesh=de(iorb)) !dispersion
     Dbands(iorb,:) = dens_bethe(Ebands(iorb,:),Wband(iorb))*de(iorb)     !DOS
  enddo
  if(master)call TB_write_Hloc(one*diag(H0))
  Hloc(1,1,:,:)=diag(H0)        !spin up
  Hloc(2,2,:,:)=diag(H0)        !spin dw = spin up, no symmetry-breaking yet

  !Allocate all Fields:
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Weiss_prev(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(dens(Norb))



  !setup solver, for the ineq sublattice only
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_prev(Nb))
  call ed_init_solver(Bath)
  call ed_break_symmetry_bath(Bath,sb_field,1d0)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     !Solve sublattice A only:
     call ed_solve(comm,Bath,Hloc)

     !Retrieve self-energy for A-lattice only:
     call ed_get_sigma_matsubara(Smats(1,:,:,:,:,:))
     call ed_get_sigma_realaxis(Sreal(1,:,:,:,:,:))   
     !Impose Neel AFM symmetry: B_sigma = A_{-sigma}
     do ispin=1,2
        Smats(2,ispin,ispin,:,:,:)=Smats(1,3-ispin,3-ispin,:,:,:)
        Sreal(2,ispin,ispin,:,:,:)=Sreal(1,3-ispin,3-ispin,:,:,:)
     enddo

     !Evaluate Gloc and update Weiss, all at once here:
     call get_delta_bethe()
     !print:
     call dmft_print_gf_matsubara(Gmats,"Gmats",iprint=4)
     call dmft_print_gf_matsubara(Greal,"Greal",iprint=4)
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=4)
     !
     !
     !
     if(mixG0)then
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_prev
        Weiss_prev=Weiss
     endif
     !
     !Perform the SELF-CONSISTENCY by fitting the new bath, only A sublattice needs to be fitted.
     if(symOrbs)then
        call ed_chi2_fitgf(comm,Weiss(1,:,:,:,:,:),bath,ispin=1,iorb=1)
        call ed_orb_equality_bath(bath,save=.true.)
     else
        call ed_chi2_fitgf(comm,Weiss(1,:,:,:,:,:),bath,ispin=1)
        call ed_chi2_fitgf(comm,Weiss(1,:,:,:,:,:),bath,ispin=2)
     endif
     !
     !MIXING:
     if(.not.mixG0)then
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
        Bath_prev=Bath
     endif
     !
     converged = check_convergence(Weiss(1,1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     call ed_get_dens(dens)
     if(nread/=0d0)call ed_search_variable(xmu,sum(dens),converged)

     call end_loop
  enddo


  call finalize_MPI()

contains


  subroutine get_delta_bethe
    integer                         :: i,j,ie,iorb
    complex(8)                      :: iw,zita(2,2),zeta(2)
    complex(8),dimension(2,2,Lmats) :: gloc
    complex(8),dimension(2,2,Lreal) :: grloc
    real(8)                         :: wm(Lmats),wr(Lreal)
    real(8)                         :: epsi(Le),dos(Le),wb
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    !

    do iorb=1,Norb              !work out each orbital independently:ACTHUNG
       do i=1,Lmats
          iw      = xi*wm(i)
          zita(1,1) = iw + xmu - Smats(1,1,1,iorb,iorb,i) !A-up,up
          zita(1,2) = iw + xmu - Smats(1,2,2,iorb,iorb,i) !A-dw,dw
          zita(2,1) = iw + xmu - Smats(2,1,1,iorb,iorb,i) !B-up,up
          zita(2,2) = iw + xmu - Smats(2,2,2,iorb,iorb,i) !B-dw,dw
          !
          zeta(1)    = zita(1,1)*zita(1,2)
          zeta(2)    = zita(2,1)*zita(2,2)
          !
          !G_{lat,sigma} = zita_{lat,sigma'} * Int de D(w)/(zita_{lat,sigma}*zita_{lat,sigma'} - e)
          Gmats(:,:,:,iorb,iorb,i)    = zero
          do ie=1,Le
             Gmats(1,1,1,iorb,iorb,i) = Gmats(1,1,1,iorb,iorb,i) + Dbands(iorb,ie)/(zeta(1) - Ebands(iorb,ie)**2)
             Gmats(2,1,1,iorb,iorb,i) = Gmats(2,1,1,iorb,iorb,i) + Dbands(iorb,ie)/(zeta(2) - Ebands(iorb,ie)**2)
          enddo
          !Update all components:
          Gmats(1,2,2,iorb,iorb,i) = zita(1,1)*Gmats(1,1,1,iorb,iorb,i) !G_{A,dw,dw}
          Gmats(1,1,1,iorb,iorb,i) = zita(1,2)*Gmats(1,1,1,iorb,iorb,i) !G_{A,up,up}
          Gmats(2,2,2,iorb,iorb,i) = zita(2,1)*Gmats(2,1,1,iorb,iorb,i) !G_{B,dw,dw}
          Gmats(2,1,1,iorb,iorb,i) = zita(2,2)*Gmats(2,1,1,iorb,iorb,i) !G_{B,up,up}
          !
          if(cg_scheme=='weiss')then
             Weiss(1,1,1,iorb,iorb,i)= one/(one/Gmats(1,1,1,iorb,iorb,i) + Smats(1,1,1,iorb,iorb,i))
             Weiss(1,2,2,iorb,iorb,i)= one/(one/Gmats(1,2,2,iorb,iorb,i) + Smats(1,2,2,iorb,iorb,i))
             Weiss(2,1,1,iorb,iorb,i)= one/(one/Gmats(2,1,1,iorb,iorb,i) + Smats(2,1,1,iorb,iorb,i))
             Weiss(2,2,2,iorb,iorb,i)= one/(one/Gmats(2,2,2,iorb,iorb,i) + Smats(2,2,2,iorb,iorb,i))
          else
             Weiss(1,1,1,iorb,iorb,i)= iw + xmu - Smats(1,1,1,iorb,iorb,i) - one/Gmats(1,1,1,iorb,iorb,i)
             Weiss(1,2,2,iorb,iorb,i)= iw + xmu - Smats(1,2,2,iorb,iorb,i) - one/Gmats(1,2,2,iorb,iorb,i)
             Weiss(2,1,1,iorb,iorb,i)= iw + xmu - Smats(2,1,1,iorb,iorb,i) - one/Gmats(2,1,1,iorb,iorb,i)
             Weiss(2,2,2,iorb,iorb,i)= iw + xmu - Smats(2,2,2,iorb,iorb,i) - one/Gmats(2,2,2,iorb,iorb,i)
          endif
       enddo
       !
       do i=1,Lreal
          iw=cmplx(wr(i),eps)
          zita(1,1) = iw + xmu - Sreal(1,1,1,iorb,iorb,i) !A-up,up
          zita(1,2) = iw + xmu - Sreal(1,2,2,iorb,iorb,i) !A-dw,dw
          zita(2,1) = iw + xmu - Sreal(2,1,1,iorb,iorb,i) !B-up,up
          zita(2,2) = iw + xmu - Sreal(2,2,2,iorb,iorb,i) !B-dw,dw
          !
          zeta(1)    = zita(1,1)*zita(1,2)
          zeta(2)    = zita(2,1)*zita(2,2)

          !G_{lat,sigma} = zita_{lat,sigma'} * Int de D(w)/(zita_{lat,sigma}*zita_{lat,sigma'} - e)
          Greal(:,:,:,iorb,iorb,i)    = zero
          do ie=1,Le
             Greal(1,1,1,iorb,iorb,i) = Greal(1,1,1,iorb,iorb,i) + Dbands(iorb,ie)/(zeta(1) - Ebands(iorb,ie)**2)
             Greal(2,1,1,iorb,iorb,i) = Greal(2,1,1,iorb,iorb,i) + Dbands(iorb,ie)/(zeta(2) - Ebands(iorb,ie)**2)
          enddo
          !Update all components:
          Greal(1,2,2,iorb,iorb,i) = zita(1,1)*Greal(1,1,1,iorb,iorb,i) !G_{A,dw,dw}
          Greal(1,1,1,iorb,iorb,i) = zita(1,2)*Greal(1,1,1,iorb,iorb,i) !G_{A,up,up}
          Greal(2,2,2,iorb,iorb,i) = zita(2,1)*Greal(2,1,1,iorb,iorb,i) !G_{B,dw,dw}
          Greal(2,1,1,iorb,iorb,i) = zita(2,2)*Greal(2,1,1,iorb,iorb,i) !G_{B,up,up}
       enddo
    enddo
  end subroutine get_delta_bethe



end program lancED

  









