program lancED
  USE EDIPACK
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  !
  implicit none
  integer                                       :: iloop,Le,Nso,iorb,ispin,ilat,jorb,jspin,io,jo
  logical                                       :: converged
  real(8),dimension(5)                          :: Wbethe,Dbethe
  integer                                       :: Nineq,Nlat,Nlso
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
  complex(8),allocatable,dimension(:,:,:) :: Self
  !
  integer                                       :: comm,rank,msize,ierr
  logical                                       :: master
  logical :: getK

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  msize= get_Size_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(Wbethe,"WBETHE",finput,default=[1d0,1d0,1d0,1d0,1d0])
  call parse_input_variable(Dbethe,"DBETHE",finput,default=[0d0,0d0,0d0,0d0,0d0])
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(getK,"GETK",finput,default=.false.)
  !
  call ed_read_input(trim(finput))

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(Nspin/=2.OR.Norb>1)stop "Wrong setup from input file: Nspin/=2 OR Norb>1"

  Nlat =2                       !two sub-lattices
  Nineq=1                       !only one is inequivalent (B=-A)
  Nso=Nspin*Norb
  Nlso=Nlat*Nso



  if(getK)then
     allocate(Self(Nlso,Nlso,Lmats));Self=zero
     if(master)call read_array("Smats",Self)
     call Bcast_MPI(comm,Self)
     allocate(Ebands(Nlso,Le))
     allocate(Dbands(Nlso,Le))
     allocate(Wband(Nlso))
     allocate(H0(Nlso))
     allocate(de(Nlso))
     do ilat=1,Nlat
        do ispin=1,Nspin
           do iorb=1,Norb
              io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin           
              Wband(io) = Wbethe(iorb)         !band width
              H0(io)    = Dbethe(iorb)         !Local energy
              Ebands(io,:) = linspace(-Wband(io),Wband(io),Le,mesh=de(io)) !dispersion
              Dbands(io,:) = dens_bethe(Ebands(io,:),Wband(io))*de(io)     !DOS
           enddo
        enddo
     enddo
     !
     call dmft_kinetic_energy_Bethe(Ebands,Dbands,H0,Self)
     !
     call finalize_MPI(comm)
     stop
  end if



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
  !
  call ed_set_hloc(Hloc)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     !Solve sublattice A only:
     call ed_solve(Bath)

     !Retrieve self-energy for A-lattice only:
     call ed_get_sigma(Smats(1,:,:,:,:,:),axis='m')
     call ed_get_sigma(Sreal(1,:,:,:,:,:),axis='r')
     !Impose Neel AFM symmetry: B_sigma = A_{-sigma}
     do ispin=1,2
        Smats(2,ispin,ispin,:,:,:)=Smats(1,3-ispin,3-ispin,:,:,:)
        Sreal(2,ispin,ispin,:,:,:)=Sreal(1,3-ispin,3-ispin,:,:,:)
     enddo

     !Evaluate Gloc and update Weiss, all at once here:
     call get_delta_bethe()
     call dmft_write_gf(Gmats,"Gloc",axis='m',iprint=4)

     !
     !Perform the SELF-CONSISTENCY by fitting the new bath, only A sublattice needs to be fitted.
     call ed_chi2_fitgf(Weiss(1,:,:,:,:,:),bath,ispin=1)
     call ed_chi2_fitgf(Weiss(1,:,:,:,:,:),bath,ispin=2)
     !
     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath
     !
     if(master)converged = check_convergence(Weiss(1,1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     call Bcast_MPI(Comm,converged)
     if(master)call end_loop
  enddo
  call ed_finalize_solver()


  call dmft_write_gf(Greal,"Gloc",axis='r',iprint=4)


  allocate(Self(Nlso,Nlso,Lmats));Self=zero
  do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
     io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
     jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
     Self(io,jo,:) = Smats(ilat,ispin,jspin,iorb,jorb,:)
  enddo
  if(master)call save_array("Smats",Self)
  deallocate(Self)

  call Barrier_MPI()
  call finalize_MPI()
  stop

contains


  subroutine get_delta_bethe
    integer                         :: i,j,ie,iorb
    complex(8)                      :: iw,zita(2,2),zeta(2)
    complex(8),dimension(2,2,Lmats) :: gloc
    complex(8),dimension(2,2,Lreal) :: grloc
    real(8)                         :: wm(Lmats),wr(Lreal)
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    !

    do iorb=1,Norb              !work out each orbital independently:ACTHUNG
       do i=1,Lmats
          iw      = xi*wm(i)
          zita(1,1) = iw + xmu - Smats(1,1,1,iorb,iorb,i) !Aup
          zita(1,2) = iw + xmu - Smats(1,2,2,iorb,iorb,i) !Adw
          zita(2,1) = iw + xmu - Smats(2,1,1,iorb,iorb,i) !Bup=Adw
          zita(2,2) = iw + xmu - Smats(2,2,2,iorb,iorb,i) !Bdw=Aup
          !
          zeta(1)    = zita(1,1)*zita(1,2) !Z_Aup.Z_Bup = Z_Aup.Z_Adw local
          zeta(2)    = zita(2,1)*zita(2,2) !Z_Adw.Z_Bdw = Z_Bup.Z_Bdw local
          !
          !G_{lat,sigma} = zita_{lat,sigma'} * Int de D(w)/(zita_{lat,sigma}*zita_{lat,sigma'} - e)
          Gmats(:,:,:,iorb,iorb,i)    = zero
          do ie=1,Le
             Gmats(1,1,1,iorb,iorb,i) = Gmats(1,1,1,iorb,iorb,i) + Dbands(iorb,ie)/(zeta(1) - Ebands(iorb,ie)**2)
             Gmats(2,1,1,iorb,iorb,i) = Gmats(2,1,1,iorb,iorb,i) + Dbands(iorb,ie)/(zeta(2) - Ebands(iorb,ie)**2)
          enddo
          !Update all components:
          Gmats(1,1,1,iorb,iorb,i) = zita(1,2)*Gmats(1,1,1,iorb,iorb,i) !G_{A,up,up}
          Gmats(1,2,2,iorb,iorb,i) = zita(1,1)*Gmats(1,1,1,iorb,iorb,i) !G_{A,dw,dw}
          Gmats(2,1,1,iorb,iorb,i) = zita(2,2)*Gmats(2,1,1,iorb,iorb,i) !G_{B,up,up}
          Gmats(2,2,2,iorb,iorb,i) = zita(2,1)*Gmats(2,1,1,iorb,iorb,i) !G_{B,dw,dw}
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







  subroutine dmft_kinetic_energy_Bethe(Ebands,Dbands,Hloc,Sigma)
    real(8),dimension(Nlso,Le),intent(in) :: Ebands    ![Nlat*Nspin*Norb][Lk]
    real(8),dimension(Nlso,Le),intent(in) :: Dbands    ![Nlat*Nspin*Norb][Lk]Lk]
    real(8),dimension(Nlso),intent(in)    :: Hloc    ![Nlat*Nspin*Norb]
    complex(8),dimension(Nlso,Nlso,Lmats) :: Sigma   ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    integer                               :: i,ie
    integer                               :: unit
    real(8),dimension(Norb)               :: Ekin
    real(8)                               :: K
    real(8),dimension(:),allocatable      :: wm
    logical                               :: dos_diag !1. T / 2. F
    complex(8)                            :: iw,zita(2,2),zeta(2)
    integer :: ilat,iorb,jorb,ispin,jspin
    !
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*dble(2*arange(1,Lmats)-1)
    !
    !
    !Get Sigma:
    allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
       jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
       Smats(ilat,ispin,jspin,iorb,jorb,:) = Sigma(io,jo,:)
    enddo


    if(master)write(*,"(A)") "Kinetic energy computation"
    if(master)call start_timer()
    

    do iorb=1,Norb              !work out each orbital independently:ACTHUNG
       do i=1,Lmats
          iw      = xi*wm(i)
          zita(1,1) = iw + xmu - Smats(1,1,1,iorb,iorb,i) !Aup
          zita(1,2) = iw + xmu - Smats(1,2,2,iorb,iorb,i) !Adw
          zita(2,1) = iw + xmu - Smats(2,1,1,iorb,iorb,i) !Bup=Adw
          zita(2,2) = iw + xmu - Smats(2,2,2,iorb,iorb,i) !Bdw=Aup
          zeta(1)    = zita(1,1)*zita(1,2) !Z_Aup.Z_Bup = Z_Aup.Z_Adw local
          zeta(2)    = zita(2,1)*zita(2,2) !Z_Adw.Z_Bdw = Z_Bup.Z_Bdw local
          !
          !G_{lat,sigma} = zita_{lat,sigma'} * Int de D(w)/(zita_{lat,sigma}*zita_{lat,sigma'} - e)
          Gmats(:,:,:,iorb,iorb,i)    = zero
          do ie=1,Le
             Gmats(1,1,1,iorb,iorb,i) = Gmats(1,1,1,iorb,iorb,i) + Dbands(iorb,ie)/(zeta(1) - Ebands(iorb,ie)**2)
             Gmats(2,1,1,iorb,iorb,i) = Gmats(2,1,1,iorb,iorb,i) + Dbands(iorb,ie)/(zeta(2) - Ebands(iorb,ie)**2)
          enddo
          !Update all components:
          Gmats(1,1,1,iorb,iorb,i) = zita(1,2)*Gmats(1,1,1,iorb,iorb,i) !G_{A,up,up}
          Gmats(1,2,2,iorb,iorb,i) = zita(1,1)*Gmats(1,1,1,iorb,iorb,i) !G_{A,dw,dw}
          Gmats(2,1,1,iorb,iorb,i) = zita(2,2)*Gmats(2,1,1,iorb,iorb,i) !G_{B,up,up}
          Gmats(2,2,2,iorb,iorb,i) = zita(2,1)*Gmats(2,1,1,iorb,iorb,i) !G_{B,dw,dw}
          !
          Ekin(iorb) = Ekin(iorb) + Gmats(1,1,1,iorb,iorb,i)*Gmats(2,1,1,iorb,iorb,i)/beta
       enddo
    enddo


    if(master)then
       unit = free_unit()
       open(unit,file="dmft_kinetic_energy.dat")
       write(unit,"(90F15.9)")sum(Ekin)
       close(unit)
    endif
    !
    deallocate(wm)
  end subroutine dmft_kinetic_energy_Bethe


  function anti_diag(V) result(M)
    real(8),dimension(:)               :: V
    real(8),dimension(size(V),size(V)) :: M
    integer :: i
    M=0d0
    forall(i=1:size(V))M(i,Nlso-i+1)=V(i)    
  end function anti_diag
  
end program

  









