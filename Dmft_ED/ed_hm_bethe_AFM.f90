program lancED
  USE EDIPACK
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
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
  integer                                       :: comm,rank
  logical                                       :: master

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(Wbethe,"WBETHE",finput,default=[1d0,1d0,1d0,1d0,1d0])
  call parse_input_variable(Dbethe,"DBETHE",finput,default=[0d0,0d0,0d0,0d0,0d0])
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
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
  call dmft_write_gf(Weiss,"Weiss",axis='m',iprint=4)

  allocate(Self(Nlso,Nlso,Lmats));Self=zero
  deallocate(Ebands,Dbands,Wband,H0,de)
  allocate(Ebands(Nlso,Le))
  allocate(Dbands(Nlso,Le))
  allocate(Wband(Nlso))
  allocate(H0(Nlso))
  allocate(de(Nlso))
  do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
     io = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
     jo = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
     Self(io,jo,:) = Smats(ilat,ispin,jspin,iorb,jorb,:)
  enddo
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

  call finalize_MPI()

  stop
  
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







  subroutine dmft_kinetic_energy_Bethe(Ebands,Dbands,Hloc,Sigma)
    real(8),dimension(Nlso,Le),intent(in) :: Ebands    ![Nlat*Nspin*Norb][Lk]
    real(8),dimension(Nlso,Le),intent(in) :: Dbands    ![Nlat*Nspin*Norb][Lk]Lk]
    real(8),dimension(Nlso),intent(in)    :: Hloc    ![Nlat*Nspin*Norb]
    complex(8),dimension(Nlso,Nlso,Le)    :: Sigma   ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][L]
    !
    integer                               :: Lk,Liw
    integer                               :: i,ik,iso
    integer                               :: ilat,unit
    real(8),dimension(Nlso,Nlso)          :: Sigma_HF
    !
    complex(8)                            :: Ak,Bk,Ck,Dk
    complex(8)                            :: Gk,Tk
    !
    complex(8),dimension(:,:),allocatable :: Ak_,Bk_,Ck_,Dk_
    complex(8),dimension(:,:),allocatable :: Gk_,Tk_
    !
    real(8),dimension(Nlso)               :: Tail0,Tail1
    real(8),dimension(Nlso)               :: Lail0,Lail1
    real(8)                               :: spin_degeneracy
    !
    real(8),dimension(Nlso)               :: H0,Hl
    real(8),dimension(Nlso)               :: H0tmp,Hltmp
    !
    real(8),dimension(Nlso)               :: Ekin_,Eloc_
    real(8),dimension(Nlat,Nso)           :: Ekin,Eloc
    !
    real(8),dimension(:),allocatable      :: wm
    logical                               :: dos_diag !1. T / 2. F
    integer                               :: mpi_rank,mpi_size
    logical                               :: mpi_master
    !
    !
    !MPI setup:
    mpi_size  = get_size_MPI()
    mpi_rank =  get_rank_MPI()
    mpi_master= get_master_MPI()
    !
    !
    Lk   = size(Ebands,2)
    Liw  = size(Sigma,3)
    !
    !Testing:
    call assert_shape(Sigma,[Nlso,Nlso,Liw],"dmft_kinetic_energy_normal_Nso_dos","Sigma")
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    wm = pi/beta*dble(2*arange(1,Liw)-1)
    !
    !Get HF part of the self-energy
    Sigma_HF = dreal(Sigma(:,:,Liw))
    !
    !
    if(mpi_master)write(*,"(A)") "Kinetic energy computation"
    if(mpi_master)call start_timer()
    H0=0d0
    Hl=0d0
    H0tmp= 0d0
    Hltmp= 0d0
    !
    Tail0=0d0
    Tail1=0d0
    Lail0=0d0
    Lail1=0d0
    !
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk
       do iso=1,Nlso
          Ak = Ebands(iso,ik)
          Bk =-Ebands(iso,ik) - Sigma_HF(iso,iso)
          do i=1+mpi_rank,Liw,mpi_size
             Gk = (xi*wm(i)+xmu) - Sigma(iso,iso,i) - Ebands(iso,ik) - Hloc(iso)
             Gk = 1d0/Gk
             Tk = 1d0/(xi*wm(i)) - Bk/(xi*wm(i))**2
             Ck = Ak*(Gk - Tk)
             Dk = Hloc(iso)*(Gk - Tk)
             H0tmp(iso) = H0tmp(iso) + Dbands(iso,ik)*Ck
             Hltmp(iso) = Hltmp(iso) + Dbands(iso,ik)*Dk
          enddo
       enddo
       if(mpi_master)call eta(ik,Lk)
    enddo
    call AllReduce_MPI(MPI_COMM_WORLD,H0tmp,H0)
    call AllReduce_MPI(MPI_COMM_WORLD,Hltmp,Hl)
    if(mpi_master)call stop_timer()
    spin_degeneracy=3d0-Nspin     !2 if Nspin=1, 1 if Nspin=2
    H0=H0/beta*2*spin_degeneracy
    Hl=Hl/beta*2*spin_degeneracy
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    do ik=1,Lk
       do iso=1,Nlso
          Ak = Ebands(iso,ik)
          Bk =-Ebands(iso,ik) - Sigma_HF(iso,iso)
          Ck= Ak*Bk
          Dk= Hloc(iso)*Bk
          Tail0(iso) = Tail0(iso) + 0.5d0*Dbands(iso,ik)*Ak
          Tail1(iso) = Tail1(iso) + 0.25d0*Dbands(iso,ik)*Ck
          Lail0(iso) = Lail0(iso) + 0.5d0*Dbands(iso,ik)*Hloc(iso)
          Lail1(iso) = Lail1(iso) + 0.25d0*Dbands(iso,ik)*Dk
       enddo
    enddo
    Tail0=Tail0*spin_degeneracy
    Tail1=Tail1*beta*spin_degeneracy
    Lail0=Lail0*spin_degeneracy
    Lail1=Lail1*beta*spin_degeneracy
    !
    Ekin_=H0+Tail0+Tail1
    Eloc_=Hl+Lail0+Lail1
    !
    if(mpi_master)then
       unit = free_unit()
       open(unit,file="dmft_kinetic_energy.info")
       write(unit,"(A1,90(A14,1X))")"#",&
            str(1)//"<K>",str(2)//"<Eloc>",&
            (str(2+i)//"<K"//str(i)//">",i=1,Nso),&
            (str(2+Nso+i)//"<Eloc"//str(i)//">",i=1,Nso)
       close(unit)
       !
       !
       do ilat=1,Nlat
          do iso=1,Nso
             i = iso + (ilat-1)*Nso
             Ekin(ilat,iso) = Ekin_(i)
             Eloc(ilat,iso) = Eloc_(i)
          enddo
       enddo
       !
       unit = free_unit()
       open(unit,file="dmft_kinetic_energy.dat")       
       do ilat=1,Nlat
          write(unit,"(100000F15.9)")sum(Ekin(ilat,:)),sum(Eloc(ilat,:)),&
               (Ekin(ilat,i),i=1,Nso),&
               (Eloc(ilat,i),i=1,Nso)
       enddo
       write(unit,"(90F15.9)")sum(Ekin)/Nlat,sum(Eloc)/Nlat
       close(unit)
    endif
    !
    deallocate(wm)
  end subroutine dmft_kinetic_energy_Bethe



end program

  









