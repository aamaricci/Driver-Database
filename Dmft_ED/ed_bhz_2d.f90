program ed_bhz
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE SF_MPI
  implicit none
  integer                                     :: iloop,Lk,Nso,Lakw
  logical                                     :: converged
  !Bath:
  integer                                     :: Nb,iorb,jorb,ispin,print_mode
  real(8),allocatable                         :: Bath(:),Bath_(:)
  !The local hybridization function:
  complex(8),allocatable                      :: Weiss(:,:,:,:,:),Weiss_(:,:,:,:,:)
  complex(8),allocatable                      :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable                      :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  complex(8),allocatable,dimension(:)         :: Gtest
  !hamiltonian input:
  complex(8),allocatable                      :: Hk(:,:,:),bhzHloc(:,:),sigmaBHZ(:,:),Zmats(:,:),impRho(:,:)
  real(8),allocatable                         :: Wtk(:)
  real(8),allocatable                         :: kxgrid(:),kygrid(:)
  integer,allocatable                         :: ik2ix(:),ik2iy(:)
  !variables for the model:
  integer                                     :: Nk,Nkpath
  real(8)                                     :: mh,lambda,wmixing,akrange
  character(len=30)                           :: Params
  character(len=16)                           :: finput
  character(len=32)                           :: hkfile
  logical                                     :: spinsym,usez,mixG0
  !
  real(8),dimension(2)                        :: Eout
  real(8),allocatable                         :: dens(:)
  complex(8),dimension(4,4)                   :: Gamma1,Gamma2,Gamma5,GammaN,GammaS,GammaJ
  complex(8),dimension(4,4)                   :: GammaE0,GammaEx,GammaEy,GammaEz
  complex(8),dimension(4,4)                   :: GammaR0,GammaRx,GammaRy,GammaRz
  real(8),dimension(:),allocatable            :: lambdasym_vector
  complex(8),dimension(:,:,:,:,:),allocatable :: Hsym_basis
  !MPI Vars:
  integer                                     :: irank,comm,rank,ierr
  logical                                     :: master,getbands,getakw,getXem,getPave

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.in')  
  call parse_input_variable(Params,"Params",finput,default="E0EzEx",&
       comment="Ex; EzEx; E0Ex; ExEy; E0Ez; E0EzEx; E0EzExEy")
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
  call parse_input_variable(usez,"USEZ",finput,default=.false.)
  call parse_input_variable(getXem,"GETXEM",finput,default=.false.)
  call parse_input_variable(getPave,"getPave",finput,default=.false.)
  call parse_input_variable(getbands,"GETBANDS",finput,default=.false.)
  call parse_input_variable(getakw,"GETAKW",finput,default=.false.)
  call parse_input_variable(Lakw,"LAKW",finput,default=250)
  call parse_input_variable(akrange,"AKRANGE",finput,default=5d0)
  !
  call ed_read_input(trim(finput))
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gtest(Lmats))
  allocate(dens(Norb))
  allocate(SigmaBHZ(Nso,Nso))
  allocate(Zmats(Nso,Nso))

  gamma1=kron_pauli( pauli_sigma_z, pauli_tau_x)
  gamma2=kron_pauli( pauli_sigma_0,-pauli_tau_y)
  gamma5=kron_pauli( pauli_sigma_0, pauli_tau_z)
  gammaN=kron_pauli( pauli_sigma_0, pauli_tau_0)
  ! gammaJ=kron_pauli( pauli_sigma_x, (pauli_tau_0+pauli_tau_3)/2d0)
  ! gammaJ=kron_pauli( pauli_sigma_x, pauli_tau_0) + kron_pauli( pauli_sigma_0, pauli_tau_x)
  gammaJ=kron_pauli( pauli_sigma_x, pauli_tau_0)
  !
  gammaE0=kron_pauli( pauli_sigma_0, pauli_tau_x )
  gammaEx=kron_pauli( pauli_sigma_x, pauli_tau_x )
  gammaEy=kron_pauli( pauli_sigma_y, pauli_tau_x )
  gammaEz=kron_pauli( pauli_sigma_z, pauli_tau_x )
  !
  gammaR0=kron_pauli( pauli_sigma_0, pauli_tau_y )
  gammaRx=kron_pauli( pauli_sigma_x, pauli_tau_y )
  gammaRy=kron_pauli( pauli_sigma_y, pauli_tau_y )
  gammaRz=kron_pauli( pauli_sigma_z, pauli_tau_y )



  if(getPave)then
     call get_Pave()
     call finalize_MPI()
     stop
  endif



  if(getXem)then
     call get_Xem(Lakw)
     call finalize_MPI()
     stop
  endif




  !Buil the Hamiltonian on a grid or on  path
  call set_sigmaBHZ()
  call build_hk(trim(hkfile))

  print_mode=1
  if(ed_mode=="nonsu2")print_mode=4

  if(getbands)then
     call read_array("Smats",Smats)
     call solve_hk_topological(so2j(Smats(:,:,:,:,1)))
     call finalize_MPI()
     stop
  endif


  if(getakw)then
     call read_array("Sreal",Sreal)
     call get_Akw(Lakw,Akrange)
     call finalize_MPI()
     stop
  endif




  !Setup solver
  if(bath_type=="replica")then
     select case(trim(Params))
     case default
        stop "Params not in [Ex; EzEx; E0Ex; ExEy; E0Ez; E0EzEx; E0EzExEy]"
     case("Ex")
        allocate(lambdasym_vector(2))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,2))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaEx) ;lambdasym_vector(2)=-sb_field
        call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
        Nb=ed_get_bath_dimension(2)!(Hsym_basis)

     case("EzEx")
        allocate(lambdasym_vector(3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaEz) ;lambdasym_vector(2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEx) ;lambdasym_vector(3)=-sb_field
        call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
        Nb=ed_get_bath_dimension(3)!(Hsym_basis)

     case("E0Ex")
        allocate(lambdasym_vector(3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0) ;lambdasym_vector(2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEx) ;lambdasym_vector(3)=-sb_field
        call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
        Nb=ed_get_bath_dimension(3)!(Hsym_basis)

     case("ExEy")
        allocate(lambdasym_vector(3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaEx) ;lambdasym_vector(2)=-sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEy) ;lambdasym_vector(3)=sb_field
        call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
        Nb=ed_get_bath_dimension(3)!(Hsym_basis)
        
     case("E0Ez")
        allocate(lambdasym_vector(3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0) ;lambdasym_vector(2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEz) ;lambdasym_vector(3)=sb_field
        call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
        Nb=ed_get_bath_dimension(3)!(Hsym_basis)
        
     case("E0EzEx")
        allocate(lambdasym_vector(4))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,4))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0) ;lambdasym_vector(2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEz) ;lambdasym_vector(3)=sb_field
        Hsym_basis(:,:,:,:,4)=j2so(GammaEx) ;lambdasym_vector(4)=-sb_field
        call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
        Nb=ed_get_bath_dimension(4)!(Hsym_basis)
        
     case("E0EzExEy")
        allocate(lambdasym_vector(5))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,5))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0) ;lambdasym_vector(2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEz) ;lambdasym_vector(3)=sb_field
        Hsym_basis(:,:,:,:,4)=j2so(GammaEx) ;lambdasym_vector(4)=-sb_field
        Hsym_basis(:,:,:,:,5)=j2so(GammaEy) ;lambdasym_vector(5)=sb_field
        call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
        Nb=ed_get_bath_dimension(5)!(Hsym_basis)
        
     case("E0EzExR0RzRx")
        allocate(lambdasym_vector(7))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,7))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0) ;lambdasym_vector(2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEz) ;lambdasym_vector(3)=sb_field
        Hsym_basis(:,:,:,:,4)=j2so(GammaEx) ;lambdasym_vector(4)=-sb_field
        Hsym_basis(:,:,:,:,5)=j2so(GammaR0) ;lambdasym_vector(5)=sb_field
        Hsym_basis(:,:,:,:,6)=j2so(GammaRz) ;lambdasym_vector(6)=sb_field
        Hsym_basis(:,:,:,:,7)=j2so(GammaRx) ;lambdasym_vector(7)=-sb_field
        call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
        Nb=ed_get_bath_dimension(7)!(Hsym_basis)
        
     end select


  else     
     Nb=ed_get_bath_dimension()
  endif


  call ed_set_hloc(j2so(bhzHloc))

  print*,Nb
  allocate(Bath(Nb))
  allocate(Bath_(Nb))
  call ed_init_solver(bath)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)
     call ed_get_sigma(Smats,axis='mats')
     call ed_get_sigma(Sreal,axis='real')
     call ed_get_dens(dens)

     !Get GLOC:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     call dmft_write_gf(Gmats,"Gloc",axis='mats',iprint=print_mode)


     !Update WeissField:
     if(cg_scheme=='delta')then        
        call dmft_self_consistency(Gmats,Smats,Weiss,j2so(bhzHloc))
     else
        call dmft_self_consistency(Gmats,Smats,Weiss)
     endif
     call dmft_write_gf(Weiss,"Weiss",axis='mats',iprint=print_mode)


     if(mixG0)then
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
        Weiss_=Weiss
     endif

     !Fit the new bath, starting from the old bath + the supplied Weiss
     select case(ed_mode)
     case default
        stop "ed_mode!=Normal/Nonsu2"
     case("normal")
        call ed_chi2_fitgf(Weiss,bath,ispin=1)
        if(.not.spinsym)then
           call ed_chi2_fitgf(Weiss,bath,ispin=2)
        else
           call ed_spin_symmetrize_bath(bath,save=.true.)
        endif
     case("nonsu2")
        call ed_chi2_fitgf(Weiss,bath)
     end select

     if(.not.mixG0)then
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
        Bath_=Bath
     endif
     !

     !Check convergence (if required change chemical potential)
     ! Gtest=zero
     ! do ispin=1,Nspin
     !    do iorb=1,Norb
     !       Gtest=Gtest+Weiss(ispin,ispin,iorb,iorb,:)/Norb/Nspin
     !    enddo
     ! enddo
     Gtest=Weiss(1,1,1,1,:)
     converged = check_convergence(Gtest,dmft_error,nsuccess,nloop)
     if(nread/=0d0)call ed_search_variable(xmu,sum(dens),converged)

     call end_loop
  enddo


  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  call dmft_write_gf(Greal,"Gloc",axis='real',iprint=print_mode)

  call dmft_kinetic_energy(Hk,Smats)

  call solve_hk_topological(so2j(Smats(:,:,:,:,1)))

  call save_array("Smats",Smats)
  call save_array("Sreal",Sreal)

  call finalize_MPI()



contains


  !---------------------------------------------------------------------
  !PURPOSE: GET BHZ HAMILTONIAN (from the NonInteracting code)
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j,ik=0
    integer                             :: ix,iy
    real(8)                             :: kx,ky    
    integer                             :: iorb,jorb
    integer                             :: isporb,jsporb
    integer                             :: ispin,jspin
    real(8)                             :: foo
    integer                             :: unit
    complex(8),dimension(Nso,Nso,Lmats) :: Gmats
    complex(8),dimension(Nso,Nso,Lreal) :: Greal
    real(8)                             :: wm(Lmats),wr(Lreal),dw
    !
    call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
    !
    if(master)write(LOGfile,*)"Build H(k) for BHZ:"
    Lk=Nk**2
    if(master)write(*,*)"# of k-points     :",Lk
    if(master)write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk)) ;Hk=zero
    !
    call TB_build_model(Hk,hk_bhz,Nso,[Nk,Nk])
    if(master.AND.present(file))then
       call TB_write_hk(Hk,trim(file),&
            Nlat=1,&
            Nspin=1,&
            Norb=Norb,&
            Nkvec=[Nk,Nk])
    endif
    allocate(bhzHloc(Nso,Nso))
    bhzHloc = zero
    bhzHloc = sum(Hk,dim=3)/Lk
    where(abs(dreal(bhzHloc))<1d-6)bhzHloc=zero
    if(master)  call TB_write_Hloc(bhzHloc)
  end subroutine build_hk






  !--------------------------------------------------------------------!
  !PURPOSE: Set the Self-Energy
  !--------------------------------------------------------------------!
  subroutine set_SigmaBHZ(sigma)
    complex(8),dimension(Nso,Nso),optional :: sigma
    sigmaBHZ = zero;if(present(sigma))sigmaBHZ=sigma
  end subroutine set_SigmaBHZ


  !--------------------------------------------------------------------!
  !PURPOSE: Solve the topological Hamiltonian
  !--------------------------------------------------------------------!
  subroutine solve_hk_topological(sigma)
    integer                                :: i,j
    integer                                :: Npts
    complex(8),dimension(Nso,Nso)          :: sigma(Nso,Nso)
    real(8),dimension(:,:),allocatable     :: kpath
    !
    if(master)then
       !This routine build the H(k) along the GXMG path in BZ, Hk(k) is constructed along this path.
       write(LOGfile,*)"Build H_TOP(k) BHZ along path:"
       !
       call set_sigmaBHZ()
       !
       Npts = 4
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_gamma
       kpath(2,:)=kpoint_x1
       kpath(3,:)=kpoint_m1
       kpath(4,:)=kpoint_gamma
       call set_sigmaBHZ(sigma)
       call TB_solve_model(hk_bhz,Nso,kpath,Nkpath,&
            colors_name=[red,blue,red,blue],&
            points_name=[character(len=20) :: "{\Symbol G}","X","M","{\Symbol G}"],&
            file="Eig_Htop.ed")
    endif
  end subroutine solve_hk_topological




  !--------------------------------------------------------------------!
  !BHZ HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_bhz(kvec,N) result(hk)
    integer                   :: N
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: ek,kx,ky
    integer                   :: ii
    if(N/=Nso)stop "hk_bhz error: N != Nspin*Norb == 4"
    kx=kvec(1)
    ky=kvec(2)
    ek = -1d0*(cos(kx)+cos(ky))
    Hk = (Mh+ek)*Gamma5 + lambda*sin(kx)*Gamma1 + lambda*sin(ky)*Gamma2
    ! !add the SigmaBHZ term to get Topologial Hamiltonian if required:
    Hk = Hk + dreal(SigmaBHZ)
    if (usez) then
       Zmats=zero
       do ii=1,Nso
          Zmats(ii,ii)  = 1.d0/abs( 1.d0 +  abs(dimag(sigmaBHZ(ii,ii))/(pi/beta)) )
       end do
       Hk = matmul(Zmats,Hk)
    endif
  end function hk_bhz






  !---------------------------------------------------------------------
  !PURPOSE: GET A(k,w)
  !---------------------------------------------------------------------
  subroutine get_Akw(Lw,aw)
    integer                                       :: Lw
    real(8)                                       :: Aw
    integer                                       :: ik
    integer                                       :: iorb,jorb
    integer                                       :: ispin,jspin
    complex(8),dimension(:,:,:,:,:),allocatable   :: Sreal_
    real(8),dimension(:,:),allocatable            :: Akreal
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Gkreal
    real(8),allocatable                           :: wr_(:),wr(:),Kgrid(:)
    real(8),dimension(:,:),allocatable            :: Kpath
    character(len=30)                             :: suffix
    !
    !
    !
    call set_SigmaBHZ()
    !
    if(master)print*,"Build A(k,w) using Sigma(w) interpolated"
    allocate(wr(Lreal),wr_(Lw))
    wr  = linspace(wini,wfin,Lreal)!In
    wr_ = linspace(-aw,aw,Lw)      !Out
    !
    wini=-aw
    wfin= aw
    call add_ctrl_var(wini,'wini')
    call add_ctrl_var(wfin,'wfin')
    !
    allocate(kpath(4,3))
    Lk=(size(kpath,1)-1)*Nkpath
    kpath(1,:)=kpoint_gamma
    kpath(2,:)=kpoint_x1
    kpath(3,:)=kpoint_m1
    kpath(4,:)=kpoint_gamma

    if(allocated(Hk))deallocate(Hk)

    allocate(Hk(Nso,Nso,Lk));Hk=zero
    allocate(kgrid(Lk) )
    !
    call TB_build_kgrid(kpath,Nkpath,kgrid)
    call TB_build_model(Hk,hk_bhz,Nso,kpath,Nkpath)
    !
    !

    allocate(Sreal_(Nspin,Nspin,Norb,Norb,Lw));Sreal_=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                call cubic_spline(wr,Sreal(ispin,jspin,iorb,jorb,:),wr_,Sreal_(ispin,jspin,iorb,jorb,:))
             enddo
          enddo
       enddo
    enddo
    call dmft_write_gf(Sreal_,"Sreal_",axis='real',iprint=print_mode)
    !
    allocate(Gkreal(Lk,Nspin,Nspin,Norb,Norb,Lw));Gkreal=zero        
    call start_timer
    do ik=1,Lk
       call dmft_gk_realaxis(Hk(:,:,ik),Gkreal(ik,:,:,:,:,:),Sreal_)
       call eta(ik,Lk)
    enddo
    call stop_timer
    !

    allocate(Akreal(Lk,Lw));Akreal=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          Akreal = Akreal - dimag(Gkreal(:,ispin,ispin,iorb,iorb,:))/pi/Nso
       enddo
    enddo

    call splot3d("Akw_real_nso.dat",kgrid,wr_,Akreal) 

  end subroutine get_Akw








  !---------------------------------------------------------------------
  !PURPOSE: GET X_em(q-->0,v)
  !---------------------------------------------------------------------
  subroutine get_Xem(Liw)
    integer                                       :: i,ik,iv,iw,Liw,Lkpath,Lk,iq
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Gkmats
    real(8),dimension(:,:),allocatable            :: Kgrid,Qgrid
    real(8),dimension(:,:),allocatable            :: Kpath
    real(8)                                       :: vm,len
    real(8),dimension(2)                          :: qvec,kvec
    complex(8)                                    :: Xem
    complex(8),dimension(Nso,Nso)                 :: Jk,Gkw,Gkqwv,A,B5,Bj,X5,Xj
    complex(8),dimension(4,4) :: GammaJ


    gammaJ=kron_pauli( pauli_sigma_0, pauli_tau_x )


    if(Liw >= Lmats) stop "get_Xem ERROR: Liw > Lmats"

    call set_SigmaBHZ()
    !

    !
    if(master)print*,"Build Xemq(q-->0,v-->0) using Sigma(iw)"
    !
    !
    call build_hk()
    call read_sigma(Smats)

    Lk = Nk**2
    allocate(Kgrid(Lk,2))
    call TB_build_kgrid([Nk,Nk],Kgrid)

    allocate(kpath(2,2))
    Lkpath    = Nkpath
    kpath(1,:)= [0d0,0d0]
    kpath(2,:)= [pi,0d0]
    ! kpath(3,:)= [pi,pi]
    ! kpath(4,:)= [0d0,0d0]
    allocate(Qgrid(Lkpath,2) )
    call TB_build_kgrid(kpath,Nkpath,Qgrid)

    call start_timer()
    do iq=1,Lkpath;iv=1
       qvec = Qgrid(iq,:)
       vm   = pi/beta*2*iv
       Xj   = zero
       do ik=1,Lk
          kvec  = Kgrid(ik,:)
          jk    = hk_dkx_bhz(kvec,Nso)
          do iw=1,Liw-iv
             Gkw   = get_simplified_gf(kvec,iw,Nso)
             Gkqwv = get_simplified_gf(kvec+qvec,iw+iv,Nso)
             A     = matmul(Gkw,jk)
             Bj    = matmul(Gkqwv,GammaJ)
             Xj    = Xj + matmul(A,Bj)
          enddo
       enddo
       Xem = trace(Xj)
       Xem = Xem/beta/(vm+1d-6)/Lk
       write(100,*)qvec(1),dreal(Xem),dimag(Xem)
       call eta(iq,Lkpath)
    enddo
    call stop_timer()
  end subroutine get_Xem








  !---------------------------------------------------------------------
  !PURPOSE: GET <P>=\sum_{k} \sum_{a b} p_{a b}(k) <\cc_{k a} \ca_{k b}>
  !---------------------------------------------------------------------
  subroutine get_Pave()
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: Gkmats
    complex(8),dimension(Nso,Nso,Lmats)               :: Gk
    real(8),dimension(:,:),allocatable                :: Kgrid
    complex(8)                                        :: P
    integer                                           :: i,ik,io,jo
    real(8),dimension(2)                              :: kvec
    complex(8),dimension(Nso,Nso)                     :: pk,rhok,nkk
    real(8),dimension(Nso) :: N

    call set_SigmaBHZ()
    !
    if(master)print*,"Build <P> using Sigma(iw)"
    !
    !

    call read_sigma(Smats)
    call build_hk()
    Lk = Nk**2
    allocate(Kgrid(Lk,2))
    call TB_build_kgrid([Nk,Nk],Kgrid)

    pk = zero
    nkk= zero
    do ik=1,Lk
       kvec = Kgrid(ik,:)
       call dmft_gk_matsubara(Hk(:,:,ik),Gkmats,Smats)
       Gk = so2j_l(Gkmats)
       !
       do io=1,Nso
          do jo=1,Nso
             rhok(io,jo) = fft_get_density(Gk(io,jo,:),beta)
          enddo
       enddo
       pk = pk + hk_dkx_bhz(kvec,Nso)*rhok/Lk
       nkk= nkk+ rhok/Lk
    enddo
    p = sum(pk)
    n = diagonal(nkk)*2.d0
    print*,p
    print*,n
    print*,dreal(pk(1,Nso))
    open(519,file="pave.dat")
    write(519,*)dreal(p),dimag(p)
    write(519,*)n
    write(519,*)dreal(pk(1,Nso))
    do io=1,Nso
       write(519,*)(pk(io,jo),jo=1,Nso)
    enddo
  end subroutine get_Pave









  !--------------------------------------------------------------------!
  !BHZ D_kx H(kx,ky):
  !--------------------------------------------------------------------!
  function hk_dkx_bhz(kvec,N) result(hk)
    integer                   :: N
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk,uk,dhk,wmn,lhk
    real(8)                   :: ek,kx,ky,Mk,xk,yk,eigk(N)
    integer                   :: i,j,ii
    if(N/=Nso)stop "hk_bhz error: N != Nspin*Norb == 4"
    uk  = hk_bhz(kvec,N)
    call eigh(uk,eigk)
    kx  = kvec(1)
    ky  = kvec(2)
    Mk  = Mh-1d0*(cos(kx)+cos(ky))
    xk  = lambda*sin(kx)
    yk  = lambda*sin(ky)
    Ek  = sqrt( Mk**2 + xk**2 + yk**2 )
    wmn = 2*Ek*kron_pauli(0.5d0*(pauli_sigma_0+pauli_sigma_x),pauli_tau_y)
    wmn = matmul(matmul(conjg(transpose(uk)),wmn),uk)
    dHk = sin(kx)*Gamma5 + lambda*cos(kx)*Gamma1
    lhk = wmn*kron_pauli(pauli_sigma_x,pauli_tau_y)!sigma_x \ocirc \tau_y
    Hk  = dHk + lHk   
  end function hk_dkx_bhz


  function get_simplified_gf(kvec,i,N) result(gk)
    real(8),dimension(:)      :: kvec
    complex(8)                :: z
    integer                   :: i,N
    complex(8),dimension(N,N) :: gk,sigma
    real(8)                   :: kx,ky
    real(8)                   :: w_,M_,x_,y,Dx,Ek
    sigma = so2j(Smats(:,:,:,:,i))
    !
    kx = kvec(1)
    ky = kvec(2)
    !
    w_ = pi/beta*(2*i-1) - dimag(sigma(1,1))
    M_ = Mh-1d0*(cos(kx)+cos(ky)) + dreal(sigma(1,1))
    x_ = lambda*sin(kx) + dreal(sigma(1,2))
    y  = lambda*sin(ky)
    Dx = dreal(sigma(1,4))
    !
    Ek = M_**2 + x_**2 + y**2 + Dx**2
    Gk = xi*w_*eye(Nso) + M_*Gamma5 + x_*Gamma1  + y*Gamma2 + Dx*GammaEx
    Gk = Gk/(-w_**2-Ek)
  end function get_simplified_gf




  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
  !--------------------------------------------------------------------!
  subroutine read_sigma(sigma)
    complex(8)        :: sigma(:,:,:,:,:)
    integer           :: iorb,jorb,ispin,jspin,i,L,unit
    real(8)           :: reS,imS,ww
    character(len=20) :: suffix
    logical :: bool
    if(size(sigma,1)/=Nspin)stop "read_sigma: error in dim 1. Nspin"
    if(size(sigma,3)/=Norb)stop "read_sigma: error in dim 3. Norb"
    L=size(sigma,5)
    if(L/=Lmats.AND.L/=Lreal)stop "read_sigma: error in dim 5. Lmats/Lreal"
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
                if(L==Lreal)&
                     suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_realw.ed"
                !
                write(*,*)"read from file=","impSigma"//reg(suffix)
                inquire(file="impSigma"//reg(suffix),exist=bool)
                if(bool)then
                   unit=free_unit()
                   open(unit,file="impSigma"//reg(suffix),status='old')
                   do i=1,L
                      read(unit,"(F26.15,6(F26.15))")ww,imS,reS
                      sigma(ispin,jspin,iorb,jorb,i)=dcmplx(reS,imS)
                   enddo
                   close(unit)
                else
                   sigma(ispin,jspin,iorb,jorb,:)=zero
                endif

             enddo
          enddo
       enddo
    enddo
  end subroutine read_sigma


  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop "error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index


  function so2j(fg) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function so2j

  function j2so(fg) result(g)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2so




  function so2j_l(fg) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: fg
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats) :: g
    integer                                     :: iw,i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                do iw=1,Lmats
                   g(i,j,iw) = fg(ispin,jspin,iorb,jorb,iw)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function so2j_l

  function j2so_l(fg) result(g)
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats) :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: g
    integer                                     :: w,i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                do w=1,Lmats
                   g(ispin,jspin,iorb,jorb,w)  = fg(i,j,w)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function j2so_l

end program ed_bhz
