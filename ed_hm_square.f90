program ed_hm_square
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                       :: iloop,Nb,Lk,Nx,Nso,ik,iorb
  logical                                       :: converged
  real(8)                                       :: wband,wmixing
  real(8),dimension(5)                          :: ts,Dband
  real(8),dimension(:),allocatable              :: dens
  !Bath:
  real(8),allocatable                           :: Bath(:),Bath_(:)
  !The local hybridization function:
  complex(8),allocatable                        :: Hloc(:,:,:,:)
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Weiss,Weiss_
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gkmats
  complex(8),allocatable,dimension(:)           :: Gtest
  !
  character(len=16)                             :: finput
  complex(8),allocatable                        :: Hk(:,:,:)
  real(8),allocatable                           :: Wt(:)
  !
  integer                                       :: comm,rank
  logical                                       :: master
  logical                                       :: mixG0,symOrbs

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(wmixing,"wmixing",finput,default=0.5d0,comment="Mixing bath parameter")
  call parse_input_variable(ts,"TS",finput,default=[0.25d0,0d0,0d0,0d0,0d0],comment="hopping parameter")
  call parse_input_variable(Dband,"Dband",finput,default=[0d0,0d0,0d0,0d0,0d0],comment="cystal field splittig (bands shift)")
  call parse_input_variable(Nx,"Nx",finput,default=100,comment="Number of kx point for 2d BZ integration")
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  call parse_input_variable(symOrbs,"symOrbs",finput,default=.false.)
  !
  call ed_read_input(trim(finput),comm)
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  if(Nspin/=1.OR.Norb>5)stop "Wrong setup from input file: Nspin/=1 OR Norb>5"
  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats),Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats),Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  allocate(dens(Norb))
  allocate(Gtest(Lmats))

  !Build Hk
  call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
  Lk = Nx*Nx
  allocate(Hk(Nso,Nso,Lk),Wt(Lk))
  call TB_build_model(Hk(:,:,:),hk_model,Nso,[Nx,Nx])
  Wt = 1d0/Lk
  Hloc   = zero
  Hloc(1,1,:,:) = sum(Hk,dim=3)/Lk
  where(abs(dreal(Hloc))<1.d-6) Hloc=0d0
  
  if(master)call TB_write_hk(Hk(:,:,:),"Hk2d.dat",Nlat=1,&
                             Nspin=1,&
                             Norb=Norb,&
                             Nkvec=[Nx,Nx])

  !setup solver
  Nb=ed_get_bath_dimension()
  allocate(bath(Nb))
  allocate(Bath_(Nb))
  call ed_init_solver(comm,bath)



  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath,Hloc) 
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)
     call ed_get_dens(dens)

     !Compute the local gfs:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     call dmft_gloc_realaxis(Hk,Greal,Sreal)
     !Print the local gfs:
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
     call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)

     !Get the Weiss field/Delta function to be fitted
     call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,SCtype=cg_scheme)
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)

     if(mixG0)then
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
        Weiss_=Weiss
     endif
     
     !Perform the SELF-CONSISTENCY by fitting the new bath
     if(symOrbs)then
        call ed_chi2_fitgf(comm,Weiss,bath,ispin=1,iorb=1)
        call ed_orb_equality_bath(bath,save=.true.)
     else
        call ed_chi2_fitgf(comm,Weiss,bath,ispin=1)
     endif

     !MIXING:
     if(.not.mixG0)then
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
        Bath_=Bath
     endif
     
     !Check convergence (if required change chemical potential)     
     Gtest=zero
     do iorb=1,Norb
        Gtest=Gtest+Weiss(1,1,iorb,iorb,:)/Norb
     enddo
     converged = check_convergence(Gtest,dmft_error,nsuccess,nloop,reset=.false.)
     if(nread/=0d0)call ed_search_variable(xmu,sum(dens),converged)

     call end_loop
  enddo

  !Get kinetic energy:
  call dmft_kinetic_energy(Hk,Smats)

!  allocate(Gkmats(Lk,Nspin,Nspin,Norb,Norb,Lmats))
!  do ik=1,Lk
!     call dmft_gk_matsubara(Hk(:,:,ik),Gkmats(ik,:,:,:,:,:),Smats)
!  enddo

  call finalize_MPI()

contains

  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Hk model for the 2d square lattice
  !-------------------------------------------------------------------------------------------
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:) :: kpoint
    integer              :: N,ih
    real(8)              :: kx,ky
    complex(8)           :: hk(N,N)
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = zero
    do ih=1,N
       Hk(ih,ih) = -one*2d0*ts(ih)*(cos(kx)+cos(ky)) + Dband(ih)
    enddo
  end function hk_model


end program ed_hm_square



