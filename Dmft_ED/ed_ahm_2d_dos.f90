program ed_ahm_2d
  USE EDIPACK
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                   :: Le
  real(8),parameter                         :: D=1
  integer                                   :: iloop,i
  integer                                   :: Nb,Nso
  logical                                   :: converged
  real(8)                                   ::dens_average
  real(8)                                   ::dens_ed,mua,mub
  !Bath:
  real(8),allocatable,dimension(:)          :: Bath,Bath_Prev
  !The local functions:
  complex(8),allocatable,dimension(:,:,:,:) :: Gmats
  complex(8),allocatable,dimension(:,:,:,:) :: Greal
  complex(8),allocatable,dimension(:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:) :: Weiss,Peiss
  complex(8),allocatable,dimension(:,:)     :: reduced_density_matrix
  !
  character(len=16)                         :: finput
  real(8)                                   :: wmixing
  real(8)                                   :: de
  real(8),dimension(:,:),allocatable        :: Ddos
  real(8),dimension(:,:,:),allocatable      :: Edos
  real(8),dimension(:,:),allocatable        :: H0
  real(8),dimension(:,:),allocatable        ::lambdasym_vector
  complex(8),dimension(:,:,:),allocatable   :: Hsym_basis
  integer                                   :: nkpath
  logical                                   :: mix_bath
  complex(8),dimension(:,:),allocatable     :: rdm
  !MPI
  integer                                       :: comm,rank
  logical                                       :: master


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)



  call parse_cmd_variable(finput,"FINPUT",default='inputAHM.conf')
  call parse_input_variable(Le,"LE",finput,default=1000)
  call parse_input_variable(mua,"mua",finput,default=0d0)
  call parse_input_variable(mub,"mub",finput,default=0d0)
  call parse_input_variable(mix_bath,"mix_bath",finput,default=.true.)
  call parse_input_variable(wmixing,"wmixing",finput,default=0.5d0,comment="Mixing bath parameter")
  !
  call ed_read_input(trim(finput))


  !Add DMFT CTRL Variables used in DMFT_TOOLS:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  Nso = Nspin*Norb;if(Nso/=1)stop "Nso != 1"


  !Allocate local functions (Nambu)
  allocate(Gmats(2,Nso,Nso,Lmats))
  allocate(Smats(2,Nso,Nso,Lmats))
  allocate(Weiss(2,Nso,Nso,Lmats))
  allocate(Greal(2,Nso,Nso,Lreal))
  allocate(Sreal(2,Nso,Nso,Lreal))
  allocate(Peiss(2,Nso,Nso,Lmats))




  !> Allocate the 2d DOS and dispersion for the Nambu structure in DMFT_TOOLS
  !  we ask for separate dispersions (or two H(k)) for the two elements on the diagonal
  allocate(Edos(2,Nso,Le))  
  Edos(1,1,:) = linspace(-D,D,Le,mesh=de)
  Edos(2,1,:) =-linspace(-D,D,Le,mesh=de)
  allocate(Ddos(Nso,Le))
  do i=1,Le
     Ddos(1,i) = dens_2dsquare(Edos(1,1,i),D)
  enddo
  Ddos= Ddos/simps(Ddos(1,:),-D,D)*de


  !> Get local Hamiltonian (used in DMFT_TOOLS)
  allocate(H0(2,1))
  H0=0d0
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_prev(Nb))


  !set Hloc
  call ed_set_Hloc(zeros(Nso,Nso))

  !> Initialize the ED solver (bath is guessed or read from file) 
  call ed_init_solver(bath)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !> Solve the impurity problem, retrieve matsubara self-energy 
     if(nread/=0.d0)then
        if(master)call fzero(get_dens,mua,mub,i,tol_rel=1d-4,tol_abs=1d-4)
        call Bcast_MPI(MPI_COMM_WORLD,mua)
        xmu = mua
        mua = xmu-1d-1
        mub = xmu+1d-1
     endif
     !
     call ed_solve(bath)     
     !> Retrieve impurity self-energies (normal, anomalous)
     call ed_get_sigma(Smats(1,:,:,:),axis='m',type='n')
     call ed_get_sigma(Smats(2,:,:,:),axis='m',type='a')


     !> Get Gloc (using DMFT_TOOLS)
     call get_gloc(Edos,Ddos,H0,Gmats,Smats,axis='m')
     call write_gf(Gmats(1,:,:,:),"Gloc",axis='mats',iprint=1)
     call write_gf(Gmats(2,:,:,:),"Floc",axis='mats',iprint=1)

     !> Update the Weiss field: (using DMFT_TOOLS).
     call dmft_self_consistency(&
          Gmats(1,:,:,:),Gmats(2,:,:,:),&
          Smats(1,:,:,:),Smats(2,:,:,:),&
          Weiss(1,:,:,:),Weiss(2,:,:,:) )
     call write_gf(Weiss(1,:,:,:),"Weiss",axis='mats',iprint=1)
     call write_gf(Weiss(2,:,:,:),"Theta",axis='mats',iprint=1)

     if((.not.mix_bath).AND.iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Peiss
     Peiss=Weiss

     !> Fit the new bath, starting from the current bath + the update Weiss field
     call ed_chi2_fitgf(Weiss(1,:,:,:),Weiss(2,:,:,:),bath,ispin=1)
     !
     !
     !> Linear mixing of the bath (this can be done in alternative of mixing Weiss)
     if(mix_bath.AND.iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath

     !Check convergence (using DMFT_TOOLS)
     converged = check_convergence(Weiss(1,1,1,:),dmft_error,nsuccess,nloop)

     call end_loop
  enddo

  !> retrieve real-axis self-energies (n, a) and build local Green's function
  call ed_get_sigma(Sreal(1,:,:,:),axis='r',type='n')
  call ed_get_sigma(Sreal(2,:,:,:),axis='r',type='a')
  call get_gloc(Edos,Ddos,H0,Greal,Sreal,axis='r')
  call write_gf(Greal(1,:,:,:),"Gloc",axis='real',iprint=1)
  call write_gf(Greal(2,:,:,:),"Floc",axis='real',iprint=1)



  !Compute the Kinetic Energy:
  call dmft_kinetic_energy(Edos(1,:,:),Ddos,H0(1,:),Smats(1,:,:,:),Smats(2,:,:,:))


  call Finalize_MPI()
  
contains

  function get_dens(mu) result(dn)
    real(8),intent(in) :: mu
    real(8)            :: dn
    xmu = mu
    call ed_solve(bath,flag_gf=.false.)
    call ed_get_dens(dens_ed)
    dn = dens_ed-nread
    write(LOGfile,*)"ED-density :", mu,dn,dens_ed
  end function get_dens


end program ed_ahm_2d



