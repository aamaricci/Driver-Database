program ed_ahm_square
  USE EDIPACK
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                       :: iloop,Lk,Nso,Nambu,Nno
  logical                                       :: converged
  !Bath:
  integer                                       :: Nb,Nsym
  real(8),allocatable                           :: Bath(:),Prev(:)
  real(8),dimension(:,:),allocatable            :: LambdaVec
  complex(8),dimension(:,:,:),allocatable       :: Hsym_basis
  !variables for the model:
  integer                                       :: Nx,Nkpath
  real(8)                                       :: delta,lambda,wmixing,ts,alpha,D
  character(len=16)                             :: finput
  logical                                       :: phsym,normal_bath
  !The local hybridization function:
  complex(8),allocatable                        :: Hloc(:,:)
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Greal,Smats,Sreal,Weiss
  complex(8),allocatable                        :: Hk(:,:,:,:)
  !MPI
  integer                                       :: comm,rank
  logical                                       :: master


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(Nsym,"NSYM",finput,default=2,comment="2: diagonal SC; 3: diag+off-diag SC")
  call parse_input_variable(ts,"TS",finput,default=0.5d0,comment="hopping parameter")
  call parse_input_variable(alpha,"ALPHA",finput,default=1d0,comment="bandwidth ratio t_2 = alpha*t_1")
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0,comment="inter-orbital hopping if Norb>1")
  call parse_input_variable(delta,"DELTA",finput,default=0.5d0,comment="crystal field splitting if Norb>1")
  call parse_input_variable(Nx,"Nx",finput,default=10,comment="Number of kx point for 2d BZ integration")
  call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
  call parse_input_variable(phsym,"phsym",finput,default=.false.,comment="Flag to enforce p-h symmetry of the bath.")
  call parse_input_variable(normal_bath,"normal",finput,default=.false.,comment="Flag to enforce no symmetry braking in the bath.")
  !
  call ed_read_input(trim(finput))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  if(Nspin/=1.OR.Norb>2)stop "This code is intended as a driver for the Norb<=2 and Nspin=1 problem"
  Nambu = 2
  Nso   = Nspin*Norb
  Nno   = Nambu*Nso


  D = 4*ts

  !Allocate Dynamical Fields:
  allocate(Gmats(2,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(2,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(2,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(2,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Weiss(2,Nspin,Nspin,Norb,Norb,Lmats))



  !> Set the basis vectors square lattice
  call TB_set_ei([1d0,0d0],[0d0,1d0])  ! real-space lattice basis vectors
  call TB_set_bk([pi2,0d0],[0d0,pi2])   ! k-space lattice basis vectors
  !
  Lk = Nx*Nx
  allocate(Hk(2,Nso,Nso,Lk))    !Need two copies to comply with Nambu basis (11,22)
  call TB_build_model(Hk(1,:,:,:),hk_model,Nso,[Nx,Nx])
  Hk(2,:,:,:) = -conjg(Hk(1,:,:,:)) !transpose?
  !
  allocate(Hloc(Nso,Nso));Hloc=zero
  Hloc = sum(Hk(1,:,:,:),dim=3)/Lk
  where(abs(dreal(Hloc))<1d-6)Hloc=zero


  !Build this any, u never know
  select case(Norb)
  case (1)
     Nsym = 2
     allocate(Hsym_basis(Nno,Nno,Nsym))
     allocate(LambdaVec(Nbath,Nsym))
     !Hsym(:,:) = pauli_Nambu
     Hsym_basis(:,:,1)=pauli_tau_z ; LambdaVec(:,1)= -D+2*D/(Nbath-1)*(arange(1,Nbath)-1)
     Hsym_basis(:,:,2)=pauli_tau_x ; LambdaVec(:,2)= sb_field
  case (2)
     allocate(Hsym_basis(Nno,Nno,Nsym))
     allocate(LambdaVec(Nbath,Nsym))
     !Hsym(:,:) = kron(pauli_Nambu,pauli_Orb)
     Hsym_basis(:,:,1)=kron( pauli_sigma_z, pauli_tau_0) ; LambdaVec(:,1)= -D+2*D/(Nbath-1)*(arange(1,Nbath)-1)
     Hsym_basis(:,:,2)=kron( pauli_sigma_x, pauli_tau_0) ; LambdaVec(:,2)= sb_field
     if(Nsym==3)then
        Hsym_basis(:,:,3)=kron( pauli_sigma_x, pauli_tau_x)
        LambdaVec(:,3)= sb_field
     endif
  case default;stop "Norb>2 not supported in this code"
  end select


  select case(bath_type)
  case default
     Nb   = ed_get_bath_dimension()
  case("replica")
     call ed_set_Hreplica(Hsym_basis,LambdaVec)
     Nb   = ed_get_bath_dimension(Nsym)
  case("general")
     call ed_set_Hgeneral(Hsym_basis,LambdaVec)
     Nb   = ed_get_bath_dimension(Nsym)
  end select

  !set Hloc
  call ed_set_Hloc(hloc)


  allocate(bath(Nb))
  allocate(prev(Nb))  !setup solver
  call ed_init_solver(bath)



  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)
     call ed_get_sigma(Smats(1,:,:,:,:,:),axis="m",type="n")
     call ed_get_sigma(Smats(2,:,:,:,:,:),axis="m",type="a")


     !Compute the local gfs:
     call get_gloc(Hk,Gmats,Smats,axis='m')
     call write_gf(Gmats(1,:,:,:,:,:),"Gloc",axis='mats',iprint=1)
     call write_gf(Gmats(2,:,:,:,:,:),"Floc",axis='mats',iprint=1)


     call dmft_self_consistency(&
          Gmats(1,:,:,:,:,:),Gmats(2,:,:,:,:,:),&
          Smats(1,:,:,:,:,:),Smats(2,:,:,:,:,:),&
          Weiss(1,:,:,:,:,:),Weiss(2,:,:,:,:,:))
     call write_gf(Weiss(1,:,:,:,:,:),"Weiss",axis='mats',iprint=1)
     call write_gf(Weiss(2,:,:,:,:,:),"Theta",axis='mats',iprint=1)


     !Perform the SELF-CONSISTENCY by fitting the new bath
     call ed_chi2_fitgf(Weiss(1,:,:,:,:,:),Weiss(2,:,:,:,:,:),bath,ispin=1)
     if(phsym)call ed_ph_symmetrize_bath(bath)



     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Prev
     Prev=Bath

     !Check convergence (if required change chemical potential)
     converged = check_convergence(Weiss(1,1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)

     call end_loop
  enddo

  !Compute the local gfs:
  call ed_get_sigma(Sreal(1,:,:,:,:,:),axis="r",type="n")
  call ed_get_sigma(Sreal(2,:,:,:,:,:),axis="r",type="a")
  call get_gloc(Hk,Greal,Sreal,axis='r')
  call write_gf(Greal(1,:,:,:,:,:),"Gloc",axis='real',iprint=1)
  call write_gf(Greal(2,:,:,:,:,:),"Floc",axis='real',iprint=1)


  ! !Compute the Kinetic Energy:
  ! call dmft_kinetic_energy(Hk(1,:,:,:),Smats(1,:,:,:,:,:),Smats(2,:,:,:,:,:))

  call Finalize_MPI()

contains




  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Hk model for the 2d square lattice
  !-------------------------------------------------------------------------------------------
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:) :: kpoint
    integer              :: N
    real(8)              :: kx,ky,t(N)
    complex(8)           :: hk(N,N)
    if(N>2)stop "hk_model error: N>2"
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = zero
    Hk(1,1) = -one*2d0*ts*(cos(kx)+cos(ky))
    if(Norb==2)then
       Hk(1,1) = -one*2d0*ts*(cos(kx)+cos(ky))       + delta 
       Hk(2,2) = -one*2d0*alpha*ts*(cos(kx)+cos(ky)) - delta
       Hk(1,2) = -one*lambda*(cos(kx)-cos(ky))
       Hk(2,1) = conjg(Hk(1,2))
    endif
  end function hk_model





end program ed_ahm_square



