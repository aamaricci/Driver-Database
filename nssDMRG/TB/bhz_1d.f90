program bhz_1d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer,parameter                           :: Norb=2
  integer,parameter                           :: Nspin=2
  integer,parameter                           :: Nso=Norb*Nspin,L=2048
  integer                                     :: Nkpath
  integer                                     :: Nlat,Nx,Nup,Ndw
  integer                                     :: i,j,k,ik,iorb,jorb,ispin,io,jo
  integer                                     :: ilat,jlat
  integer                                     :: ix
  real(8)                                     :: kx
  real(8),dimension(:,:),allocatable          :: kgrid,kpath,ktrims,Rgrid
  integer,dimension(:,:),allocatable          :: Links
  complex(8),dimension(:,:,:),allocatable     :: Hk
  complex(8),dimension(:,:),allocatable       :: Hij,Hloc
  complex(8),dimension(:,:,:,:),allocatable   :: Hlat
  real(8),dimension(:),allocatable            :: Eij
  real(8),dimension(:),allocatable            :: rhoDiag
  complex(8),dimension(:,:),allocatable       :: rhoH
  complex(8),dimension(:,:,:,:,:),allocatable :: Gmats,Greal
  real(8),dimension(:,:),allocatable          :: dens
  real(8)                                     :: mh,lambda,E0
  real(8)                                     :: xmu,beta
  character(len=20)                           :: file
  logical                                     :: pbc
  complex(8),dimension(Nso,Nso)               :: GammaX,Gamma0,GammaZ

  call parse_input_variable(nx,"NX","input.conf",default=100)
  call parse_input_variable(mh,"MH","input.conf",default=0d0)
  call parse_input_variable(lambda,"LAMBDA","input.conf",default=0.3d0)
  call parse_input_variable(pbc,"pbc","input.conf",default=.false.)
  call parse_input_variable(beta,"beta","input.conf",default=1000d0)
  call save_input_file("input.conf")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(0d0,"xmu")
  call add_ctrl_var(-5d0,"wini")
  call add_ctrl_var(5d0,"wfin")
  call add_ctrl_var(0.02d0,"eps")


  !SETUP THE GAMMA MATRICES:
  gammaX=kron( pauli_sigma_z, pauli_tau_x)
  gamma0=kron( pauli_sigma_0, pauli_tau_0)
  gammaZ=kron( pauli_sigma_0, pauli_tau_z)


  call TB_set_ei([1d0,0d0])
  call TB_set_bk([pi2,0d0])


  !First solve the H(k) problem (use PBC here)
  allocate(Hk(Nso,Nso,Nx))
  call TB_build_model(Hk,hk_model,Nso,[Nx])
  Hloc=sum(Hk,dim=3)/Nx
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc)
  deallocate(Hloc)
  !SOLVE ALONG A PATH IN THE BZ.
  allocate(kpath(3,3))
  kpath(1,:)=-kpoint_X1
  kpath(2,:)= kpoint_Gamma
  kpath(3,:)= kpoint_X1
  call TB_Solve_model(hk_model,Nso,kpath,500,&
       colors_name=[red1,blue1,red1,blue1],&
       points_name=[character(len=20) :: '-X', 'G', 'X'],&
       file="Eigenband.dat")


  !Build the local GF:
  allocate(Gmats(Nspin,Nspin,Norb,Norb,L))
  allocate(Greal(Nspin,Nspin,Norb,Norb,L))
  call get_gloc(Hk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L),axis='mats')
  call get_gloc(Hk,Greal,zeros(Nspin,Nspin,Norb,Norb,L),axis='real')
  call write_gf(Gmats,"Gloc",axis='mats',iprint=1)
  call write_gf(Greal,"Gloc",axis='real',iprint=1)




  allocate(Links(2,1))          !Links: right,left
  Links(1,:) = [1]
  Links(2,:) =-[1]


  if(pbc)then
     open(100,file="bhz_obsVSnlat_PBC.dat")
  else
     open(100,file="bhz_obsVSnlat_OBC.dat")
  endif

  do Nlat=2,Nx,2

     if(allocated(Hlat))deallocate(Hlat)
     allocate(Hlat(Nso,Nso,Nlat,Nlat))
     Hlat=zero

     !>Build direct Lattice Hamiltonian
     call TB_build_model(Hlat,ts_model,Nso,[Nlat],Links,pbc=pbc)
     !
     allocate(Hij(Nlat*Nso,Nlat*Nso))
     allocate(rhoH(Nlat*Nso,Nlat*Nso))
     allocate(Eij(Nlat*Nso))
     allocate(rhoDiag(Nlat*Nso))
     allocate(dens(Nlat,Nso))
     !
     Hij = zero
     do concurrent(ilat=1:Nlat,jlat=1:Nlat,io=1:Nso,jo=1:Nso)
        i = io + (ilat-1)*Nso        
        j = jo + (jlat-1)*Nso
        Hij(i,j) = Hlat(io,jo,ilat,jlat)
     enddo
     !
     call eigh(Hij,Eij)
     !
     rhoDiag = fermi(Eij,beta)
     rhoH    = matmul(Hij , matmul(diag(rhoDiag), conjg(transpose(Hij))) ) 
     do concurrent(ilat=1:Nlat,io=1:Nso)
        i = io + (ilat-1)*Nso
        dens(ilat,io) = rhoH(i,i)
     enddo
     !     
     Nup = Nlat*Nso/2
     Ndw = Nlat*Nso-Nup
     !
     E0 = sum(Eij(:Nup)) + sum(Eij(:Ndw))
     E0 = E0/Nlat/Nso
     !
     write(*,*)Nlat,E0,(sum(dens(:,i))/Nlat,i=1,4),sum(dens)/Nlat
     write(100,*)Nlat,E0,(sum(dens(:,i))/Nlat,i=1,4),sum(dens)/Nlat
     !
     deallocate(Hij,Eij,rhoDiag,rhoH,dens)
  enddo

  close(100)


contains

  

  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: ek
    real(8)                   :: kx
    complex(8),dimension(N,N) :: hk
    kx = kpoint(1)
    ek = -1d0*cos(kx)
    Hk = (Mh+ek)*GammaZ + lambda*sin(kx)*GammaX
  end function hk_model


  function ts_model(link,Nso) result(Hts)
    integer                       :: link
    integer                       :: Nso
    complex(8),dimension(Nso,Nso) :: Hts
    select case(link)
    case (0) !LOCAL PART
       Hts =  Mh*GammaZ
    case (1) !RIGHT HOPPING
       Hts = -0.5d0*GammaZ + xi*0.5d0*lambda*GammaX
    case (2) !LEFT HOPPING
       Hts = -0.5d0*GammaZ - xi*0.5d0*lambda*GammaX
    case default 
       stop "ts_model ERROR: link != [0:2]"
    end select
  end function ts_model










end program bhz_1d


