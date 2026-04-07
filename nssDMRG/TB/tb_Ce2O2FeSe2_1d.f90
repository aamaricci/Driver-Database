program Ce2O2FeSe2_1d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer,parameter                           :: Norb=3
  integer,parameter                           :: Nspin=1
  integer,parameter                           :: Nso=Norb*Nspin,L=2048
  integer                                     :: Nkpath
  integer                                     :: Nlat,Nx,Nup,Ndw
  integer                                     :: i,j,k,ik,iorb,jorb,ispin,io,jo
  integer                                     :: ilat,jlat
  integer                                     :: ix
  real(8)                                     :: kx
  real(8),dimension(:,:),allocatable          :: kgrid,kpath
  integer,dimension(:,:),allocatable          :: Links
  complex(8),dimension(:,:,:),allocatable     :: Hk
  complex(8),dimension(:,:),allocatable       :: Hij,Hloc
  complex(8),dimension(:,:,:,:),allocatable   :: Hlat
  real(8),dimension(:),allocatable            :: Eij
  real(8),dimension(:),allocatable            :: rhoDiag
  complex(8),dimension(:,:),allocatable       :: rhoH
  complex(8),dimension(:,:,:,:,:),allocatable :: Gmats,Greal
  real(8),dimension(:,:),allocatable          :: dens
  real(8)                                     :: mh,lambda,Ef
  real(8)                                     :: xmu,eps,beta,wmin,wmax,filling
  character(len=20)                           :: file
  logical                                     :: pbc
  character(len=20)                       :: w90file

  call parse_input_variable(w90file,"w90file","input.conf",default="Ce2O2FeSe2_hr.w90")
  call parse_input_variable(filling,"FILLING","input.conf",default=dble(Nso))
  call parse_input_variable(nx,"NX","input.conf",default=200)
  call parse_input_variable(nkpath,"NKPATH","input.conf",default=500)
  call parse_input_variable(pbc,"pbc","input.conf",default=.true.)
  call parse_input_variable(beta,"beta","input.conf",default=1000d0)
  call parse_input_variable(wmin,"WMIN","input.conf",default=-2d0)
  call parse_input_variable(wmax,"WMAX","input.conf",default=2d0)
  call parse_input_variable(eps,"EPS","input.conf",default=0.001d0)
  call save_input_file("input.conf")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(0d0,"xmu")
  call add_ctrl_var(wmin,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")

  call TB_set_ei([1d0])
  call TB_build_bk(verbose=.true.)

  call TB_w90_setup(reg(w90file),Nlat=[1],norb=[Norb],verbose=.true.)
  if(filling/=0d0)then
    call TB_w90_FermiLevel([Nx],filling,Ef)
    print*,"Ef=",Ef
  endif

  !First solve the H(k) problem (use PBC here)
  call TB_build_model(Hk,Nso,[Nx]) !W90 interface
  Hloc=sum(Hk,dim=3)/Nx
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc,"w90Hloc.dat")
  deallocate(Hloc)

  !SOLVE ALONG A PATH IN THE BZ.
  allocate(kpath(3,3))
  kpath(1,:)=-[0.5d0]
  kpath(2,:)= [0d0]
  kpath(3,:)= [0.5d0]
  call TB_Solve(Nso,kpath,Nkpath,&
       colors_name=[red,green,blue],&
       points_name=[character(len=20) :: '-X', 'G', 'X'],&
       file="W90Eigenband.dat")


  !Build the local GF:
  allocate(Gmats(Nspin,Nspin,Norb,Norb,L))
  allocate(Greal(Nspin,Nspin,Norb,Norb,L))
  call get_gloc(Hk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L),axis='mats')
  call get_gloc(Hk,Greal,zeros(Nspin,Nspin,Norb,Norb,L),axis='real')
  call write_gf(Gmats,"Gloc",axis='mats',iprint=1)
  call write_gf(Greal,"Gloc",axis='real',iprint=1)

  !Compute the local density:
  allocate(dens(Nspin,Norb))
  do ispin=1,Nspin
    do iorb=1,Norb
       dens(ispin,iorb) = fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta)
    enddo
 enddo
!plot observables
  open(10,file="w90observables.nint")
  write(10,"(20F20.12)")(dens(1,iorb)*2,iorb=1,Nso),sum(dens)*2
  close(10)
  write(*,"(A,20F14.9)")"w90 Occupations =",(dens(1,iorb)*2,iorb=1,Nso),sum(dens)*2


call TB_w90_delete()

end program Ce2O2FeSe2_1d


