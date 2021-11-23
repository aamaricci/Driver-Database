program testMF
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                          :: dim,info,unit,Len,i
  real(8)                          :: mh,b,tol,v,tz,tzdens,sigmadens
  real(8)                          :: gint,sigma,x(1),dx(1)
  real(8)                          :: integral
  logical                          :: minimize,withgf,withbands
  logical                          :: bool
  character(len=20)                :: finput
  real(8)                          :: smin,smax
  real(8),dimension(:),allocatable :: sarray
  integer,parameter                           :: Norb=2,Nspin=1,Nso=Nspin*Norb
  integer                                     :: Nkx,L,Nktot,Npts
  integer                                     :: iorb,ispin,io
  complex(8),dimension(:,:,:),allocatable     :: Hk
  complex(8),dimension(:,:,:,:,:),allocatable :: Greal
  real(8)                                     :: dens(Nso)
  real(8),dimension(:,:),allocatable          :: kpath

  call parse_cmd_variable(finput,"FINPUT",default="inputMF.conf")
  call parse_input_variable(b,"b",finput,default=0.2d0)
  call parse_input_variable(gint,"gint",finput,default=0.5d0)
  call parse_input_variable(mh,"MH",finput,default=0.05d0)
  call parse_input_variable(sigma,"sigma",finput,default=0.05d0)
  call parse_input_variable(minimize,"minimize",finput,default=.false.)
  call parse_input_variable(len,"len",finput,default=1000)
  call parse_input_variable(smin,"smin",finput,default=-1d0)
  call parse_input_variable(smax,"smax",finput,default=1d0)
  !
  call parse_input_variable(Nkx,"Nkx",finput,default=20)
  call parse_input_variable(L,"L",finput,default=512)
  call parse_input_variable(dim,"dim",finput,default=2)
  call parse_input_variable(v,"v",finput,default=0.3d0)
  call parse_input_variable(tol,"tol",finput,default=1d-6)
  call parse_input_variable(withgf,"withgf",finput,default=.false.)
  call parse_input_variable(withbands,"withbands",finput,default=.false.)
  !
  inquire(file="sigma.restart",exist=bool)
  if(bool)then
     open(free_unit(unit),file="sigma.restart")
     read(unit,*)sigma
     close(unit)
  endif
  call save_input(trim(finput))

  call add_ctrl_var(1000d0,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(0d0,"xmu")
  call add_ctrl_var(-1d0,"wini")
  call add_ctrl_var(1d0,"wfin")
  call add_ctrl_var(0.001d0,"eps")

  x(1)=sigma
  dx(1)=0.01d0
  call fmin(prb_f,x,lambda=dx)
  sigma=x(1)
  tz = sigma/gint



  if(.not.minimize)then
     open(free_unit(unit),file="fsigma.dat")
     allocate(sarray(len))
     sarray = linspace(smin,smax,len)
     do i=1,len
        x(1) = sarray(i)
        write(unit,*)x(1),prb_f(x)
     enddo
     close(unit)
  endif



  Nktot= Nkx*Nkx
  call TB_set_bk([1d0,0d0],[0d0,1d0])


  allocate(Hk(Nso,Nso,Nktot))
  call TB_build_model(Hk,hk_model,Nso,[Nkx,Nkx])
  dens = get_dens(0d0)
  open(10,file="observables.nint")
  write(10,"(20F20.12)")(dens(iorb),iorb=1,Nso),sum(dens)
  close(10)
  write(*,"(A,20F14.9)")"Occupations =",(dens(iorb),iorb=1,Nso),sum(dens)

  tzdens = dens(1)-dens(2)
  sigmadens = gint*tzdens

  open(free_unit(unit),file="tz.dat")
  write(unit,*)mh,gint,tz,tzdens
  close(unit)


  open(free_unit(unit),file="sigma.dat")
  write(unit,*)mh,gint,sigma,sigmadens
  close(unit)
  
  open(free_unit(unit),file="sigma.restart")
  write(unit,*)sigma
  close(unit)


  !SOLVE ALONG A PATH IN THE BZ.
  if(withbands)then
     Npts=3
     allocate(kpath(Npts,3))
     kpath(1,:)=-kpoint_X1
     kpath(2,:)=kpoint_Gamma
     kpath(3,:)=kpoint_X1
     call TB_Solve_model(Hk_model,Nso,kpath,100,&
          colors_name=[red1,blue1],&
          points_name=[character(len=20) :: '-X', 'G', 'X'],&
          file="Eigenband.nint")
  endif

  if(withgf)then
     !Build the local GF:
     allocate(Greal(Nspin,Nspin,Norb,Norb,L))
     call dmft_gloc_realaxis(Hk,Greal,zeros(Nspin,Nspin,Norb,Norb,L))
     call dmft_print_gf_realaxis(Greal,"Gloc",1)
  endif




contains

  function prb_f(a) result(f)
    real(8),dimension(:) :: a
    real(8)              :: f
    real(8)              :: integral
    sigma = a(1)
    call quad(Fintegrand,0d0,1d0,result=integral)
    f = (sigma**2)/2d0/gint - 2d0*integral
  end function prb_f


  function Fintegrand(x) result(fint)
    real(8) :: x
    real(8) :: fint
    fint = (x**(dim-1))*sqrt(x**2 + (sigma + mh + b*x**2)**2)
    return
  end function Fintegrand



  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: ek
    real(8)                   :: kx,ky,meff
    complex(8),dimension(N,N) :: hk
    if(N/=2)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = (sigma + mh + b*(kx**2+ky**2))*pauli_tau_z + v*kx*pauli_tau_x + v*ky*pauli_tau_y
  end function hk_model

  function get_dens(mu) result(ndens)
    real(8),intent(in)            :: mu
    real(8),dimension(Nso)        :: ndens,Evec
    complex(8),dimension(Nso,Nso) :: Rho
    real(8),dimension(Nso)        :: Rtmp,Efvec
    integer                       :: io,ik
    ndens     = 0d0
    do ik=1,Nktot
       Rho = Hk(:,:,ik)
       call eigh(Rho,Evec)
       Efvec = fermi(Evec-mu,1000d0)
       forall(io=1:Nso)Rtmp(io) = sum( abs(Rho(io,:))**2*Efvec )
       ndens = ndens + Rtmp
    enddo
    ndens=ndens/Nktot
  end function get_dens


end program testMF
