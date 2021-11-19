!sigma = tgz*Tz
program testMF
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                          :: info,unit,Len,i,N
  real(8)                          :: mh,t,v,tol,tz,tzdens
  real(8)                          :: gsigma,sigma,x(1),dx(1)
  real(8)                          :: integral
  logical                          :: minimize,withgf,withbands
  logical                          :: bool
  character(len=20)                :: finput
  real(8)                          :: smin,smax,ky_g
  real(8),dimension(:),allocatable :: sarray
  !
  ! 
  integer,parameter                           :: Norb=2,Nspin=2,Nso=Nspin*Norb
  integer                                     :: Nkx,L,Nktot,Npts
  integer                                     :: iorb,ispin,io
  complex(8),dimension(Nso,Nso)               :: Gamma1,Gamma2,Gamma5
  complex(8),dimension(:,:,:),allocatable     :: Hk
  complex(8),dimension(:,:,:,:,:),allocatable :: Greal
  real(8)                                     :: dens(Nso)
  real(8),dimension(:,:),allocatable          :: kpath

  call parse_cmd_variable(finput,"FINPUT",default="inputMF.conf")
  call parse_input_variable(t,"t",finput,default=0.2d0)
  call parse_input_variable(gsigma,"gsigma",finput,default=0.5d0)
  call parse_input_variable(mh,"MH",finput,default=1.d0)
  call parse_input_variable(v,"v",finput,default=0.3d0)
  call parse_input_variable(sigma,"sigma",finput,default=0.05d0)
  call parse_input_variable(minimize,"minimize",finput,default=.false.)
  call parse_input_variable(len,"len",finput,default=100)
  call parse_input_variable(smin,"smin",finput,default=-2d0)
  call parse_input_variable(smax,"smax",finput,default=2d0)
  !
  call parse_input_variable(Nkx,"Nkx",finput,default=31)
  call parse_input_variable(N,"N",finput,default=50)
  call parse_input_variable(L,"L",finput,default=1024)
  call parse_input_variable(tol,"tol",finput,default=1d-4)
  call parse_input_variable(withgf,"withgf",finput,default=.false.)
  call parse_input_variable(withbands,"withbands",finput,default=.false.)
  !
  inquire(file="sigma.restart",exist=bool)
  if(bool)then
     print*,"reading sigma"
     open(free_unit(unit),file="sigma.restart")
     read(unit,*)sigma
     close(unit)
  endif
  call save_input(trim(finput))

  call add_ctrl_var(1000d0,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(0d0,"xmu")
  call add_ctrl_var(-10d0,"wini")
  call add_ctrl_var(10d0,"wfin")
  call add_ctrl_var(0.04d0,"eps")


  x(1)=sigma
  dx(1)=0.01d0
  call fmin(bhz_f,x,lambda=dx)
  sigma=x(1)

  open(free_unit(unit),file="sigma.restart")
  write(unit,*)sigma
  close(unit)


  if(.not.minimize)then
     open(free_unit(unit),file="fsigma.dat")
     allocate(sarray(len))
     sarray = linspace(smin,smax,len)
     do i=1,len
        x(1) = sarray(i)
        write(unit,*)x(1),bhz_f(x)
     enddo
     close(unit)
  endif


  Nktot= Nkx*Nkx
  !SETUP THE GAMMA MATRICES:
  gamma1=kron_pauli( pauli_tau_z, pauli_sigma_x)
  gamma2=kron_pauli( pauli_tau_0,-pauli_sigma_y)
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z)
  call TB_set_ei([1d0,0d0],[0d0,1d0])
  call TB_set_bk([pi2,0d0],[0d0,pi2])


  allocate(Hk(Nso,Nso,Nktot))
  call TB_build_model(Hk,hk_model,Nso,[Nkx,Nkx])
  dens = get_dens(0d0)
  open(10,file="observables.nint")
  write(10,"(20F20.12)")(dens(iorb),iorb=1,Nso),sum(dens)
  close(10)
  write(*,"(A,20F14.9)")"Occupations =",(dens(iorb),iorb=1,Nso),sum(dens)

  tz = dens(1)-dens(2)


  open(free_unit(unit),file="sigma.dat")
  write(unit,*)mh,gsigma,sigma,tz*gsigma
  close(unit)

  open(free_unit(unit),file="meff.dat")
  write(unit,*)mh,gsigma,mh-sigma,mh-gsigma*tz
  close(unit)

contains

  function bhz_f(a) result(f)
    real(8),dimension(:) :: a
    real(8)              :: f
    real(8)              :: integral
    sigma = a(1)
    call quad2d(N,integral)
    f = (sigma**2)/gsigma - 2d0*integral
  end function bhz_f

  function hk_bhz(kvec) result(hk)
    real(8),dimension(:) :: kvec
    real(8)              :: hk
    real(8)              :: ek,x2,y2,kx,ky,meff
    kx=kvec(1)
    ky=kvec(2)
    ek = -2d0*t*(cos(kx)+cos(ky))
    x2 =  v*sin(kx);x2=x2**2
    y2 =  v*sin(ky);y2=y2**2
    meff= mh + sigma!*gsigma
    Hk = sqrt( (Meff + ek)**2 + (x2+y2) )
  end function hk_bhz

  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: ek
    real(8)                   :: kx,ky,meff
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    ek = -2d0*t*(cos(kx)+cos(ky))
    Meff= mh + sigma!*gsigma
    Hk = (Meff + ek)*Gamma5 + v*sin(kx)*Gamma1 + v*sin(ky)*Gamma2
  end function hk_model
  
  function hk_x(kx) result(fx)
    real(8) :: kx
    real(8) :: fx
    fx =  hk_bhz([kx,ky_g])
  end function hk_x

  subroutine quad2d(Ly,int)
    integer,optional                 :: Ly
    integer                          :: Ly_
    real(8)                          :: int,inty
    real(8),dimension(:),allocatable :: Yarray
    Ly_=50;if(present(Ly))Ly_=Ly
    !
    allocate(Yarray(Ly_))
    Yarray = linspace(0d0,pi2,Ly_)
    int=0d0
    do i=1,Ly_
       ky_g = Yarray(i)
       call quad(hk_x,0d0,pi2,result=inty)
       int=int+inty/Ly_/pi2
    enddo
  end subroutine quad2d



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
end program





