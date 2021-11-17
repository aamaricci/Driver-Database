program testMF
  USE SCIFOR

  implicit none
  integer                          :: info,unit,Len,i,N
  real(8)                          :: mh,t,v,tol
  real(8)                          :: gsigma,sigma,x(1),dx(1)
  real(8)                          :: integral
  logical                          :: minimize
  logical                          :: bool
  character(len=20)                :: finput
  real(8)                          :: smin,smax,ky_g
  real(8),dimension(:),allocatable :: sarray
  !

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
  call parse_input_variable(N,"N",finput,default=50)
  call parse_input_variable(tol,"tol",finput,default=1d-4)
  !
  inquire(file="sigma.restart",exist=bool)
  if(bool)then
     print*,"reading sigma"
     open(free_unit(unit),file="sigma.restart")
     read(unit,*)sigma
     close(unit)
  endif
  call save_input(trim(finput))



  x(1)=sigma
  dx(1)=0.01d0
  call fmin(bhz_f,x,lambda=dx)
  sigma=x(1)

  open(free_unit(unit),file="sigma.dat")
  write(unit,*)mh,gsigma,sigma
  close(unit)

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


contains

  function bhz_f(a) result(f)
    real(8),dimension(:) :: a
    real(8)              :: f
    real(8)              :: integral
    integer              :: N0
    sigma = a(1)
    call quad2d(N0,integral)
    f = (sigma**2)/2/gsigma - 2d0*integral
  end function bhz_f

  function hk_bhz(kvec) result(hk)
    real(8),dimension(:) :: kvec
    real(8)              :: hk
    real(8)              :: ek,x2,y2,kx,ky
    kx=kvec(1)
    ky=kvec(2)
    ek = -2d0*t*(cos(kx)+cos(ky))
    x2 =  v*sin(kx);x2=x2**2
    y2 =  v*sin(ky);y2=y2**2
    Hk = sqrt( (sigma + Mh + ek)**2 + (x2+y2) )
  end function hk_bhz

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




end program testMF







! N0 = N
! call simps2d_recursive(hk_bhz,[0d0,pi2],[0d0,pi2],&
!      N0=N0,threshold=tol,int=integral)
! integral = integral/N0/N0
! subroutine simps2d_recursive(func,xrange,yrange,N0,Nstep,threshold,iter,int)
!   interface
!      function func(x)
!        real(8),dimension(:) :: x
!        real(8)              :: func
!      end function func
!   end interface
!   real(8),dimension(2)      :: xrange,yrange
!   integer                   :: N,icount
!   real(8)                   :: eps,int0
!   integer,optional          :: N0,Nstep,iter
!   integer                   :: N0_,Nstep_
!   real(8),optional          :: threshold
!   real(8)                   :: threshold_
!   real(8)                   :: int
!   N0_=50;if(present(N0))N0_=N0
!   Nstep_=10;if(present(Nstep))Nstep_=Nstep
!   threshold_=0.1d0;if(present(threshold))threshold_=threshold
!   !
!   N=N0_
!   eps=1d0
!   icount=1
!   int=simps2d(func,xrange,yrange,N,N)
!   do while (eps>threshold_)
!      icount= icount+1
!      int0  = int
!      N     = N+Nstep_
!      int=simps2d(func,xrange,yrange,N,N)
!      eps=abs(int-int0)/abs(int)
!      if(present(iter))iter=icount
!      if(present(N0))N0=N
!   enddo
! end subroutine simps2d_recursive
