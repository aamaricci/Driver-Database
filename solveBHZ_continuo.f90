program testMF
  USE SCIFOR
  implicit none
  integer                          :: dim,info,unit,Len,i
  real(8)                          :: mh,b,tol,v
  real(8)                          :: gsigma,sigma,x(1),dx(1)
  real(8)                          :: integral
  logical                          :: minimize
  logical                          :: bool
  character(len=20)                :: finput
  real(8)                          :: smin,smax
  real(8),dimension(:),allocatable :: sarray

  call parse_cmd_variable(finput,"FINPUT",default="inputMF.conf")
  call parse_input_variable(b,"b",finput,default=0.2d0)
  call parse_input_variable(gsigma,"gsigma",finput,default=0.5d0)
  call parse_input_variable(mh,"MH",finput,default=0.05d0)
  call parse_input_variable(sigma,"sigma",finput,default=0.05d0)
  call parse_input_variable(minimize,"minimize",finput,default=.false.)
  call parse_input_variable(len,"len",finput,default=1000)
  call parse_input_variable(smin,"smin",finput,default=-1d0)
  call parse_input_variable(smax,"smax",finput,default=1d0)
  !
  call parse_input_variable(dim,"dim",finput,default=2)
  call parse_input_variable(v,"v",finput,default=0.3d0)
  call parse_input_variable(tol,"tol",finput,default=1d-6)
  !
  inquire(file="sigma.restart",exist=bool)
  if(bool)then
     open(free_unit(unit),file="sigma.restart")
     read(unit,*)sigma
     close(unit)
  endif
  call save_input(trim(finput))

  
  x(1)=sigma
  dx(1)=0.01d0
  call fmin(prb_f,x,lambda=dx)
  sigma=x(1)


  open(free_unit(unit),file="sigma.dat")
  write(unit,*)mh,gsigma,sigma,prb_f([sigma])
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
        write(unit,*)x(1),prb_f(x)
     enddo
     close(unit)
  endif


contains

  function prb_f(a) result(f)
    real(8),dimension(:) :: a
    real(8)              :: f
    real(8)              :: integral
    sigma = a(1)
    call quad(Fintegrand,0d0,1d0,result=integral)
    f = (sigma**2)/2d0/gsigma - 2d0*integral
  end function prb_f


  function Fintegrand(x) result(fint)
    real(8) :: x
    real(8) :: fint
    fint = (x**(dim-1))*sqrt(x**2 + (sigma + mh + b*x**2)**2)
    return
  end function Fintegrand



end program testMF
