!       \                 /
!          \           / 
!             \  _  / 
!       1--2--3--4--5--6--7
!                o
!                I 
!
program Kondo1d
  USE EDLAT
  USE SCIFOR
#ifdef _MPI
  USE MPI
#endif
  implicit none
  character(len=16)   :: finput
  real(8)             :: ts,vtrap,atrap
  integer             :: N,N1,i
  logical             :: pbc
  !
#ifdef _MPI
  call init_MPI
#endif
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="chain hopping parameter")
  call parse_input_variable(atrap,"atrap",finput,default=0d0,comment="optical trap bottom energy ")
  call parse_input_variable(vtrap,"vtrap",finput,default=0d0,comment="optical trap potential amplitude: 1/2.V.x**2")
  call parse_input_variable(pbc,"PBC",finput,default=.false.,comment="T: PBC, F: OBC")
  call ed_read_input(trim(finput))
  !

  if(any(Nsites(1:Norb)==0))stop "This driver is for Kondo problem only: Nsites=[1,N]"


  call ed_Hij_init(Nsites(1:Norb))

  N  = Nsites(1)                !odd
  N1 = (N+1)/2                  !N1%2==0



  call ed_Hij_add_link(1,2,1,1,1,one*ts)
  do i=2,N-1
     call ed_Hij_add_link(i,i-1,1,1,1,one*ts)
     call ed_Hij_add_link(i,i+1,1,1,1,one*ts)
  enddo
  call ed_Hij_add_link(N,N-1,1,1,1,one*ts)
  !PBC
  if(pbc)then
     call ed_Hij_add_link(1,N,1,1,1,one*ts)
     call ed_Hij_add_link(N,1,1,1,1,one*ts)
  end if
  do i=1,N
     write(100,*)i,0,e_local(i)
     call ed_Hij_add_link(i,i,1,1,1,one*e_local(i))
  enddo

  !
  call ed_Hij_info()
  Jkindx = N1


  call ed_init_solver()
  call ed_solve()


#ifdef _MPI
  call finalize_MPI()
#endif

contains
  function e_local(i) result(e)
    integer :: i
    real(8) :: e
    real(8) :: r2
    r2 = dble(i-N1)**2
    e  = atrap + 0.5d0*Vtrap*r2
  end function e_local

end program Kondo1d













! contains

!   subroutine solve_drude_u0()
!     real(8)             :: kx,ek,Icurrent,beta
!     integer             :: N,N1,i,Tlen,it
!     logical             :: tbool
!     real(8),allocatable :: temperature_list(:)
!     integer,allocatable :: Tord(:)
!     !
!     N  = Nsites(2)                !odd
!     N1 = (N+1)/2                  !N1%2==0
!     !
!     inquire(file="temperature.restart",exist=Tbool)
!     if(Tbool)then
!        write(LOGfile,"(A)")"Reading temperature list from file temperature.restart"
!        Tlen = file_length("temperature.restart")
!        open(100,file="temperature.restart")
!        allocate(temperature_list(Tlen),Tord(Tlen))
!        do i=1,Tlen
!           read(100,*)temperature_list(i)
!        enddo
!        close(100)
!        call sort(temperature_list,Tord)                !sort from smallest to largest
!        temperature_list = temperature_list(Tlen:1:-1) !invert order
!     else
!        Tlen=1
!        allocate(temperature_list(Tlen))
!        temperature_list = temp
!     endif
!     !
!     do N=6,16
!        open(100,file="drude_u0_N"//str(N)//".nint")
!        do it=1,Tlen
!           beta = 1d0/temperature_list(it)
!           Icurrent=0d0
!           do i=1,N
!              kx = (i-1)*pi2/N
!              ek = -2*ts*cos(kx+gamma)
!              Icurrent = Icurrent - 2*ts*sin(kx+gamma)*fermi(ek,beta)/N*pi2
!           enddo
!           write(100,*)temperature_list(it),Icurrent/gamma,Icurrent,gamma
!        enddo
!        close(100)
!     enddo
!   end subroutine solve_drude_u0








