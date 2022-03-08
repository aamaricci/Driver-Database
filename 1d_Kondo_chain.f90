program Kondo1d
  USE EDLAT
  USE SCIFOR
#ifdef _MPI
  USE MPI
#endif
  implicit none
  character(len=16)   :: finput
  real(8)             :: ts,gamma
  integer             :: N,N1,i
  logical             :: pbc
  !
#ifdef _MPI
  call init_MPI
#endif
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="chain hopping parameter")
  call parse_input_variable(pbc,"PBC",finput,default=.false.,comment="T: PBC, F: OBC")
  call parse_input_variable(gamma,"GAMMA",finput,default=0.01d0,comment="field amplitude")
  call ed_read_input(trim(finput))
  !

  if(any(Nsites(1:Norb)==0))stop "This driver is for Kondo problem only: Nsites=[1,N]"


  call ed_Hij_init(Nsites(1:Norb))

  call solve_drude_u0()
  !
  N  = Nsites(1)                !odd
  N1 = (N+1)/2                  !N1%2==0
  !                I 
  !                o
  !       1--2--3--4--5--6--7--8
  call ed_Hij_add_link(1,2,1,1,1,one*ts)
  do i=2,N-1
     call ed_Hij_add_link(i,i-1,1,1,1,one*ts)
     call ed_Hij_add_link(i,i+1,1,1,1,one*ts)
  enddo
  call ed_Hij_add_link(N,N-1,1,1,1,one*ts)
  !
  !PBC
  if(pbc)then
     call ed_Hij_add_link(1,N,1,1,1,one*ts)
     call ed_Hij_add_link(N,1,1,1,1,one*ts)
  end if
  !
  call ed_Hij_info()
  Jkindx = N1


  call ed_init_solver()
  call ed_solve()


#ifdef _MPI
  call finalize_MPI()
#endif


contains

  subroutine solve_drude_u0()
    real(8)             :: kx,ek,Icurrent,beta
    integer             :: N,N1,i,Tlen,it
    logical             :: tbool
    real(8),allocatable :: temperature_list(:)
    integer,allocatable :: Tord(:)
    !
    N  = Nsites(2)                !odd
    N1 = (N+1)/2                  !N1%2==0
    !
    inquire(file="temperature.restart",exist=Tbool)
    if(Tbool)then
       write(LOGfile,"(A)")"Reading temperature list from file temperature.restart"
       Tlen = file_length("temperature.restart")
       open(100,file="temperature.restart")
       allocate(temperature_list(Tlen),Tord(Tlen))
       do i=1,Tlen
          read(100,*)temperature_list(i)
       enddo
       close(100)
       call sort(temperature_list,Tord)                !sort from smallest to largest
       temperature_list = temperature_list(Tlen:1:-1) !invert order
    else
       Tlen=1
       allocate(temperature_list(Tlen))
       temperature_list = temp
    endif
    !
    do N=6,16
       open(100,file="drude_u0_N"//str(N)//".nint")
       do it=1,Tlen
          beta = 1d0/temperature_list(it)
          Icurrent=0d0
          do i=1,N
             kx = (i-1)*pi2/N
             ek = -2*ts*cos(kx+gamma)
             Icurrent = Icurrent - 2*ts*sin(kx+gamma)*fermi(ek,beta)/N*pi2
          enddo
          write(100,*)temperature_list(it),Icurrent/gamma,Icurrent,gamma
       enddo
       close(100)
    enddo
  end subroutine solve_drude_u0

end program Kondo1d













! N  = Nsites(2)                !odd
! N1 = (N+1)/2                  !N1%2==0


! !                         1 
! !                         |
! ![(N1+1=]6--7--8--9[==N]--1--2--3--4--5[=N1]
! !x       1--2--3--4--   --5--6--7--8--9
! !imp-chain
! if(t0/=0d0)then
!    call ed_Hij_add_link(1,1,1,2,1,one*t0)
!    call ed_Hij_add_link(1,1,2,1,1,one*t0)
! endif
! !
! select case(N)     
! case default
!    !                         1 
!    !                         |
!    ![(N1+1=]6--7--8--9[==N]--1--2--3--4--5[=N1]
!    !x       1--2--3--4--   --5--6--7--8--9
!    call ed_Hij_add_link(1,2,2,2,1,one*ts)
!    call ed_Hij_add_link(1,N,2,2,1,one*ts)
!    do i=2,N1-1
!       call ed_Hij_add_link(i,i-1,2,2,1,one*ts)
!       call ed_Hij_add_link(i,i+1,2,2,1,one*ts)
!    enddo
!    call ed_Hij_add_link(N1,N1-1,2,2,1,one*ts)
!    !
!    call ed_Hij_add_link(N1+1,N1+2,2,2,1,one*ts)
!    do i=N1+2,N-1
!       call ed_Hij_add_link(i,i-1,2,2,1,one*ts)
!       call ed_Hij_add_link(i,i+1,2,2,1,one*ts)
!    enddo
!    call ed_Hij_add_link(N,N-1,2,2,1,one*ts)
!    call ed_Hij_add_link(N,1,2,2,1,one*ts)
!    !
!    if(pbc)then
!       call ed_Hij_add_link(N1,N1+1,2,2,1,one*ts)
!       call ed_Hij_add_link(N1+1,N1,2,2,1,one*ts)
!    end if
!    !
! case (3)
!    !                1 
!    !                |
!    !        3[==N]--1--2
!    call ed_Hij_add_link(1,2,2,2,1,one*ts)
!    call ed_Hij_add_link(2,1,2,2,1,one*ts)
!    call ed_Hij_add_link(1,3,2,2,1,one*ts)
!    call ed_Hij_add_link(3,1,2,2,1,one*ts)
!    if(pbc)then
!       call ed_Hij_add_link(2,3,2,2,1,one*ts)
!       call ed_Hij_add_link(3,2,2,2,1,one*ts)
!    end if
! end select
! call ed_Hij_info()

! write(100,*)N1,2
! do i=N1,N
!    write(100,*)i,1
! enddo
! do i=1,N1-1
!    write(100,*)i,1
! enddo

