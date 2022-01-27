!                             1==imp 
!                             |
![(N+1)/2+1=]6--7--8--9[==N]--1--2--3--4--5[=N+1/2]
program Kondo1d
  USE EDLAT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  character(len=16)    :: finput
  real(8)              :: ts,t0,gamma,kx,ek,Icurrent,beta
  integer              :: N,N1,i
  logical              :: pbc
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(t0,"T0",finput,default=0.d0,comment="imp-chain coupling")
  call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="chain hopping parameter")
  call parse_input_variable(pbc,"PBC",finput,default=.false.,comment="T: PBC, F: OBC")
  call parse_input_variable(gamma,"GAMMA",finput,default=0.01d0,comment="field amplitude")
  call ed_read_input(trim(finput))
  !
  beta = 1d0/temp
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  if(Nspin/=1.OR.Norb/=2)stop "This driver is for Kondo problem only: Norb=2, Nspin=1"
  if(any(Nsites(1:Norb)==0))stop "This driver is for Kondo problem only: Nsites=[1,N]"
  if(Nsites(1)/=1)stop "This driver is for Kondo problem only: Nsites=[1,N]"
  if(mod(Nsites(2),2)==0)stop "This driver is for Kondo problem only: Nsites=[1,N], with N odd to ensure chain symmetry"


  !1d chain with OBC and staggered energies
  call ed_Hij_init(Nsites(1:Norb))
  !chain

  N  = Nsites(2)                !odd
  N1 = (N+1)/2                  !N1%2==0

  Icurrent=0d0
  do i=1,N
     kx = i*pi2/N
     ek = -2*abs(ts)*cos(kx+gamma)
     Icurrent = Icurrent + 2*abs(ts)*sin(kx+gamma)*fermi(ek,beta)/N*pi2
  enddo
  open(free_unit(i),file="drude_u0.nint")
  write(i,*)Icurrent,gamma,Icurrent/gamma
  close(i)

  !                         1 
  !                         |
  ![(N1+1=]6--7--8--9[==N]--1--2--3--4--5[=N1]
  !x       1--2--3--4--   --5--6--7--8--9
  !imp-chain
  ! if(t0/=0d0)then
  call ed_Hij_add_link(1,1,1,2,1,one*t0)
  call ed_Hij_add_link(1,1,2,1,1,one*t0)
  ! endif
  !
  select case(N)     
  case default
     !                         1 
     !                         |
     ![(N1+1=]6--7--8--9[==N]--1--2--3--4--5[=N1]
     !x       1--2--3--4--   --5--6--7--8--9
     call ed_Hij_add_link(1,2,2,2,1,one*ts)
     call ed_Hij_add_link(1,N,2,2,1,one*ts)
     do i=2,N1-1
        call ed_Hij_add_link(i,i-1,2,2,1,one*ts)
        call ed_Hij_add_link(i,i+1,2,2,1,one*ts)
     enddo
     call ed_Hij_add_link(N1,N1-1,2,2,1,one*ts)
     !
     call ed_Hij_add_link(N1+1,N1+2,2,2,1,one*ts)
     do i=N1+2,N-1
        call ed_Hij_add_link(i,i-1,2,2,1,one*ts)
        call ed_Hij_add_link(i,i+1,2,2,1,one*ts)
     enddo
     call ed_Hij_add_link(N,N-1,2,2,1,one*ts)
     call ed_Hij_add_link(N,1,2,2,1,one*ts)
     !
     if(pbc)then
        call ed_Hij_add_link(N1,N1+1,2,2,1,one*ts)
        call ed_Hij_add_link(N1+1,N1,2,2,1,one*ts)
     end if
     !
  case (3)
     !                1 
     !                |
     !        3[==N]--1--2
     call ed_Hij_add_link(1,2,2,2,1,one*ts)
     call ed_Hij_add_link(2,1,2,2,1,one*ts)
     call ed_Hij_add_link(1,3,2,2,1,one*ts)
     call ed_Hij_add_link(3,1,2,2,1,one*ts)
     if(pbc)then
        call ed_Hij_add_link(2,3,2,2,1,one*ts)
        call ed_Hij_add_link(3,2,2,2,1,one*ts)
     end if
  end select
  call ed_Hij_info()

  write(100,*)N1,2
  do i=N1,N
     write(100,*)i,1
  enddo
  do i=1,N1-1
     write(100,*)i,1
  enddo


  call ed_init_solver()
  call ed_solve()



end program Kondo1d






