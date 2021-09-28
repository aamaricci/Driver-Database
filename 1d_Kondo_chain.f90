!                             1 
!                             |
![(N+1)/2+1=]6--7--8--9[==N]--1--2--3--4--5[=N+1/2]
program Kondo1d
  USE EDLAT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  character(len=16)    :: finput
  real(8)              :: ts,t0
  integer              :: N,N1,i
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(t0,"T0",finput,default=0.d0,comment="imp-chain coupling")
  call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="chain hopping parameter")
  call ed_read_input(trim(finput))
  !
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
  call ed_Hij_init(Nsites)
  !chain

  N  = Nsites(2)                !odd
  N1 = (N+1)/2                  !N1%2==0

  !                             1 
  !                             |
  ![(N+1)/2+1=]6--7--8--9[==N]--1--2--3--4--5[=N+1/2]
  !imp-chain
  call ed_Hij_add_link(1,1,1,2,1,one*t0)
  call ed_Hij_add_link(1,1,2,1,1,one*t0)
  !
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
  call ed_Hij_info()
  call ed_Hij_write()


  call ed_init_solver()
  call ed_solve()



end program Kondo1d






