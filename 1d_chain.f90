!1--2--..--N
program chain1d
  USE EDLAT
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  character(len=16)                  :: finput
  real(8)                            :: ts,t0,gamma,kx,ek,Icurrent
  integer                            :: N,N1,i
  logical                            :: gflag,pbc
  integer                            :: comm,rank
  logical                            :: master  
  !
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="chain hopping parameter")
  call parse_input_variable(gamma,"GAMMA",finput,default=0.01d0,comment="field amplitude")
  call parse_input_variable(pbc,"PBC",finput,default=.true.,comment="T: PBC, F: OBC")
  call ed_read_input(trim(finput))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  if(Nspin/=1.OR.Norb/=1)stop "This driver is for 1d chain problem only: Norb=1, Nspin=1"
  if(any(Nsites(1:Norb)==0))stop "This driver is for 1d chain problem only: Nsites=[1,N]"

  N  = Nsites(1)
  
  Icurrent=0d0
  do i=1,N
     kx = i*pi2/N
     ek = -2*abs(ts)*cos(kx+gamma)
     Icurrent = Icurrent + 2*abs(ts)*sin(kx+gamma)*fermi(ek,beta)/N*pi2
  enddo
  open(free_unit(i),file="drude_u0.nint")
  write(i,*)Icurrent,gamma,Icurrent/gamma
  close(i)
  !

  !1d chain with PBC
  call ed_Hij_init(Nsites(1:Norb))
  call ed_Hij_add_link(1,2,1,1,1,one*ts)
  do i=2,N-1
     call ed_Hij_add_link(i,i-1,1,1,1,one*ts)
     call ed_Hij_add_link(i,i+1,1,1,1,one*ts)
  enddo
  call ed_Hij_add_link(N,N-1,1,1,1,one*ts)
  if(pbc)then
     call ed_Hij_add_link(1,N,1,1,1,one*ts)
     call ed_Hij_add_link(N,1,1,1,1,one*ts)
  end if
  if(master)call ed_Hij_info()
  !
  !
  call ed_init_solver()
  call ed_solve()
  !
  call finalize_MPI()
  !

  
end program chain1d






