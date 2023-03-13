program hubbard_dimer
  USE EDLAT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none


  character(len=16)    :: finput
  real(8)              :: ts,eloc
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(ts,"TS",finput,default=-0.25d0,comment="hopping parameter")
  call parse_input_variable(eloc,"ELOC",finput,default=0d0,comment="local energies")
  call ed_read_input(trim(finput))
  !
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  if(Nspin/=1.OR.Norb/=1.OR.Nsites(1)/=2)stop "This test driver is for a nonmagnetic, single-orbital, dimer only"


  !2d square with OBC and staggered energies
  call ed_Hij_init(Nsites)
  !> ionic potential
  call ed_Hij_add_link(1,1,1,1,1, one*eloc)
  call ed_Hij_add_link(2,2,1,1,1,-one*eloc)
  !> hoppings
  call ed_Hij_add_link(1,2,1,1,1,one*ts)
  call ed_Hij_add_link(2,1,1,1,1,one*ts)
  !
  call ed_Hij_info()
  call ed_Hij_write()


  call ed_init_solver()
  call ed_solve()



end program hubbard_dimer






