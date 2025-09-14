program dmrg_spin_1d
  USE SCIFOR
  USE DMRG
#ifdef _MPI
  USE MPI
#endif
  implicit none
  character(len=64)                   :: finput
  integer                             :: i,SUN,Unit,pos
  real(8)                             :: Hvec(3),Noise,R,Sij
  type(site),dimension(:),allocatable :: MyDot
  type(sparse_matrix)                 :: bSz,bSp,SiSj
  real(8),dimension(:,:),allocatable  :: Hlr
  integer                             :: irank,comm,rank,ierr
  logical                             :: master

#ifdef _MPI  
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
#endif
  
  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(SUN,"SUN",finput,default=2,&
       comment="Spin SU(N) value. 2=> spin 1/2, 3=> spin 1")
  call parse_input_variable(Noise,"NOISE",finput,default=0d0,&
       comment="Magnetic field noise amplitude")
  call parse_input_variable(Hvec,"Hvec",finput,default=[0d0,0d0,0d0],&
       comment="Magnetic field direction")
  call read_input(finput)


  !Init DMRG

  ! allocate(MyDot(2*Ldmrg))
  ! call mersenne_init(12345)
  ! do i=1,2*Ldmrg
  !    R      = Noise*mersenne()
  !    MyDot(i) = spin_site(sun=SUN,Hvec=R*Hvec)
  ! enddo
  allocate(MyDot(1))
  MyDot(1) = spin_site(sun=SUN,Hvec=Hvec)

  ! if(allocated(Hlr))deallocate(Hlr)
  ! allocate(Hlr(Nspin*Norb,Nspin*Norb))
  Hlr = diag([Jp,Jx/2d0])
  call init_dmrg(Hlr,ModelDot=MyDot)


  !Run DMRG algorithm
  call run_DMRG()


  !Post-processing and measure quantities:
  !Measure <Sz(i)>
  call Measure_DMRG(MyDot(1)%operators%op(key="S_z"),file="SzVSj")

  


  !Measure <S(i).S(i+1)>
  if(master)unit=fopen("SiSjVSsite"//str(label_DMRG('u')),append=.true.)
  call Init_measure_dmrg("SiSjVSsite")
  do pos=1,Ldmrg-1
     bSz = Build_Op_DMRG(MyDot(1)%operators%op("S_z"),pos,set_basis=.true.)
     bSp = Build_Op_DMRG(MyDot(1)%operators%op("S_p"),pos,set_basis=.true.)
     SiSj= get_SiSj(bSz,bSp,MyDot(1)%operators%op("S_z"),MyDot(1)%operators%op("S_p"))
     SiSj= Advance_Corr_DMRG(SiSj,pos)
     Sij = Average_Op_DMRG(SiSj,pos)
     if(master)write(unit,*)pos,Sij
     if(master)call eta(pos,Ldmrg-1)
  enddo
  call End_measure_dmrg()
  if(Master)close(unit)


  !Finalize DMRG
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif
  
contains


  function get_SiSj(Sz1,Sp1,Sz2,Sp2) result(sisj)
    type(sparse_matrix) :: sisj
    type(sparse_matrix) :: Sz1,Sp1,Sz2,Sp2
    SiSj = 0.5d0*(Sp1.x.Sp2%dgr()) +  0.5d0*(Sp1%dgr().x.Sp2)  + (Sz1.x.Sz2)
  end function get_SiSj


end program dmrg_spin_1d





