program dmrg_spin_1d
  USE SCIFOR
  USE DMRG
#ifdef _MPI
  USE MPI
#endif
  implicit none
  character(len=64)                  :: finput
  integer                            :: i,j,SUN,Unit,pos,Nsites
  real(8)                            :: Hvec,Noise,R,Sij
  type(site)                         :: MyDot
  type(sparse_matrix)                :: Sz,Sz2
  real(8),dimension(:,:),allocatable :: Hlr
  integer                            :: irank,comm,rank,ierr
  logical                            :: master=.true.,irun,imeasure
  character(len=:),allocatable       :: key_Sz
  
#ifdef _MPI  
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
#endif

  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(imeasure,"imeasure",finput,default=.true.,&
       comment="Bool to perform measurements. T for post-processing.")
  call parse_input_variable(irun,"irun",finput,default=.true.,&
       comment="Bool to run DMRG. F for post-processing")       
  call parse_input_variable(SUN,"SUN",finput,default=2,&
       comment="Spin SU(N) value. 2=> spin 1/2, 3=> spin 1")
  call parse_input_variable(Noise,"NOISE",finput,default=0d0,&
       comment="Magnetic field noise amplitude")
  call parse_input_variable(Hvec,"Hvec",finput,default=0d0,&
       comment="Magnetic field direction")

  call read_input(finput)

  if(Imeasure)then
     save_block=.true.
     save_umat=.true.
  endif

  Nsites=2*Ldmrg

  MyDot = spin_site(sun=SUN,Hz=Hvec)
  Hlr   = diag([Jp,Jx/2d0])

  !Init DMRG
  call init_dmrg(Hlr,ModelDot=[MyDot])


  !Run DMRG algorithm
  if(Irun)call run_DMRG()



  if(imeasure)then
     !Post-processing and measure quantities:
     !Measure <Sz(i)>, <Sz(i).Sz(i)>
     key_Sz="S"//mydot%okey(0,1,ilink="n")
     Sz =MyDot%operators%op(key_Sz)
     Sz2=matmul(Sz,Sz)
     !
     call Measure_DMRG([Sz,Sz2],file="Sz_Sz2VSj",pos=arange(1,Nsites))

     !Nearest-neighbour reference correlations: i, <S_i.S_(i+1)>.
     if(master)unit=fopen("spin_nnVSj"//str(label_DMRG('u')),append=.false.)
     do i=1,Nsites-1
        Sij=Measure_SpinSpin_DMRG(i,i+1)
        if(master)write(unit,*)i,Sij
     enddo
     if(master)close(unit)

     !Long-range reference correlations: j, <S_1.S_j>.
     if(master)unit=fopen("spin_1jVSj"//str(label_DMRG('u')),append=.false.)     
     do j=1,Nsites
        Sij=Measure_SpinSpin_DMRG(1,j)
        if(master)write(unit,*)j,Sij
     enddo
     if(master)close(unit)
     call End_Measure_DMRG()
     call Sz%free()
     call Sz2%free()
  endif

  !Finalize DMRG
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif

end program dmrg_spin_1d



