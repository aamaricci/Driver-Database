program anderson_hubbard_1d
  USE SCIFOR
  USE DMRG
#ifdef _MPI
  USE MPI
#endif
  implicit none

  integer                                        :: Nso
  character(len=64)                              :: finput
  integer                                        :: i,unit,iorb,ispin,idum
  real(8)                                        :: ts,Wdis
  type(site),dimension(:),allocatable            :: MyDot
  type(sparse_matrix),dimension(:,:),allocatable :: N,C
  type(sparse_matrix),dimension(:),allocatable   :: dens,docc,sz,s2z,Mvec
  real(8),dimension(:,:),allocatable             :: Hloc,Hlr
  real(8),dimension(:),allocatable               :: Erandom
  integer                                        :: irank,comm,rank,ierr
  logical                                        :: master,imeasure,bool

#ifdef _MPI  
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
#endif



  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"TS",finput,default=-1d0,&
       comment="Hopping amplitude")
  call parse_input_variable(Wdis,"WDIS",finput,default=0d0,&
       comment="local disorder amplitude")
  call parse_input_variable(Idum,"IDUM",finput,default=123456,&
       comment="disorder MT-seed")
  call parse_input_variable(imeasure,"imeasure",finput,default=.true.,&
       comment="Bool to include measurements directly in the run.")
  call read_input(finput)


  if(Norb>1)stop "This code is for Norb==1. STOP"
  Nso = Nspin*Norb
  !
  ! How do you own disorder? Disorder
  ! Now somewhere between the sacred silence
  ! Sacred silence and sleep
  ! Somewhere between the sacred silence and sleep
  ! Disorder, disorder, disorder
  allocate(Erandom(iNlat))
  call mersenne_init(idum)
  call mt_random(erandom)
  erandom=(2d0*erandom-1d0)*Wdis/2d0
  !
  inquire(file='erandom_'//str(idum)//'.restart',exist=bool)
  if(bool)then
     if(file_length('erandom_'//str(idum)//'.restart')/=iNlat)&
          stop "setup_disorder error: size(erandom_"//str(idum)//".restart) != iNlat"
     call read_array('erandom_'//str(idum)//'.restart',erandom)
  endif
  if(master)call save_array('erandom_'//str(idum)//'.used',erandom)



  allocate(Hloc(Nso,Nso))
  allocate(Hlr(Nso,Nso))
  !
  allocate(MyDot(iNlat))
  do i=1,iNlat
     Hloc     = Erandom(i)*pauli_sigma_0 !charge-noise
     MyDot(i) = electron_site(Hloc)
  enddo
  Hlr = ts*pauli_sigma_0          !this can be disordered too
  !
  call init_dmrg(Hlr,ModelDot=MyDot)


  !Run DMRG algorithm
  call run_DMRG()


  if(imeasure)then
     !Post-processing and measure quantities:
     !Measure <Sz(i)>
     allocate(C(Norb,Nspin),N(Norb,Nspin))
     do ispin=1,Nspin
        do iorb=1,Norb
           C(iorb,ispin) = myDot(1)%operators%op(key="C"//myDot(1)%okey(iorb,ispin))
           N(iorb,ispin) = matmul(C(iorb,ispin)%dgr(),C(iorb,ispin))
        enddo
     enddo
     allocate(Mvec(3*Norb),sz(Norb))
     do iorb=1,Norb
        sz(iorb)          = n(iorb,1)-n(iorb,2)
        Mvec(iorb)        = n(iorb,1)+n(iorb,2)
        Mvec(iorb+Norb)   = matmul(n(iorb,1),n(iorb,2))
        Mvec(iorb+2*Norb) = matmul(sz(iorb),sz(iorb))
     enddo
     !
     call Measure_DMRG(Mvec,file="n_d_s2z_l1VSj", pos=arange(1,Ldmrg))
  endif


  !Finalize DMRG
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif


end program anderson_hubbard_1d







