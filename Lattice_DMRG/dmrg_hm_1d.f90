program hubbard_1d
  USE SCIFOR
  USE DMRG
#ifdef _MPI
  USE MPI
#endif
  implicit none

  integer                                        :: Nso
  character(len=64)                              :: finput
  integer                                        :: i,unit,iorb,ispin,pos
  real(8)                                        :: ts(2),Mh,lambda,val
  type(site)                                     :: MyDot
  type(sparse_matrix) :: P,Cl,Pl,Tij
  type(sparse_matrix),dimension(:,:),allocatable :: N,C
  type(sparse_matrix),dimension(:),allocatable   :: pair,dens,docc,sz,s2z,Mvec
  real(8),dimension(:,:),allocatable             :: Hloc,Hlr
  integer                                        :: irank,comm,rank,ierr
  logical                                        :: master,imeasure

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
  call parse_input_variable(ts,"TS",finput,default=(/( -0.5d0,i=1,2 )/),&
       comment="Hopping amplitudes")
  call parse_input_variable(Mh,"MH",finput,default=0d0,&
       comment="Crystal field splittings")
  call parse_input_variable(lambda,"LAMBDA",finput,default=0d0,&
       comment="off-diagonal amplitude")

  call read_input(finput)


  Nso = Nspin*Norb
  allocate(Hloc(Nso,Nso))
  allocate(Hlr(Nso,Nso))
  select case(Norb)
  case(1)
     Hloc = Mh*pauli_z          !use it as a Zeeman field
     Hlr  = diag([ts(1:Norb),ts(1:Norb)])
  case(2)
     Hloc = Mh*kron(pauli_0,pauli_z)
     Hlr  = diag([ts(1:Norb),ts(1:Norb)]) &
          + lambda*kron(pauli_0,pauli_x)
  case default;stop "This code is for Norb<=2. STOP"
  end select
  if(master)then
     call print_matrix(Hloc,"Hloc.dmrg")
     call print_matrix(Hlr,"Hlr.dmrg")
  endif

  !Setup Dot basis:
  MyDot = electron_site()

  !Init DMRG
  call init_dmrg(Hlr,ModelDot=[MyDot])


  !Run DMRG algorithm
  call run_DMRG()


  if(imeasure)then
     !Post-processing and measure quantities:
     allocate(C(Norb,Nspin),N(Norb,Nspin))
     P=myDot%operators%op(key="P")
     do ispin=1,Nspin
        do iorb=1,Norb
           C(iorb,ispin) = myDot%operators%op(key="C"//myDot%okey(iorb,ispin))
           N(iorb,ispin) = matmul(C(iorb,ispin)%dgr(),C(iorb,ispin))
        enddo
     enddo
     allocate(Mvec(4*Norb),sz(Norb))
     do iorb=1,Norb
        sz(iorb)          = n(iorb,1)-n(iorb,2)
        Mvec(iorb)        = n(iorb,1)+n(iorb,2)
        Mvec(iorb+Norb)   = matmul(n(iorb,1),n(iorb,2))
        Mvec(iorb+2*Norb) = matmul(sz(iorb),sz(iorb))
        Mvec(iorb+3*Norb) = matmul(C(iorb,1),C(iorb,2))
     enddo
     !
     call Measure_DMRG(Mvec,file="n_d_s2z_pairVSj", pos=arange(1,Ldmrg))



     ! !Measure <S(i).S(i+1)>
     ! if(master)unit=fopen("PijVSsite"//str(label_DMRG('u')),append=.true.)
     ! call Init_measure_dmrg("PijVSsite")
     ! do pos=1,Ldmrg-1
     !    Pl  = Build_Op_DMRG(P,pos,set_basis=.true.)
     !    Cl  = Build_Op_DMRG(C(1,1),pos,set_basis=.true.)
     !    Tij = get_Tij()
     !    Tij = Advance_Corr_DMRG(Tij,pos)
     !    val = Average_Op_DMRG(Tij,pos)
     !    if(master)write(unit,*)pos,val
     !    if(master)call eta(pos,Ldmrg-1)
     ! enddo
     ! call End_measure_dmrg()
     ! if(Master)close(unit)

  endif


  !Finalize DMRG
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif

contains

  function get_Tij() result(Tij)
    type(sparse_matrix) :: Tij
    Tij = matmul(Cl,Pl).x.C(1,2) 
  end function get_Tij


end program hubbard_1d







