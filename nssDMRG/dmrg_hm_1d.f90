program hubbard_1d
  USE SCIFOR
  USE DMRG
#ifdef _MPI
  USE MPI
#endif
  implicit none

  integer                                        :: Nso
  character(len=64)                              :: finput
  integer                                        :: i,j,unit,iorb,ispin,pos,Nsites
  real(8)                                        :: ts,Mh,lambda,val,K,alpha
  real(8)                                        :: Eloc,Etotal
  type(site)                                     :: MyDot
  type(sparse_matrix)                            :: Kij,Hi
  type(sparse_matrix),dimension(:,:),allocatable :: N,C
  type(sparse_matrix),dimension(:),allocatable   :: pair,dens,docc,sz,s2z,Mvec
  real(8),dimension(:,:),allocatable             :: Hloc,Hlr,avLocal,corr
  real(8),dimension(:),allocatable               :: corr_values
  integer                                        :: irank,comm,rank,ierr
  logical                                        :: master=.true.,imeasure,irun
  
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
  call parse_input_variable(ts,"TS",finput,default=-0.5d0,comment="Hopping amplitude")
  call parse_input_variable(alpha,"alpha",finput,default=1d0,comment="bandwidth ratio")
  call parse_input_variable(Mh,"MH",finput,default=0d0,comment="Crystal field splittings")
  call parse_input_variable(lambda,"LAMBDA",finput,default=0d0,comment="off-diagonal amplitude")

  call read_input(finput)
  
  if(Nspin/=2)stop "Use of this code requires Nspin=2. STOP"


  if(Imeasure)then
     save_block=.true.
     save_umat=.true.
  endif

  Nsites = 2*Ldmrg
  Nso    = Nspin*Norb
  allocate(Hloc(Nso,Nso))
  allocate(Hlr(Nso,Nso))
  select case(Norb)
  case(1)
     Hloc = Mh*pauli_z          !use it as a Zeeman field
     Hlr  = ts*pauli_0
  case(2)                       !spin x orbital <= ext x int
     Hloc = one*Mh*kron(pauli_0,pauli_z)
     Hlr  = one*ts*kron(pauli_0,diag([1d0,alpha])) + one*lambda/2d0*kron(pauli_0,pauli_x)
  case default;stop "This code is for Norb<=2. STOP"
  end select
  if(master)then
     call print_matrix(Hloc,"Hloc.dmrg")
     call print_matrix(Hlr,"Hlr.dmrg")
  endif

  !Setup Dot basis:
  MyDot = electron_site(Hloc)

  !Init DMRG
  call init_dmrg(Hlr,ModelDot=[MyDot])


  !Run DMRG algorithm
  if(Irun)call run_DMRG()


  if(imeasure)then      
     !Post-processing and measure quantities:
     allocate(C(Norb,Nspin),N(Norb,Nspin))
     do ispin=1,Nspin
        do iorb=1,Norb
           C(iorb,ispin) = myDot%operators%op(key="C"//myDot%okey(iorb,ispin,ilink='n'))
           N(iorb,ispin) = matmul(C(iorb,ispin)%dgr(),C(iorb,ispin))
        enddo
     enddo
     allocate(Mvec(3*Norb),sz(Norb))
     do iorb=1,Norb
        sz(iorb)          = 0.5d0*(n(iorb,1)-n(iorb,2))
        Mvec(iorb)        = n(iorb,1)+n(iorb,2)
        Mvec(iorb+Norb)   = matmul(n(iorb,1),n(iorb,2))
        Mvec(iorb+2*Norb) = matmul(sz(iorb),sz(iorb))
     enddo
     !
     !
     call Measure_DMRG(Mvec,file="n_d_s2zVSj", pos=arange(1,Nsites))
     !
     !Measure <K>
     if(master)unit=fopen("K"//str(label_DMRG('u')),append=.true.)
     call Measure_Energy_DMRG(Hlr,K,Eloc,Etotal,Kij,Hi)
     if(master)write(unit,*)K,Eloc,Etotal
     if(master)close(unit)
     !
     !The density matrix is flattened with io outermost and jo innermost.
     if(master)print*,"measure density.density nn"
     allocate(corr_values(Nso*Nso))
     if(master)unit=fopen("nn_nnVSj"//str(label_DMRG('u')),append=.false.)
     do i=1,Nsites-1
         corr=Measure_DensityDensity_DMRG(i,i+1)
         if(master)write(unit,*)i,flatten_correlation(corr)
     enddo
     if(master)close(unit)
 
     if(master)print*,"measure density.density 1n"
     if(master)unit=fopen("nn_1jVSj"//str(label_DMRG('u')),append=.false.)  
     do j=1,Nsites-1
         corr=Measure_DensityDensity_DMRG(1,j)
         if(master)write(unit,*)j,flatten_correlation(corr)
     enddo
     if(master)close(unit)
     !
     call End_measure_dmrg()
     do ispin=1,Nspin
      do iorb=1,Norb
        call C(iorb,ispin)%free()
        call N(iorb,ispin)%free()
      enddo
    enddo
  endif


  !Finalize DMRG
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif

  contains

    function flatten_correlation(Cij) result(values)
      real(8),intent(in) :: Cij(:,:)
      real(8)            :: values(size(Cij))
      integer            :: io,jo,k
      k=0
      do io=1,size(Cij,1)
         do jo=1,size(Cij,2)
            k=k+1
            values(k)=Cij(io,jo)
         enddo
      enddo
    end function flatten_correlation


end program hubbard_1d




