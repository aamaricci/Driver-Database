program Ce2O2FeSe2_1d
  USE SCIFOR
  USE DMRG
#ifdef _MPI
  USE MPI
#endif
  implicit none

integer                                        :: Nso
character(len=64)                              :: finput
integer                                        :: i,unit,iorb,ispin
type(site)                                     :: Dot
real(8),dimension(:,:),allocatable             :: Hloc,Hlr,T
type(sparse_matrix),dimension(:,:),allocatable :: N,C
type(sparse_matrix),dimension(:),allocatable   :: dens,docc,sz,s2z,Mvec
integer                                        :: irank,comm,rank,ierr
logical                                        :: master,imeasure,irun


#ifdef _MPI  
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
#endif


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(irun,"irun",finput,default=.true.,&
       comment="Bool to run DMRG. F for post-processing")
  call parse_input_variable(imeasure,"imeasure",finput,default=.true.,&
       comment="Bool to perform measurements. T for post-processing.")
  call read_input(finput)



  if(Nspin/=2.OR.Norb/=3)&
       stop "Wrong setup from input file. We require Nspin=2; Norb=3"
  Nso=Nspin*Norb


  !>Local Hamiltonian:
  allocate(Hloc(Nso,Nso))
  Hloc = eye(Nspin).kx.diag([-0.399538d0,-0.325538d0,-0.842538d0])
  Dot  = electron_site(Hloc)

  !>Hopping Hamiltonian (i->i+1, right hop direction)

  allocate(T(Norb,Norb));T   = zero
  T(1,1) =  0.187000 
  T(2,1) =  0.054000
  T(3,1) =  0.020000
  T(1,2) = -0.054000
  T(2,2) =  0.351000
  T(3,2) =  0.349000
  T(1,3) =  0.020000
  T(2,3) = -0.349000
  T(3,3) = -0.433000
  allocate(Hlr(Nso,Nso));Hlr = zero
  Hlr = eye(Nspin).kx.T

  !Init DMRG
  call init_dmrg(Hlr,ModelDot=[Dot])

  !Run DMRG algorithm
  if(Irun)call run_DMRG()

  if(Imeasure)then
     !Post-processing and measure quantities:
     allocate(C(Norb,Nspin),N(Norb,Nspin))
     do ispin=1,Nspin
        do iorb=1,Norb
           C(iorb,ispin) = dot%operators%op(key="C"//dot%okey(iorb,ispin))
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




end program Ce2O2FeSe2_1d



