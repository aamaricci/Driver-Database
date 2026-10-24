program BHZ_1d
  USE SF_SPIN
  USE SF_CONSTANTS
  USE SF_ARRAYS
  USE SF_LINALG
  USE SF_IOTOOLS
  USE DMRG
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none

  integer                                        :: Nso
  character(len=64)                              :: finput
  integer                                        :: i,j,unit,iorb,ispin,Nsites,N3d
  real(8)                                        :: eh,mh,lambda
  real(8)                                        :: Eloc,Etot,Product,Ekin,E0loc,Eint,Eshift
  type(site)                                     :: Dot
  complex(8),dimension(4,4)                      :: GammaZ,GammaX
  complex(8),dimension(:,:),allocatable          :: Hloc,Hlr
  type(sparse_matrix),dimension(:,:),allocatable :: N,C
  type(sparse_matrix),dimension(:),allocatable   :: dens,docc,sz,s2z,Mvec
  type(sparse_matrix)                            :: Tz,Tx,O4f,H0loc,Hint,Hshift
  complex(8)                                     :: Corr
  integer                                        :: irank,comm,rank,ierr
  logical                                        :: master,imeasure,irun,ienergy,getIJcorr,getAVcorr,get1jcorr
  real(8),dimension(:,:),allocatable             :: dqs

#ifdef _MPI  
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
#endif


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(irun,"irun",finput,default=.true.,comment="Bool to run DMRG. F for post-processing")
  call parse_input_variable(imeasure,"imeasure",finput,default=.true.,comment="Bool to perform measurements. T for post-processing.")
  call parse_input_variable(ienergy,"ienergy",finput,default=.true.,comment="Bool to get energy measures. T for post-processing.")
  call parse_input_variable(mh,"MH",finput,default=0.5d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call read_input(finput)


  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"

  if(Imeasure.OR.Ienergy)then
    save_block=.true.
    save_umat=.true.
  endif       

  Nsites=2*Ldmrg
  Nso=Nspin*Norb
  !
  gammaX=kron( pauli_sigma_z, pauli_tau_x)
  gammaZ=kron( pauli_sigma_0, pauli_tau_z)

  !>Local Hamiltonian:
  allocate(Hloc(Nso,Nso))
  Hloc = Mh*GammaZ
  Dot  = electron_site(Hloc,H0loc=H0loc,Hint=Hint,Hshift=Hshift)

  !>Hopping Hamiltonian (i->i+1, right hop direction)
  if(allocated(Hlr))deallocate(Hlr)
  allocate(Hlr(Nso,Nso))
  Hlr = -0.5d0*GammaZ + 0.5d0*xi*lambda*GammaX
  if(master)then
    call print_matrix(Hloc,"Hloc.dmrg")
    call print_matrix(Hlr,"Hlr.dmrg")
  endif  

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
     allocate(Mvec(4*Norb),sz(Norb))
     do iorb=1,Norb
        sz(iorb)          = n(iorb,1)-n(iorb,2)
        Mvec(iorb)        = n(iorb,1)
        Mvec(iorb+Norb)   = n(iorb,2)
        Mvec(iorb+2*Norb) = matmul(n(iorb,1),n(iorb,2))
        Mvec(iorb+3*Norb) = matmul(sz(iorb),sz(iorb))
     enddo
     Tz = (n(2,1)+n(2,2)-n(1,1)-n(1,2))/2d0
     !
     !Measure <Tz>,<N>,<D>,<Sz^2> (write to file integrated)
     call Measure_DMRG([Tz,Mvec],file="tz_nup_ndw_d_s2zVSj", pos=arange(1,Nsites))
   
     do ispin=1,Nspin
      do iorb=1,Norb
        call C(iorb,ispin)%free()
        call N(iorb,ispin)%free()
      enddo
     enddo
    !
     call End_Measure_DMRG()
  endif

  if(ienergy)then
    !Measure energies: <K>,<Hloc>
    call Init_Measure_DMRG()
    if(master)print*,"measure energies"
    if(master)unit=fopen("Ekin_Eloc_Etot"//str(label_DMRG('u')),append=.false.)
    call Measure_Energy_DMRG(Hlr,Ekin,Eloc,Etot,&
       H0loc=H0loc,Hint=Hint,Hshift=Hshift,&
       E0loc=E0loc,Eint=Eint,Eshift=Eshift)
    if(master)write(unit,*)Ekin,Eloc,Etot,E0loc,Eint,Eshift
    if(master)close(unit)
    call End_Measure_DMRG()
  endif




  !Finalize DMRG
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif



contains


  subroutine get_correlations(label,OpA,dqA,OpB,dqB)
    character(len=*)                   :: label
    type(sparse_matrix)                :: OpA,OpB
    real(8),dimension(:)               :: dqA,dqB
    integer                            :: i,j,r,ic,count,cut
    real(8),dimension(:,:),allocatable :: Cij
    real(8),dimension(:),allocatable   :: Cavg
    real(8)                            :: sum_corr
    !
    allocate(Cij(Nsites,Nsites))
    Cij=zero
    !
    !Plot <O_i.O_j> with |j-i| as Manhattan distance
    if(master) unit = fopen(str(label)//"_r"//str(label_DMRG('u')), append=.false.)
    if(master)call start_timer("<Oi.Oj>_|i-j|")
    ic = Nsites / 2      
    do r = 1, ic - 1
      i = ic - r / 2
      j = ic + (r + 1) / 2
      if(Cij(i,j)==zero)Cij(i,j) = dreal(Measure_Corr_DMRG(OpA,dqA,OpB,dqB,i,j,connected=.true.))
      if(master) write(unit,*) abs(i-j), Cij(i,j)
      if(master)call eta(r,ic)
    enddo
    if(master)call stop_timer()
    if(master)close(unit)
    ! 
  end subroutine get_correlations



end program BHZ_1d



