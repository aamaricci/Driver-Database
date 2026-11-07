! Post-processing driver for the interacting 1D BHZ model.
! Use with irun=.false. in DMRG.conf: the ground state is read from restart.
program BHZ_1d_post
  USE SF_SPIN
  USE SF_CONSTANTS
  USE SF_ARRAYS
  USE SF_LINALG
  USE SF_IOTOOLS
  USE SF_TIMER
  USE DMRG
  USE RESTART_DIR, only: prepare_restart,cleanup_restart
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif

  implicit none
  integer                                        :: Nso,Nsites,iorb,ispin,unit,comm,rank
  character(len=64)                              :: finput
  real(8)                                        :: Eloc,Etot,Ekin,E0loc,Eint,Eshift
  real(8)                                        :: mh,lambda
  type(site)                                     :: Dot
  complex(8),dimension(4,4)                      :: GammaZ,GammaX
  complex(8),dimension(:,:),allocatable          :: Hloc,Hlr
  type(sparse_matrix),dimension(:,:),allocatable :: C,N
  type(sparse_matrix)                            :: Sz
  type(sparse_matrix)                            :: Splus,Eplus
  real(8),allocatable                            :: mz(:)
  integer                                        :: corr_rmax
  logical                                        :: corr_full_distance
  logical                                        :: master=.true.
#ifdef _MPI
  call init_MPI()
  comm=MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank=get_Rank_MPI(comm)
  master=get_Master_MPI(comm)
#endif
  call parse_cmd_variable(finput,"FINPUT",default='used.DMRG.conf')
  call parse_input_variable(mh,"MH",finput,default=0.5d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_input_variable(corr_full_distance,"CORR_FULL_DISTANCE",finput,default=.false., &
       comment="Extend central-pair correlations to r=L-1 instead of L/2-1")
  call read_input(finput)
  !
  if(Nspin/=2 .or. Norb/=2) stop "BHZ post driver requires Nspin=Norb=2"
  if(QNdim/=2.or.to_lower(DMRG_MODE)/="normal") &
    error stop "Transverse correlators require normal mode and QNDIM=2"
  Nsites=2*Ldmrg
  corr_rmax=Nsites/2-1
  if(corr_full_distance)corr_rmax=Nsites-1
  Nso=Nspin*Norb
  !
  GammaX=kron(pauli_sigma_z,pauli_tau_x)
  GammaZ=kron(pauli_sigma_0,pauli_tau_z)
  !
  allocate(Hloc(Nso,Nso),Hlr(Nso,Nso))
  Hloc=mh*GammaZ
  Dot=electron_site(Hloc)
  !
  Hlr=-0.5d0*GammaZ+0.5d0*xi*lambda*GammaX
  !
  call prepare_restart()
  call init_dmrg(Hlr,ModelDot=[Dot])
  !
  allocate(C(Norb,Nspin),N(Norb,Nspin))
  do ispin=1,Nspin
    do iorb=1,Norb
      C(iorb,ispin)=Dot%operators%op(key="C"//Dot%okey(iorb,ispin))
      N(iorb,ispin)=matmul(C(iorb,ispin)%dgr(),C(iorb,ispin))
    enddo
  enddo
  Sz= N(1,1)+N(2,1)-N(1,2)-N(2,2)
  !
  if(master)print*,"<Sz_i.Sz_j>_c"
  call get_correlations("sz.sz",Sz)
  call get_Sz_disconnected_correlations()
  !
  ! Definite QN shifts: plus=(+1,-1), minus=(-1,+1); even fermion parity.
  Splus=matmul(C(1,1)%dgr(),C(1,2))+matmul(C(2,1)%dgr(),C(2,2))
  call get_transverse_correlations('sx.sx',Splus)
  !
  !Free the operators:
  do ispin=1,Nspin
    do iorb=1,Norb
      call C(iorb,ispin)%free()
      call N(iorb,ispin)%free()
    enddo
  enddo
  call Splus%free()
  !
  !
  call End_Measure_DMRG()
  call finalize_dmrg()
  call cleanup_restart()
#ifdef _MPI
  call finalize_MPI()
#endif
contains

  subroutine get_transverse_correlations(label,Plus)
    character(len=*),intent(in) :: label
    type(sparse_matrix),intent(in) :: Plus
    type(sparse_matrix) :: Minus
    complex(8) :: value
    integer :: i,j,r,ic,iu
    Minus=Plus%dgr()
    ic=Nsites/2
    if(master)then
      iu=fopen(str(label)//"_r"//str(label_DMRG('u')),append=.false.)
      write(iu,*) '# r, Re Cxx, Im Cxx, i; j=i+r; Cxx=(C+- + C-+)/4'
    endif
    do r=1,corr_rmax
      i=ic-r/2
      j=ic+(r+1)/2
      value=0.25d0*( &
        Measure_Corr_DMRG(Plus,[1d0,-1d0],Minus,[-1d0,1d0],i,j)+ &
        Measure_Corr_DMRG(Minus,[-1d0,1d0],Plus,[1d0,-1d0],i,j))
      if(master)write(iu,*) r,dreal(value),dimag(value),i
    enddo
    if(master)close(iu)
    call Minus%free()
  end subroutine get_transverse_correlations

  subroutine get_Sz_disconnected_correlations()
    integer :: i,j,r,ic,iu
    ! Compute the profile from the same checkpoint; no dependency on average output.
    allocate(mz(Nsites))
    do i=1,Nsites
      mz(i)=Measure_Op_DMRG(Sz,i)
    enddo
    ic=Nsites/2
    if(master)then
      iu=fopen('sz.sz_disconnected_r'//str(label_DMRG('u')),append=.false.)
      write(iu,*) '# r, Dzz(i,r), i; Dzz=mz(i)*mz(i+r); mz=n_up-n_down'
      do r=1,corr_rmax
        i=ic-r/2
        j=ic+(r+1)/2
        write(iu,*) r,mz(i)*mz(j),i
      enddo
      close(iu)
    endif
    deallocate(mz)
  end subroutine get_Sz_disconnected_correlations

  subroutine get_correlations(label,Op)
    character(len=*),intent(in)    :: label
    type(sparse_matrix),intent(in) :: Op
    integer                        :: i,j,r,ic
    real(8)                        :: value
    if(master) unit=fopen(str(label)//"_r"//str(label_DMRG('u')),append=.false.)
    ic=Nsites/2
    if(master)call start_timer("<"//str(label)//">_|i-j|")
    do r=1,corr_rmax
      i=ic-r/2
      j=ic+(r+1)/2
      value=dreal(Measure_Corr_DMRG(Op,[0d0,0d0],Op,[0d0,0d0],i,j,connected=.true.))
      if(master) write(unit,*) abs(i-j),value
      if(master) call eta(r,corr_rmax)
    enddo
    if(master)call stop_timer()
    if(master) close(unit)
  end subroutine get_correlations

end program BHZ_1d_post
