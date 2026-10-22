! Post-processing driver for the interacting 1D BHZ model.
! Use with irun=.false. in DMRG.conf: the ground state is read from restart.
program wick_1d_post
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
  integer                                        :: Nso,Nsites,iorb,ispin,unit,comm,rank
  character(len=64)                              :: finput
  real(8)                                        :: Eloc,Etot,Ekin,E0loc,Eint,Eshift
  real(8)                                        :: mh,lambda
  type(site)                                     :: Dot
  complex(8),dimension(4,4)                      :: GammaZ,GammaX
  complex(8),dimension(:,:),allocatable          :: Hloc,Hlr
  type(sparse_matrix),dimension(:,:),allocatable :: C,N
  type(sparse_matrix)                            :: Sz,Tx,Ty,Tz,Tx2,Ty2,Tz2,Hint,H0loc,Hshift
  logical                                        :: master
#ifdef _MPI
  call init_MPI()
  comm=MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank=get_Rank_MPI(comm) 
  master=get_Master_MPI(comm)
#else
  master=.true.
#endif
  call parse_cmd_variable(finput,"FINPUT",default='used.DMRG.conf')
  call read_input(finput)
  !
  if(Nspin/=2 .or. Norb/=2) stop "BHZ post driver requires Nspin=Norb=2"
  Nsites=2*Ldmrg
  Nso=Nspin*Norb
  !
  GammaX=kron(pauli_sigma_z,pauli_tau_x)
  GammaZ=kron(pauli_sigma_0,pauli_tau_z)
  !
  allocate(Hloc(Nso,Nso),Hlr(Nso,Nso)); Hloc=mh*GammaZ
  Dot=electron_site(Hloc,H0loc=H0loc,Hint=Hint,Hshift=Hshift)
  !
  Hlr=-0.5d0*GammaZ+0.5d0*xi*lambda*GammaX
  !
  call init_dmrg(Hlr,ModelDot=[Dot])
  !
  allocate(C(Norb,Nspin),N(Norb,Nspin))
  do ispin=1,Nspin; do iorb=1,Norb
    C(iorb,ispin)=Dot%operators%op(key="C"//Dot%okey(iorb,ispin))
    N(iorb,ispin)=matmul(C(iorb,ispin)%dgr(),C(iorb,ispin))
  enddo; enddo
  Sz= N(1,1)+N(2,1)-N(1,2)-N(2,2)
  Tz=(N(2,1)+N(2,2)-N(1,1)-N(1,2))/2d0
  !
  Tx=0.5d0*(matmul(C(1,1)%dgr(),C(2,1))+matmul(C(2,1)%dgr(),C(1,1)))
  Tx=Tx+0.5d0*(matmul(C(1,2)%dgr(),C(2,2))+matmul(C(2,2)%dgr(),C(1,2)))
  !
  Ty=cmplx(0d0,-0.5d0,8)*(matmul(C(1,1)%dgr(),C(2,1))-matmul(C(2,1)%dgr(),C(1,1)))
  Ty=Ty+cmplx(0d0,-0.5d0,8)*(matmul(C(1,2)%dgr(),C(2,2))-matmul(C(2,2)%dgr(),C(1,2)))
  !
  Tx2=matmul(Tx,Tx)
  Ty2=matmul(Ty,Ty)
  Tz2=matmul(Tz,Tz)
  !
  call get_wick_residuals(C)
  !
  !Free the operators:
  do ispin=1,Nspin
    do iorb=1,Norb
      call C(iorb,ispin)%free(); 
      call N(iorb,ispin)%free()
    enddo 
  enddo
  !
  call End_Measure_DMRG()
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif


contains


  ! For a number-conserving Gaussian state Wick's theorem gives
  ! <c1^dag c2^dag c2 c1> = n1*n2
  !   - <c1^dag c2><c2^dag c1>.
  ! The residual below is zero for a Gaussian state.  We sample local
  ! inter-orbital density channels; a complete tensor scales as O(L^4).
  subroutine get_wick_residuals(Cop)
    type(sparse_matrix),dimension(:,:),intent(in) :: Cop
    integer                                       :: i,unit,ispin,jspin
    real(8)                                       :: exact,wick,residual,n1,n2,coh12,coh21
    real(8),dimension(2,4)                        :: dqs
    character(len=16),dimension(4)                :: kinds
    dqs=0d0;
    kinds=["none","none","none","none"]
    if(master) unit=fopen("wick_residual"//str(label_DMRG('u')),append=.false.)
    if(master) write(unit,*)"# site sigma sigma_prime exact wick residual"
    do i=1,Nsites; 
      do ispin=1,Nspin; 
        do jspin=1,Nspin
            exact=dreal(Measure_Product_DMRG([Cop(1,ispin)%dgr(),Cop(2,jspin)%dgr(),&
            Cop(2,jspin),Cop(1,ispin)],dqs,kinds,[i,i,i,i]))
            n1=dreal(Measure_Corr_DMRG(Cop(1,ispin)%dgr(),[0d0,0d0],Cop(1,ispin),[0d0,0d0],i,i,connected=.false.))
            n2=dreal(Measure_Corr_DMRG(Cop(2,jspin)%dgr(),[0d0,0d0],Cop(2,jspin),[0d0,0d0],i,i,connected=.false.))
            coh12=dreal(Measure_Corr_DMRG(Cop(1,ispin)%dgr(),[0d0,0d0],Cop(2,jspin),[0d0,0d0],i,i,connected=.false.))
            coh21=dreal(Measure_Corr_DMRG(Cop(2,jspin)%dgr(),[0d0,0d0],Cop(1,ispin),[0d0,0d0],i,i,connected=.false.))
            wick=n1*n2-coh12*coh21
            !
            residual=exact-wick
            !
            if(master) write(unit,*)i,ispin,jspin,exact,wick,residual
        enddo; 
      enddo; 
    enddo
    if(master) close(unit)
  end subroutine get_wick_residuals










end program wick_1d_post
