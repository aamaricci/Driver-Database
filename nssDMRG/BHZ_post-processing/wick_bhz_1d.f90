! Post-processing driver for the interacting 1D BHZ model.
! Reads used.DMRG.conf, or FINPUT explicitly.
! Prepares restart/ from restart_*.tgz, or copies an uncompressed restart_* directory;
! removes only that temporary directory after successful post-processing.
program wick_1d_post
  USE SF_SPIN
  USE SF_CONSTANTS
  USE SF_ARRAYS
  USE SF_LINALG
  USE SF_IOTOOLS
  USE DMRG
  use DMRG_GLOBAL, only: current_L
  USE RESTART_DIR, only: prepare_restart,cleanup_restart
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif

  implicit none
  integer                                        :: Nso,Nsites,iorb,ispin,unit,comm=0,rank
  character(len=256)                             :: finput
  real(8)                                        :: Eloc,Etot,Ekin,E0loc,Eint,Eshift
  real(8)                                        :: mh,lambda
  type(site)                                     :: Dot
  complex(8),dimension(4,4)                      :: GammaZ,GammaX
  complex(8),dimension(:,:),allocatable          :: Hloc,Hlr
  type(sparse_matrix),dimension(:,:),allocatable :: C,N
  type(sparse_matrix)                            :: Hint,H0loc,Hshift
  logical                                        :: master=.true.
#ifdef _MPI
  call init_MPI()
  comm=MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank=get_Rank_MPI(comm)
  master=get_Master_MPI(comm)
#else
  comm=0
  master=.true.
#endif
  call parse_cmd_variable(finput,"FINPUT",default='used.DMRG.conf')
  call parse_input_variable(mh,"MH",finput,default=0.5d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
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
  call prepare_restart()
  call init_dmrg(Hlr,ModelDot=[Dot])
  !
  allocate(C(Norb,Nspin),N(Norb,Nspin))
  do ispin=1,Nspin; do iorb=1,Norb
      C(iorb,ispin)=Dot%operators%op(key="C"//Dot%okey(iorb,ispin))
      N(iorb,ispin)=matmul(C(iorb,ispin)%dgr(),C(iorb,ispin))
    enddo; enddo
  !
  call Init_Measure_DMRG()
  !
  call get_wick_residuals(C)
  !
  call get_one_particle_density(C)
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
  call cleanup_restart()
#ifdef _MPI
  call finalize_MPI()
#endif


contains



  ! Full-chain one-particle matrix; spin off-diagonal blocks vanish in normal mode.
  ! R(mu,nu)=<c_mu^dagger c_nu>, mu=(i-1)*Norb+a within each spin block.
  subroutine get_one_particle_density(Cop)
    type(sparse_matrix),intent(in) :: Cop(:,:)
    complex(8),allocatable         :: R(:,:),work(:,:)
    real(8),allocatable            :: occupations(:),dq(:)
    real(8)                        :: entropy(Nspin),defect(Nspin),traces(Nspin),lo(Nspin),hi(Nspin)
    complex(8)                     :: realpart,imagpart
    real(8)                        :: n,diagonal_imag(Nspin),trace_error(Nspin),occupation_tol
    type(sparse_matrix)            :: creator
    character(len=:),allocatable   :: key
    integer                        :: L,K,s,mu,nu,i,j,a,b,iu,iv,iw
    !
    ! Endpoint clipping only repairs floating-point excursions outside [0,1].
    occupation_tol=1d-8
    if(QNdim/=2.or.to_lower(DMRG_MODE)/='normal') &
      error stop 'One-particle non-freeness requires normal mode and QNDIM=2'
    L=current_L
    K=L*Norb
    allocate(R(K,K),work(K,K),occupations(K))
    !
    entropy      = 0d0
    defect       = 0d0
    traces       = 0d0
    diagonal_imag= 0d0
    !
    if(master)then
      iu=fopen('one_particle_density'//str(label_DMRG('u')),append=.false.)
      iv=fopen('natural_occupations'//str(label_DMRG('u')),append=.false.)
      write(iu,*) '# spin i a j b Re R Im R'
      write(iv,*) '# spin alpha occupation (raw eigenvalue, ascending)'
    endif
    !
    do s=1,Nspin
      key='C'//Dot%okey(1,s)
      dq =Dot%operators%dq(key)
      !
      do mu=1,K
        i=(mu-1)/Norb+1
        a=mod(mu-1,Norb)+1
        creator=Cop(a,s)%dgr()
        do nu=mu,K
          j=(nu-1)/Norb+1
          b=mod(nu-1,Norb)+1
          !
          realpart=Measure_Corr_DMRG(creator,-dq,Cop(b,s),dq,i,j, &
            typeA='fermionic',typeB='fermionic',connected=.false.)
          ! Same-block Average_Op returns a real expectation in this library.
          ! Recover Im <A> as Re <-i A>, also valid for complex LR contractions.
          imagpart=Measure_Corr_DMRG(dcmplx(0d0,-1d0)*creator,-dq,Cop(b,s),dq,i,j, &
            typeA='fermionic',typeB='fermionic',connected=.false.)
          !
          R(mu,nu)=dcmplx(dreal(realpart),dreal(imagpart))
          if(mu==nu)diagonal_imag(s)=max(diagonal_imag(s),abs(dimag(R(mu,nu))))
          R(nu,mu)=conjg(R(mu,nu))
          !
        enddo
        call creator%free()
      enddo
      !
      ! Export both triangles before diagonalization overwrites the workspace.
      if(master)then
        do mu=1,K
          i=(mu-1)/Norb+1
          a=mod(mu-1,Norb)+1
          do nu=1,K
            j=(nu-1)/Norb+1
            b=mod(nu-1,Norb)+1
            write(iu,*) s,i,a,j,b,dreal(R(mu,nu)),dimag(R(mu,nu))
          enddo
          traces(s)=traces(s)+dreal(R(mu,mu))
          R(mu,mu) =dcmplx(dreal(R(mu,mu)),0d0)
        enddo
        !
        !Get natural occupations:
        if(diagonal_imag(s)>occupation_tol)error stop '1RDM: complex diagonal exceeds tolerance'
        work=R
        call eigh(work,occupations)
        lo(s)=minval(occupations)
        hi(s)=maxval(occupations)
        if(lo(s)<-occupation_tol.or.hi(s)>1d0+occupation_tol) &
          error stop '1RDM: occupation outside [0,1] beyond tolerance'
        trace_error(s)=abs(sum(occupations)-traces(s))
        !
        do mu=1,K
          write(iv,*) s,mu,occupations(mu)
          n=max(0d0,min(1d0,occupations(mu)))
          if(n>0d0)entropy(s)=entropy(s)-n*log(n)
          if(n<1d0)entropy(s)=entropy(s)-(1d0-n)*log(1d0-n)
          defect(s)=defect(s)+n*(1d0-n)
        enddo
      endif
      deallocate(dq)
    enddo
    !
    if(master)then
      close(iu);close(iv)
      iw=fopen('non_freeness'//str(label_DMRG('u')),append=.false.)
      write(iw,*) '# Pure full-chain state, natural logarithms (nats); L=',L
      write(iw,*) '# C_total C_per_site C_up C_down Tr[R(1-R)]_total Tr[R(1-R)]_per_site'
      write(iw,*) sum(entropy),sum(entropy)/L,entropy,sum(defect),sum(defect)/L
      close(iw)
      iw=fopen('one_particle_checks'//str(label_DMRG('u')),append=.false.)
      write(iw,*) '# spin trace_R min_eigenvalue max_eigenvalue max_diagonal_Im trace_eigenvalue_error'
      write(iw,*) '# Hermitian lower triangle reconstructed; cross-spin blocks zero by QN selection.'
      write(iw,*) '# Eigenvalues outside [0,1] by more than ',occupation_tol,' are rejected.'
      do s=1,Nspin
        write(iw,*) s,traces(s),lo(s),hi(s),diagonal_imag(s),trace_error(s)
      enddo
      close(iw)
    endif
    !
    deallocate(R,work,occupations)
  end subroutine get_one_particle_density






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
