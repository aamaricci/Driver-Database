! Post-processing driver for the interacting 1D BHZ model.
! Use with irun=.false. in DMRG.conf: the ground state is read from restart.
program BHZ_1d_post
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
  call Measure_DMRG([Tx2,Ty2,Tz2],file="orbital_fluctuations",pos=arange(1,Nsites))
  call Measure_DMRG([Tx,Ty],file="in_plane_orbital_polarization",pos=arange(1,Nsites))
  !
  if(master)print*,"<Sz_i.Sz_j>_c"
  call get_correlations("sz.sz",Sz)
  if(master)print*,"<Tz_i.Tz_j>_c"     
  call get_correlations("tz.tz",Tz)
  if(master)print*,"<Tx_i.Tx_j>_c"
  call get_correlations("tx.tx",Tx)
  if(master)print*,"<Ty_i.Ty_j>_c"
  call get_correlations("ty.ty",Ty)
  !
  !Get T_\perp(ij)=<T_x(i).T_x(j)>_c + <T_y(i).T_y(j)>_c
  if(master)print*,"<Tperp_i.Tperp_j>_c"
  call get_perp_correlations(Tx,Ty)
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

  subroutine get_correlations(label,Op)
    character(len=*),intent(in)    :: label
    type(sparse_matrix),intent(in) :: Op
    integer                        :: i,j,r,ic
    real(8)                        :: value
    if(master) unit=fopen(str(label)//"_r"//str(label_DMRG('u')),append=.false.)
    ic=Nsites/2
    if(master)call start_timer("<"//str(label)//">_|i-j|")
    do r=1,ic-1
      i=ic-r/2
      j=ic+(r+1)/2
      value=dreal(Measure_Corr_DMRG(Op,[0d0,0d0],Op,[0d0,0d0],i,j,connected=.true.))
      if(master) write(unit,*) abs(i-j),value
      if(master) call eta(r,ic-1)
    enddo
    if(master)call stop_timer()
    if(master) close(unit)
  end subroutine get_correlations

  
  ! C_perp(r)=C_xx(r)+C_yy(r), connected.  The final line stores the
  ! integrated susceptibility sum_{ij} C_perp(i,j)/Nsites.
  subroutine get_perp_correlations(OpX,OpY)
    type(sparse_matrix),intent(in) :: OpX,OpY
    integer                        :: i,j,r,ic
    real(8)                        :: value,chi
    if(master) unit=fopen("tperp.tperp_r"//str(label_DMRG('u')),append=.false.)
    if(master)call start_timer("<T_perp.Tperp>_|i-j|")
    ic=Nsites/2; 
    do r=1,ic-1
      i=ic-r/2; j=ic+(r+1)/2
      value=dreal(Measure_Corr_DMRG(OpX,[0d0,0d0],OpX,[0d0,0d0],i,j,connected=.true.))
      value=value+dreal(Measure_Corr_DMRG(OpY,[0d0,0d0],OpY,[0d0,0d0],i,j,connected=.true.))
      if(master) write(unit,*) abs(i-j),value
      if(master) call eta(r,ic-1)
    enddo
    if(master)call stop_timer()
    if(master) close(unit)
    !
  end subroutine get_perp_correlations



end program BHZ_1d_post
