! Post-processing of a saved BHZ ground state; no DMRG optimization is run.
! Build against the CMPLX nssDMRG library, with matching MPI flags.
program BHZ_1d_structure
  USE SF_SPIN
  USE SF_CONSTANTS
  USE SF_ARRAYS
  USE SF_LINALG
  USE SF_IOTOOLS
  USE DMRG
  use DMRG_GLOBAL, only: current_L
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
  complex(8),dimension(4,4)                      :: Hloc,Hlr
  type(sparse_matrix)                            :: C,charge,spin
  integer                                        :: i,j
  integer                                        :: cut,k
  integer,allocatable                            :: positions(:)
  real(8),allocatable                            :: q(:)
  complex(8),allocatable                         :: F(:)
  logical                                        :: connected
  type(sparse_matrix),dimension(:,:),allocatable :: N
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
  call parse_input_variable(cut,'BULK_CUT',finput,default=0,comment="Index to cut out extrema of the 1d lattice 1+bulk_cut:L-bulk_cut")
  call parse_input_variable(connected,'STRUCTURE_CONNECTED',finput,default=.false.,comment="Whether to include connected correlations")
  call read_input(finput)
  !
  if(Nspin/=2 .or. Norb/=2) stop "BHZ post driver requires Nspin=Norb=2"
  Nsites=2*Ldmrg
  Nso=Nspin*Norb
  !
  GammaX=kron(pauli_sigma_z,pauli_tau_x)
  GammaZ=kron(pauli_sigma_0,pauli_tau_z)
  !
  Hloc=mh*GammaZ
  !
  Dot=electron_site(Hloc,H0loc=H0loc,Hint=Hint,Hshift=Hshift)
  !
  Hlr=-0.5d0*GammaZ+0.5d0*xi*lambda*GammaX
  !
  call init_dmrg(Hlr,ModelDot=[Dot])
  !
  allocate(N(Norb,Nspin))
  do j=1,Nspin
     do i=1,Norb
        C=Dot%operators%op('C'//Dot%okey(i,j))
        N(i,j)=matmul(C%dgr(),C)
        call C%free()
     enddo
  enddo
  !
  charge = N(1,1)+N(2,1)+N(1,2)+N(2,2)
  spin   =(N(1,1)+N(2,1)-N(1,2)-N(2,2))/2d0
  

  call Init_Measure_DMRG()
  if(cut<0.OR.2*cut>=current_L)error stop 'Invalid BULK_CUT'

  !
  positions=[(i,i=cut+1,current_L-cut)]
  q        =[(2d0*pi*dble(k)/dble(size(positions)),k=0,size(positions)-1)]
  !
  ! F(q)=sum_ij exp[i*q*(j-i)] <A_i A_j> / number of selected sites.
  ! Equal-site terms and all pair multiplicities are retained, also in OBC.
  if(Master)write(*,*)"Computing spin structure factor S(q)" 
  call Measure_Structure_DMRG(spin,[0d0,0d0],spin,[0d0,0d0],q,F, &
       positions=positions,&
       connected=connected,&
       same_operator=.true.,hermitian=.true.,commuting=.true.)
  unit=fopen('Sq'//str(label_DMRG('u')),append=.false.)
  write(unit,*) '# q [rad/site], Re F(q), Im F(q)', 'connected=',connected,' L=',current_L,' BULK_CUT=',cut
  do k=1,size(q)
     write(unit,*) q(k),dreal(F(k)),dimag(F(k))
  enddo
  close(unit)
  !
  if(Master)write(*,*)"Computing charge structure factor N(q)" 
  call Measure_Structure_DMRG(charge,[0d0,0d0],charge,[0d0,0d0],q,F, &
       positions=positions,&
       connected=connected,&
       same_operator=.true.,hermitian=.true.,commuting=.true.)
  unit=fopen('Nq'//str(label_DMRG('u')),append=.false.)
  write(unit,*) '# q [rad/site], Re F(q), Im F(q)', 'connected=',connected,' L=',current_L,' BULK_CUT=',cut
  do k=1,size(q)
    write(unit,*) q(k),dreal(F(k)),dimag(F(k))
  enddo
  close(unit)     
  !
  !
  call End_Measure_DMRG()
  !
  !Free memory
  do j=1,Nspin
     do i=1,Norb
        call N(i,j)%free()
     enddo
  enddo
  call charge%free()
  call spin%free()
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif

end program BHZ_1d_structure
