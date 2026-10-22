! Post-processing of a saved BHZ ground state; no DMRG optimization is run.
! Build against the CMPLX nssDMRG library, with matching MPI flags.
program BHZ_1d_polarization
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
  integer                            :: Nso,Nsites,iorb,ispin,unit,comm,rank
  character(len=64)                  :: finput
  real(8)                            :: Eloc,Etot,Ekin,E0loc,Eint,Eshift
  real(8)                            :: mh,lambda
  type(site)                         :: Dot
  complex(8),dimension(4,4)          :: GammaZ,GammaX
  complex(8),dimension(4,4)          :: Hloc,Hlr
  integer                            :: i,j
  type(sparse_matrix),dimension(2,2) :: N
  type(sparse_matrix)                :: C,charge,spin
  type(sparse_matrix)                :: Sz,Tx,Ty,Tz,Tx2,Ty2,Tz2,Hint,H0loc,Hshift
  complex(8)                         :: z(2)
  real(8)                            :: phase
  logical                            :: master=.true.

#ifdef _MPI
  call init_MPI()
  comm=MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  master=get_Master_MPI(comm)
#endif
  call parse_cmd_variable(finput,'FINPUT',default='used.DMRG.conf')
  call parse_input_variable(mh,'MH',finput,default=0.5d0)
  call parse_input_variable(lambda,'LAMBDA',finput,default=0.3d0)
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
  do j=1,Nspin
     do i=1,Norb
        C=Dot%operators%op('C'//Dot%okey(i,j))
        N(i,j)=matmul(C%dgr(),C)
        call C%free()
     enddo
  enddo
  charge = N(1,1)+N(2,1)+N(1,2)+N(2,2)
  spin   =(N(1,1)+N(2,1)-N(1,2)-N(2,2))/2d0
  !
  call Init_Measure_DMRG()
  ! z_sigma=<exp[i*2*pi/L*sum_j j*n_j,sigma]>, j=1,...,L.
  ! Use the full chain and keep the raw origin/parity convention.
  z(1)=Measure_Polarization_DMRG(N(1,1)+N(2,1))
  z(2)=Measure_Polarization_DMRG(N(1,2)+N(2,2))
  !
  if(master)then
     unit=fopen('z_sigma'//str(label_DMRG('u')),append=.false.)
     write(unit,*) '# spin (1=up,2=down), Re z, Im z, abs(z), arg(z) [rad], P=arg(z)/(2*pi)'
     do j=1,2
        phase=atan2(aimag(z(j)),real(z(j),8))
        write(unit,*) j,dreal(z(j)),dimag(z(j)),abs(z(j)),phase,phase/(2d0*pi)
     enddo
     close(unit)
  endif
  !
  !
  call End_Measure_DMRG()
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

end program BHZ_1d_polarization
