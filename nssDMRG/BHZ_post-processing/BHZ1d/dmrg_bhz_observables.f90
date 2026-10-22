! Build with the CMPLX library. Run in a fresh directory with FINPUT=DMRG.conf.
program bhz_observables
  use SCIFOR
  use DMRG
  use DMRG_GLOBAL, only: current_L
#ifdef _MPI
  use MPI
#endif
  implicit none
  type(site) :: model(1)
  type(sparse_matrix) :: c,n(2,2),density,spin
  complex(8) :: Hloc(4,4),Hlr(4,4),zup,zdown,cs,cp
  complex(8),allocatable :: F(:)
  real(8),allocatable :: q(:)
  real(8) :: mass,lambda
  integer :: i,j,L,cut,unit,comm,r,norigins,rmax
  integer,allocatable :: positions(:)
  character(len=256) :: finput
  logical :: master=.true.,get_structure,get_z,connected,irun,get_string
#ifdef _MPI
  call init_MPI()
  comm=MPI_COMM_WORLD
  master=get_Master_MPI(comm)
#endif
  call parse_cmd_variable(finput,'FINPUT',default='DMRG.conf')
  call parse_input_variable(mass,'MH',finput,default=0.5d0)
  call parse_input_variable(lambda,'LAMBDA',finput,default=0.3d0)
  call parse_input_variable(get_structure,'GET_STRUCTURE',finput,default=.true.)
  call parse_input_variable(get_z,'GET_Z',finput,default=.true.)
  call parse_input_variable(get_string,'GET_STRING',finput,default=.false.)
  call parse_input_variable(connected,'STRUCTURE_CONNECTED',finput,default=.false.)
  call parse_input_variable(irun,'IRUN',finput,default=.true.)
  call parse_input_variable(cut,'BULK_CUT',finput,default=0)
  call read_input(finput)
  if(Nspin/=2.or.Norb/=2)error stop 'BHZ requires Nspin=2, Norb=2'
  if(QNdim/=2.or.to_lower(DMRG_MODE)/='normal')error stop 'This example requires normal mode and QNDIM=2'
  save_umat=.true.;block_umat_cache=.true.
  Hloc=mass*kron(pauli_sigma_0,pauli_tau_z)
  Hlr=-0.5d0*kron(pauli_sigma_0,pauli_tau_z)&
       +0.5d0*xi*lambda*kron(pauli_sigma_z,pauli_tau_x)
  model(1)=electron_site(Hloc)
  do j=1,2
     do i=1,2
        c=model(1)%operators%op('C'//model(1)%okey(i,j))
        n(i,j)=matmul(c%dgr(),c)
     enddo
  enddo
  density=n(1,1)+n(2,1)+n(1,2)+n(2,2)
  spin=(n(1,1)+n(2,1)-n(1,2)-n(2,2))/2d0
  call init_dmrg(Hlr,ModelDot=model)
  if(irun)call run_DMRG()
  call Init_Measure_DMRG()
  if(get_structure)then
     L=current_L
     if(cut<0.or.2*cut>=L)error stop 'Invalid BULK_CUT'
     positions=[(i,i=cut+1,L-cut)]
     q=2d0*pi*arange(0,size(positions)-1)/dble(size(positions))
     call Measure_Structure_DMRG(density,[0d0,0d0],density,[0d0,0d0],q,F,&
          positions=positions,connected=connected)
     call write_factor('Nq.dat')
     call Measure_Structure_DMRG(spin,[0d0,0d0],spin,[0d0,0d0],q,F,&
          positions=positions,connected=connected)
     call write_factor('Sq.dat')
  endif
  if(get_z)then
     zup=Measure_Polarization_DMRG(n(1,1)+n(2,1))
     zdown=Measure_Polarization_DMRG(n(1,2)+n(2,2))
     if(master)then
        open(newunit=unit,file='z_sigma.dat',status='replace')
        write(unit,*) '# Eq.28: spin (1=up,2=down), Re z, Im z, |z|, arg z [rad]'
        write(unit,*)1,real(zup,8),aimag(zup),abs(zup),atan2(aimag(zup),real(zup,8))
        write(unit,*)2,real(zdown,8),aimag(zdown),abs(zdown),atan2(aimag(zdown),real(zdown,8))
        close(unit)
     endif
  endif
  if(get_string)then
     L=current_L
     if(PBCdmrg)error stop 'GET_STRING averages non-wrapping OBC intervals only'
     if(cut<0.or.2*cut>=L-1)error stop 'Invalid BULK_CUT for string correlations'
     if(master)then
        open(newunit=unit,file='string_order.dat',status='replace')
        write(unit,*) '# r, Re Cstring, Im Cstring, Re Cparity, Im Cparity, origins'
     endif
     call parse_input_variable(rmax,'STRING_RMAX',finput,default=(L-2*cut)/2)
     if(rmax<1.or.rmax>L-2*cut-1)error stop 'Invalid STRING_RMAX'
     do r=1,rmax
        cs=0d0;cp=0d0;norigins=L-2*cut-r
        do i=cut+1,L-cut-r
           cs=cs+Measure_String_DMRG(spin,i,i+r)
           cp=cp+Measure_Parity_DMRG(density,i,i+r)
        enddo
        cs=cs/norigins;cp=cp/norigins
        if(master)write(unit,*)r,real(cs,8),aimag(cs),real(cp,8),aimag(cp),norigins
     enddo
     if(master)close(unit)
  endif
  call End_Measure_DMRG()
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif
contains
  subroutine write_factor(filename)
    character(len=*),intent(in) :: filename
    integer :: k
    if(.not.master)return
    open(newunit=unit,file=filename,status='replace')
    write(unit,*) '# q [rad/site], Re F(q), Im F(q); normalized by window size'
    write(unit,*) '# connected = ',connected
    do k=1,size(q)
       write(unit,*)q(k),real(F(k),8),aimag(F(k))
    enddo
    close(unit)
  end subroutine
end program
