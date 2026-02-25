program ssh1d
  USE SPIN_ED
  USE SCIFOR
#ifdef _MPI
  USE MPI
#endif
  implicit none
  character(len=16)   :: finput
  real(8)             :: bz,j,j0,dj,spin_gap,E0
  real(8),allocatable :: egs(:)
  integer             :: N,N1,i,Len
  logical             :: pbc,master=.true.
  !
#ifdef _MPI
  call init_MPI
  master = get_Master_MPI()
#endif
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputSSH.conf')
  call parse_input_variable(Bz,"Bz",finput,default=0d0,comment="local effective field")
  call parse_input_variable(J0,"J0",finput,default=1d0,comment="spin hopping base")
  call parse_input_variable(dJ,"dJ",finput,default=0d0,comment="spin hopping fluctation")
  call parse_input_variable(pbc,"PBC",finput,default=.false.,comment="T: PBC, F: OBC")
  call ed_read_input(trim(finput))
  !
  if( Norb>1 )stop "This driver is for SSH problem with: Norb == 1"
  !
  !> INIT H(Ri,Rj) matrix with Nsites[1:Norb,1==Imp]
  call ed_Hij_init(Nsites)
  !> BUILD H(Ri,Rj)_aa
  N   = Nsites(1)                !odd
  N1  = (N+1)/2                  !N1%2==0
  do i=2,N
     call ed_Hij_add_link(i-1,i,1,1,jspin(i-1))
     call ed_Hij_add_link(i,i-1,1,1,jspin(i-1))
  enddo
  if(pbc)then
     call ed_Hij_add_link(N,1,1,1,jspin(N))
     call ed_Hij_add_link(1,N,1,1,jspin(N))
  end if
  !
  !PRINT INFO H(Ri,Rj)
  call ed_Hij_info()
  !
  !> SOLVE THE SPIN PROBLEM
  call ed_init_solver()
  call ed_solve()
  !
  if(master)then
     len=file_length("list_gs_energy.ed")
     allocate(egs(len))
     open(100,file="list_gs_energy.ed")
     do i=1,len
        read(100,*)egs(i)
     enddo
     close(100)
     E0 = minval(egs)
     open(100,file="EgsVSL.dat")
     write(100,*)N,E0/N
     close(100)
     !
     egs = egs-E0
     spin_gap = minval(egs-E0,mask=egs>=0d0)
     open(100,file="spin_gap.dat")
     write(*,*)"Spin Gap=",spin_gap
     write(100,*)spin_gap
     close(100)
  endif
  !
  !
#ifdef _MPI
  call finalize_MPI()
#endif

contains

  function jspin(i) result(j)
    integer    :: i
    complex(8) :: j
    J=one*(J0-dJ)
    if(mod(i,2)==0)J=one*(J0+dJ)
  end function jspin

end program Ssh1d












