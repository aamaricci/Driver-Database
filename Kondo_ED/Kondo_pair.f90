!       o 
!       I
program KondoPair
  USE KONDO_ED
  USE SCIFOR
#ifdef _MPI
  USE MPI
#endif
  implicit none
  character(len=16)                :: finput
  real(8)                          :: ts,timp,mue,mug,Be,sg
  integer                          :: N,N1,i,indi,j,ispin
  logical                          :: pbc

#ifdef _MPI
  call init_MPI
#endif
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(mug,"MUG",finput,default=0d0,comment="local energies of the bath at impurity sites")
  call parse_input_variable(mue,"MUE",finput,default=0d0,comment="local energies at impurity sites")
  call ed_read_input(trim(finput))
  ! allocate(Mindx(Nsites(Norb+1)))
  ! allocate(Bindx(Nsites(Norb+1))) 
  ! call parse_input_variable(Mindx,"Mindx",finput,&
  !      default=(/( 0,i=1,size(Mindx) )/),&
  !      comment="labels of the sites where to apply potential shifts")
  ! call parse_input_variable(Bindx,"Bindx",finput,&
  !      default=(/( 0d0,i=1,size(Bindx) )/),&
  !      comment="local Magnetic field along Z at impurity sites")

  if(any(Nsites(1:Norb)==0))stop "This driver is for Kondo problem only: Nsites=[Ns,..,iNs]"  
  if(size(Nsites)/=Norb+1)stop "This driver is for Kondo problem only: Nsites=[Ns,..,iNs]"
  if( Nsites(Norb+1) > Nsites(1) )stop "This driver is for Kondo problem only: iNs > eNs"
  if( Norb>1 )stop "This driver is for Kondo problem with: Norb == 1"
  !
  !> INIT H(Ri,Rj) matrix with Nsites[1:Norb,1==Imp]
  call ed_Hij_init(Nsites,Nspin)

  !> BUILD BATH PART: H(Ri,Rj)_aa
  do ispin=1,Nspin
     do i=1,Nsites(2)
        if(Mindx(i)==0)cycle
        j=Mindx(i)
        call ed_Hij_add_link(j,j,1,1,ispin,one*mug)
        call ed_Hij_add_link(j,j,2,2,ispin,one*mue)
     enddo
     do i=1,N
        if(Bindx(i)==0d0)cycle
        sg = (-1d0)**ispin
        Be = Bindx(i)*sg
        call ed_Hij_add_link(i,i,2,2,ispin,one*Be)
     enddo
  enddo
  !
  !
  !PRINT INFO H(Ri,Rj)
  call ed_Hij_info()
  !
  !> SOLVE THE KONDO PROBLEM
  call ed_init_solver()
  call ed_solve()
  !
  !
  !
#ifdef _MPI
  call finalize_MPI()
#endif



end program KondoPair












