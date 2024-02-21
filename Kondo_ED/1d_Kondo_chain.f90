!       1--2--3--4--5--6--7
!      _o _o _o  o _o _o _o
!      _I-_I-_I--I-_I-_I-_I
!_==optional
program Kondo1d
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

  !
#ifdef _MPI
  call init_MPI
#endif
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="chain hopping parameter")
  call parse_input_variable(timp,"TIMP",finput,default=0.d0,comment="impurity hopping parameter")
  call parse_input_variable(mue,"MUE",finput,default=0d0,comment="local energies at impurity sites")
  call parse_input_variable(mug,"MUG",finput,default=0d0,comment="local energies of the bath at impurity sites")
  call parse_input_variable(pbc,"PBC",finput,default=.false.,comment="T: PBC, F: OBC")
  call ed_read_input(trim(finput))
  !
  !
  if(any(Nsites(1:Norb)==0))stop "This driver is for Kondo problem only: Nsites=[Ns,..,iNs]"  
  if(size(Nsites)/=Norb+1)stop "This driver is for Kondo problem only: Nsites=[Ns,..,iNs]"
  if( Nsites(Norb+1) > Nsites(1) )stop "This driver is for Kondo problem only: iNs > eNs"
  if( Norb>1 )stop "This driver is for Kondo problem with: Norb == 1"
  !
  !> INIT H(Ri,Rj) matrix with Nsites[1:Norb,1==Imp]
  call ed_Hij_init(Nsites,Nspin)
  !
  !> BUILD BATH PART
  N   = Nsites(1)  !odd
  N1  = (N+1)/2    !N1%2==0
  do ispin=1,Nspin
     call ed_Hij_add_link(1,2,1,1,ispin,one*ts)
     do i=2,N-1
        call ed_Hij_add_link(i,i-1,1,1,ispin,one*ts)
        call ed_Hij_add_link(i,i+1,1,1,ispin,one*ts)
     enddo
     call ed_Hij_add_link(N,N-1,1,1,ispin,one*ts)
     if(pbc)then
        call ed_Hij_add_link(1,N,1,1,ispin,one*ts)
        call ed_Hij_add_link(N,1,1,1,ispin,one*ts)
     end if
     do i=1,Nsites(2)
        if(Mindx(i)==0)cycle;j=Mindx(i)
        call ed_Hij_add_link(j,j,1,1,ispin,one*mug)
     enddo
  enddo
  !
  !
  !> BUILD THE IMPURITY PART
  N   = Nsites(2)      !odd
  do ispin=1,Nspin
     if(N>1)then
        N1  = (N+1)/2 !N1%2==0
        call ed_Hij_add_link(1,2,2,2,ispin,one*timp)
        do i=2,N-1
           call ed_Hij_add_link(i,i-1,2,2,ispin,one*timp)
           call ed_Hij_add_link(i,i+1,2,2,ispin,one*timp)
        enddo
        call ed_Hij_add_link(N,N-1,2,2,ispin,one*timp)
        if(pbc)then
           call ed_Hij_add_link(1,N,2,2,ispin,one*timp)
           call ed_Hij_add_link(N,1,2,2,ispin,one*timp)
        end if
     endif
     do i=1,N
        if(Mindx(i)==0)cycle
        j=Mindx(i)
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

end program Kondo1d




!       \                 /
!          \           / 
!             \  _  / 

! call parse_input_variable(atrap,"atrap",finput,default=0d0,comment="optical trap bottom energy ")
! call parse_input_variable(vtrap,"vtrap",finput,default=0d0,comment="optical trap potential amplitude: 1/2.V.x**2")

! if(any([atrap,vtrap]/=0d0))then
!    do i=1,N
!       call ed_Hij_add_link(i,i,1,1,ispin,one*e_local(i,N1))
!    enddo
! endif

! if(any([atrap,vtrap]/=0d0))then
!    do i=1,N
!       call ed_Hij_add_link(i,i,2,2,ispin,one*e_local(i,N1))
!    enddo
! endif

! contains

!   function e_local(i,N) result(e)
!     integer :: i,N
!     real(8) :: e
!     real(8) :: r2
!     r2 = dble(i-N)**2
!     e  = atrap + 0.5d0*Vtrap*r2
!   end function e_local

