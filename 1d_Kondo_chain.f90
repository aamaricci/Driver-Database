!       \                 /
!          \           / 
!             \  _  / 
!       1--2--3--4--5--6--7
!      _o _o _o  o _o _o _o
!      _I-_I-_I--I-_I-_I-_I
!_==optional
program Kondo1d
  USE EDLAT
  USE SCIFOR
#ifdef _MPI
  USE MPI
#endif
  implicit none
  character(len=16)   :: finput
  real(8)             :: ts,vtrap,atrap
  integer             :: N,N1,i,indi
  logical             :: pbc
  !
#ifdef _MPI
  call init_MPI
#endif
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="chain hopping parameter")
  call parse_input_variable(atrap,"atrap",finput,default=0d0,comment="optical trap bottom energy ")
  call parse_input_variable(vtrap,"vtrap",finput,default=0d0,comment="optical trap potential amplitude: 1/2.V.x**2")
  call parse_input_variable(pbc,"PBC",finput,default=.false.,comment="T: PBC, F: OBC")
  call ed_read_input(trim(finput))
  !
  if(any(Nsites(1:Norb)==0))stop "This driver is for Kondo problem only: Nsites=[Ns,..,iNs]"  
  if(size(Nsites)/=Norb+1)stop "This driver is for Kondo problem only: Nsites=[Ns,..,iNs]"
  if( Nsites(Norb+1) > Nsites(1) )stop "This driver is for Kondo problem only: iNs > eNs"
  call ed_Hij_init(Nsites(1:Norb))

  N  = Nsites(1)                !odd
  N1 = (N+1)/2                  !N1%2==0
  call ed_Hij_add_link(1,2,1,1,1,one*ts)
  do i=2,N-1
     call ed_Hij_add_link(i,i-1,1,1,1,one*ts)
     call ed_Hij_add_link(i,i+1,1,1,1,one*ts)
  enddo
  call ed_Hij_add_link(N,N-1,1,1,1,one*ts)
  !PBC
  if(pbc)then
     call ed_Hij_add_link(1,N,1,1,1,one*ts)
     call ed_Hij_add_link(N,1,1,1,1,one*ts)
  end if
  if(atrap/=0d0.OR.vtrap/=0d0)then
     do i=1,N
        call ed_Hij_add_link(i,i,1,1,1,one*e_local(i))
     enddo
  endif
  !
  call ed_Hij_info()
  !
  do i=1,Nsites(Norb+1)
     indi=N1-Nsites(Norb+1)/2 + (i-1)
     if(Nsites(Norb+1)==Nsites(1))indi=i
     Jkindx(i)=indi
  enddo

  call ed_init_solver()
  call ed_solve()


#ifdef _MPI
  call finalize_MPI()
#endif

contains

  function e_local(i) result(e)
    integer :: i
    real(8) :: e
    real(8) :: r2
    r2 = dble(i-N1)**2
    e  = atrap + 0.5d0*Vtrap*r2
  end function e_local

end program Kondo1d












