program testEDkron
  USE SCIFOR
  USE DMRG, id=>sp_eye
#ifdef _MPI
  USE MPI
#endif
  implicit none

  integer                                      :: Nso
  character(len=64)                            :: finput
  integer                                      :: i,j,Nsb,SUN,iter,Niter
  real(8),dimension(:),allocatable             :: Vec,Hvec,Hvec_
  type(site),dimension(:),allocatable :: My_Dot
  type(sparse_matrix)                          :: C,N,Cl,Cr,P
  type(sparse_matrix),allocatable,dimension(:) :: Hsb
  real(8),dimension(:,:),allocatable  :: Hlr
  real(8)                                      :: t0_
  real(8)                                        :: ts,Mh,lambda,val,K,alpha
  integer                                      :: m_sb
  real(8),dimension(:,:),allocatable           :: Evecs,Rho,Hloc
  real(8),dimension(:),allocatable             :: Evals
  integer                                      :: Neigen=2,vecDim
  integer                             :: irank,comm,rank,ierr
  logical                             :: master

#ifdef _MPI  
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
#endif


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(ts,"TS",finput,default=-0.5d0,comment="Hopping amplitude")
  call parse_input_variable(alpha,"alpha",finput,default=1d0,comment="bandwidth ratio")
  call parse_input_variable(Mh,"MH",finput,default=0d0,comment="Crystal field splittings")
  call parse_input_variable(lambda,"LAMBDA",finput,default=0d0,comment="off-diagonal amplitude")
  call parse_input_variable(Niter,"Niter",finput,default=5,&
       comment="number of blocks enlarge call == size of the system to solve: 2^Niter")
  call read_input(finput)



  Nso = Nspin*Norb
  allocate(Hloc(Nso,Nso))
  allocate(Hlr(Nso,Nso))
  select case(Norb)
  case(1)
     Hloc = Mh*pauli_z          !use it as a Zeeman field
     Hlr  = ts*pauli_0
  case(2)                       !spin x orbital <= ext x int
     Hloc = one*Mh*kron(pauli_0,pauli_z)
     Hlr  = one*ts*kron(pauli_0,diag([1d0,alpha])) + one*lambda/2d0*kron(pauli_0,pauli_x)
  case default;stop "This code is for Norb<=2. STOP"
  end select
  if(master)then
     call print_matrix(Hloc,"Hloc.dmrg")
     call print_matrix(Hlr,"Hlr.dmrg")
  endif

  allocate(My_Dot(1))
  My_Dot(1) = electron_site(Hloc)


  !Init DMRG
  call init_dmrg(Hlr,ModelDot=[My_Dot])


  target_qn = DMRG_qn
  !


  print*,""
  print*,""
  print*,"######################################"
  print*,"   o->o + o<-o"
  print*,"######################################"
  print*,""
  print*,""

  write(LOGfile,"(A22,2I12)")"Blocks Length (L-R) = ",left%length,right%length
  left  = block(My_Dot(1))
  right = block(My_Dot(1))
  do iter=1,Niter
     if(master)write(*,*)"Enlarge to:",2*iter
     call enlarge_block(left,My_Dot(1),label='left')
     call enlarge_block(right,My_Dot(1),label='right')
  enddo



  !#################################
  !    Build SUPER-BLOCK Sector
  !#################################
  current_L         = left%length + right%length
  current_target_QN = int(target_qn*current_L*Norb)
  if(master)then
     write(LOGfile,"(A22,I12)")"SuperBlock Length = ",current_L
     write(LOGfile,"(A22,"//str(size(current_target_QN))//"F12.7)")"Target_QN = ",current_target_QN
     write(LOGfile,"(A22,G12.7)")"Total         = ",sum(current_target_QN)
     write(LOGfile,"(A22,G12.7)")"Filling       = ",sum(current_target_QN)/current_L
     write(LOGfile,"(A22,G12.7)")"Filling/Norb  = ",sum(current_target_QN)/current_L/Norb
  endif
  call sb_get_states()
  if(master)print*,size(sb_states)



  sparse_H = .false.
  call sb_build_Hv()
  allocate(Evals(Neigen))
  ! allocate(Evecs(size(sb_states),Neigen))
  vecDim = sb_vecDim_Hv()
  allocate(Evecs(vecDim,Neigen))
  if(master)then
     t0_=t_start()
     call start_timer()
  endif
  call sp_eigh(MPI_COMM_WORLD,spMatVec_MPI_direct_main,evals,evecs,&
       3*Neigen,&
       500,&
       tol=1d-12,&
       iverbose=.false.)
  if(master)then
     print*,"sp_peigh Time:",t_stop()
     call stop_timer()
     do i=1,Neigen
        print*,i,Evals(i)/2/left%length/Norb
        write(11,*)i,Evals(i)/2/left%length/Norb
     enddo
  endif
  deallocate(evals,evecs)
  call sb_delete_Hv()
  print*,""



  !Finalize DMRG
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif





end program testEDkron
