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
  call parse_input_variable(SUN,"SUN",finput,default=2,&
       comment="Spin SU(N) value. 2=> spin 1/2, 3=> spin 1")
  call parse_input_variable(DMRGtype,"DMRGtype",finput,default="infinite",&
       comment="DMRG algorithm: Infinite, Finite")
  call parse_input_variable(Niter,"Niter",finput,default=5,&
       comment="number of blocks enlarge call == size of the system to solve: 2^Niter")
  call read_input(finput)



  allocate(My_Dot(1))
  My_Dot(1) = spin_site(sun=SUN,Hvec=Hvec)


  !Init DMRG
  Hlr = diag([Jp,Jx/2d0])
  call init_dmrg(Hlr,ModelDot=My_Dot)


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
