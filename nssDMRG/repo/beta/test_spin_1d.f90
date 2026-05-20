program test_DMRG
  USE SCIFOR
  USE DMRG, id=>sp_eye
#ifdef _MPI
  USE MPI
#endif
  implicit none

  integer                                      :: Nso
  character(len=64)                            :: finput
  integer                                      :: i,j,Nsb,SUN,iter
  real(8),dimension(:),allocatable             :: Vec,Hvec,Hvec_
  integer                                      :: current_L
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
  logical                             :: master,bool_left,bool_right

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
  call read_input(finput)



  allocate(My_Dot(1))
  My_Dot(1) = spin_site(sun=SUN,Hvec=Hvec)


  !Init DMRG
  Hlr = diag([Jp,Jx/2d0])
  call init_dmrg(Hlr,ModelDot=My_Dot)


  target_qn = DMRG_qn
  !

  inquire(file=str(block_file//suffix_dmrg('left')//".restart"), exist=bool_left)
  inquire(file=str(block_file//suffix_dmrg('right')//".restart"), exist=bool_right)
  if(bool_left.AND.bool_right)then
     call left%load(str(block_file//suffix_dmrg('left')//".restart"))
     call right%load(str(block_file//suffix_dmrg('right')//".restart"))
  else
     left  = block(My_Dot(1))
     right = block(My_Dot(1))
  endif


  write(LOGfile,"(A22,2I12)")"Blocks Length (L-R) = ",left%length,right%length
  ! do iter=1,Ldmrg
  call enlarge_block(left,My_Dot(1),grow='left')
  call enlarge_block(right,My_Dot(1),grow='right')
  ! enddo

  !#################################
  !    Build SUPER-BLOCK Sector
  !#################################
  current_L         = left%length + right%length
  current_target_QN = int(target_qn*current_L*Norb)
  write(LOGfile,"(A22,I12)")"SuperBlock Length = ",current_L
  write(LOGfile,"(A22,"//str(size(current_target_QN))//"F12.7)")"Target_QN = ",current_target_QN
  write(LOGfile,"(A22,G12.7)")"Total         = ",sum(current_target_QN)
  write(LOGfile,"(A22,G12.7)")"Filling       = ",sum(current_target_QN)/current_L
  write(LOGfile,"(A22,G12.7)")"Filling/Norb  = ",sum(current_target_QN)/current_L/Norb
  call sb_get_states()
  print*,size(sb_states)

  call sb_diag()

  if(MpiMaster)then
     write(LOGfile,*)"- - - - - - - - - - - - - - - - - - - - -"
     select case(left%type())
     case ("fermion","f")
        write(LOGfile,"(A,"//str(Lanc_Neigen)//"F24.15)")&
             "Energies/N                           :",gs_energy/sum(current_target_QN)
     case ("spin","s")
        write(LOGfile,"(A,"//str(Lanc_Neigen)//"F24.15)")&
             "Energies/L                           :",gs_energy/current_L
     end select
     write(LOGfile,*)"- - - - - - - - - - - - - - - - - - - - -"
  endif


  !Finalize DMRG
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif














contains



  !This is user defined Function to be passed to the SYSTEM
  function spin_1d_hmodel(left,right) result(Hlr)
    type(block)                        :: left,right
    real(8),dimension(:,:),allocatable :: Hlr
    if(allocated(Hlr))deallocate(Hlr)
    allocate(Hlr(Nspin*Norb,Nspin*Norb))
    !
    !workout local part, like random local field
    !if(left%Dim==1 AND right%Dim==1) then operate over local H
    !if(left%Dim==1 OR right%Dim==1) then operate over local H
    Hlr(1,1) = Jp
    Hlr(2,2) = Jx/2d0
  end function spin_1d_hmodel


end program test_DMRG
