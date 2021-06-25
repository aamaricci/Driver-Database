program ss_DFT
  USE SLAVE_SPINS
  !
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI
  USE MPI
#endif
  implicit none

  integer                                 :: Nktot,Nkpath,Nkvec(3),Npts,Nlso
  integer                                 :: i,j,k,ik,ilat,iorb,jorb,io,ispin
  real(8),dimension(3)                    :: e1,e2,e3
  real(8),dimension(:,:),allocatable      :: kpath
  complex(8),dimension(:,:,:),allocatable :: Hk
  complex(8),dimension(:,:),allocatable   :: Hloc,Hk_f,Z_f,S_f
  real(8),allocatable                     :: Dens(:),Zeta(:), Self(:),Tmp(:)
  character(len=60)                       :: InputFile,latfile,ineqfile,hkfile,OrderFile
  character(len=40),allocatable           :: points_name(:)
  real(8)                                 :: ef
  logical                                 :: zHkflag
  logical                                 :: master=.true.,bool
  logical                                 :: bool_hk
  logical                                 :: bool_lat
  logical                                 :: bool_ineq
  logical                                 :: bool_order
  integer                                 :: unit
  integer,allocatable,dimension(:)        :: ineq_sites
  integer,dimension(3)                    :: Nin_w90
  character(len=5),dimension(3)           :: OrderIn_w90

#ifdef _MPI
  call init_MPI
  master = get_master_MPI()
#endif

  call parse_cmd_variable(InputFile,"INPUTFILE",default="inputDFT.conf")
  call parse_input_variable(hkfile,"hkfile",InputFile,default="hk.conf")
  call parse_input_variable(latfile,"latfile",InputFile,default="lat.conf")
  call parse_input_variable(ineqfile,"ineqfile",InputFile,default="ineq.conf")
  call parse_input_variable(orderfile,"Orderfile",InputFile,default="order.conf")
  call parse_input_variable(Nin_w90,"Nin_w90",InputFile,default=[Nspin,Norb,Nlat])
  call parse_input_variable(zHkflag,"zHkflag",InputFile,default=.true.)
  call parse_input_variable(Nkvec,"NKVEC",InputFile,default=[10,10,10])
  call parse_input_variable(nkpath,"NKPATH",InputFile,default=500)
  call ss_read_input(reg(InputFile))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  Nlso = Nlat*Nspin*Norb

  allocate(Zeta(Nlso))
  allocate(Self(Nlso))

  inquire(file=reg(hkfile),exist=bool_hk)
  inquire(file=reg(latfile),exist=bool_lat)
  inquire(file=reg(ineqfile),exist=bool_ineq)
  inquire(file=reg(orderfile),exist=bool_order)


  !Get/Set Wannier ordering:
  if(bool_order)then
     open(free_unit(unit),file=reg(orderfile))
     read(unit,*)OrderIn_w90(1),OrderIn_w90(2),OrderIn_w90(3)
     close(unit)
  else
     OrderIn_w90=[character(len=5)::"Nspin","Norb","Nlat"]
  endif

  !Setup inequivalent sites in the unit cell
  allocate(ineq_sites(Nlat));ineq_sites=1
  if(bool_ineq)then
     open(free_unit(unit),file=reg(ineqfile))
     do i=1,Nlat
        read(unit,*)ineq_sites(i)
        write(*,"(A,I5,A,I5)")"Site",i,"corresponds to ",ineq_sites(i)
     enddo
     close(unit)
  endif


  !Set basis vectors:
  if(bool_lat)then
     open(free_unit(unit),file=reg(latfile))
     read(unit,*)e1
     read(unit,*)e2
     read(unit,*)e3
     close(unit)
  else
     e1 = [1d0,0d0,0d0]
     e2 = [0d0,1d0,0d0]
     e3 = [0d0,0d0,1d0]
  endif
  call TB_set_ei(e1,e2,e3)
  call TB_build_bk(verbose=.true.)


  if(bool_hk)then
     call TB_read_hk(Hk,reg(hkfile),Nkvec)
     if(size(Hk,1)/=Nlat*Nspin*Norb)stop "ss_DFT error: wrong size in Hk as read from file"
  else
     write(*,*)"ERROR: file",trim(hkfile)," not found. Stop"
     stop
  endif

  write(*,*)"Using Nk_total="//str(size(Hk,3))
  allocate(Hloc(Nlso,Nlso))
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  if(master)call TB_write_Hloc(Hloc,"Hloc.dat")

  !########################################
  !SOLVE SS:
  call start_timer
  call ss_solve(Hk,ineq_sites=ineq_sites)
  call stop_timer("SS SOLUTION")
  !Retrieve Zeta and ReSigma(0)=lambda0-lambda
  call ss_get_zeta(zeta)
  call ss_get_Self(self)
  call save_array("renorm.save",[zeta,self])
  !########################################


  !Store the renormalized Hamiltonian:
  if(zHkflag)then
     call start_timer
     allocate(Hk_f(Nlso,Nlso),Z_f(Nlso,Nlso),S_f(Nlso,Nlso))
     Z_f = diag(sqrt(zeta))
     S_f = diag(self)
     do ik=1,Nktot
        Hk_f       = Hk(:,:,ik) - Hloc
        Hk_f       = (Z_f .x. Hk_f) .x. Z_f
        Hk(:,:,ik) = Hk_f + Hloc + S_f
     enddo
     deallocate(Hk_f,S_f,Z_f)
     call TB_write_hk("z"//reg(hkfile),Nkvec)
     call stop_timer("SS get zHk")
  endif


  deallocate(zeta,self)
  call TB_w90_delete()


#ifdef _MPI
  call finalize_MPI()
#endif

end program ss_DFT
