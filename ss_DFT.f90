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
  real(8),dimension(:,:),allocatable      :: kpath,kgrid
  complex(8),dimension(:,:,:),allocatable :: Hk
  complex(8),dimension(:,:),allocatable   :: Hloc
  real(8),allocatable                     :: Dens(:),Zeta(:), Self(:),Tmp(:)
  character(len=60)                       :: w90file,InputFile,latfile,kpathfile,ineqfile,hkfile,OrderFile
  character(len=40),allocatable           :: points_name(:)
  real(8)                                 :: ef
  logical                                 :: FSflag,Spinor,Bandsflag,zHkflag
  logical                                 :: master=.true.,bool
  logical                                 :: bool_hk
  logical                                 :: bool_lat
  logical                                 :: bool_kpath
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
  call parse_input_variable(w90file,"w90file",InputFile,default="hij.conf")
  call parse_input_variable(hkfile,"hkfile",InputFile,default="hk.conf")
  call parse_input_variable(latfile,"latfile",InputFile,default="lat.conf")
  call parse_input_variable(kpathfile,"kpathfile",InputFile,default="kpath.conf")
  call parse_input_variable(ineqfile,"ineqfile",InputFile,default="ineq.conf")
  call parse_input_variable(orderfile,"Orderfile",InputFile,default="order.conf")
  call parse_input_variable(Spinor,"Spinor",InputFile,default=.false.)
  call parse_input_variable(Bandsflag,"Bandsflag",InputFile,default=.true.)
  call parse_input_variable(zHkflag,"zHkflag",InputFile,default=.false.)
  call parse_input_variable(FSflag,"FSflag",InputFile,default=.false.)
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
  Nktot=product(Nkvec)

  allocate(Zeta(Nlso))
  allocate(Self(Nlso))

  inquire(file=reg(hkfile),exist=bool_hk)
  inquire(file=reg(latfile),exist=bool_lat)
  inquire(file=reg(kpathfile),exist=bool_kpath)
  inquire(file=reg(ineqfile),exist=bool_ineq)
  inquire(file=reg(orderfile),exist=bool_order)


  !Get/Set Wannier ordering:
  if(bool_order)then
     open(free_unit(unit),file=reg(orderfile))
     read(unit,*)Nin_w90(1),Nin_w90(2),Nin_w90(3)
     read(unit,*)OrderIn_w90(1),OrderIn_w90(2),OrderIn_w90(3)
     close(unit)
  else
     write(*,"(A)")"Using default order for W90 input: [Nspin, Norb, Nlat]"
     Nin_w90    =[Nspin,Norb,Nlat]
     OrderIn_w90=[character(len=5)::"Nspin","Norb","Nlat"]
  endif

  !Setup the path in the BZ.
  if(bool_kpath)then
     Npts = file_length(reg(kpathfile))
     allocate(kpath(Npts,3))
     allocate(points_name(Npts))
     open(free_unit(unit),file=reg(kpathfile))
     do i=1,Npts
        read(unit,*)points_name(i),kpath(i,:)
     enddo
     close(unit)
  else
     write(*,"(A)")"Using default path for 3d BZ: [M,R,Gm,X,M,Gm,Z,A,R]"
     Npts = 9
     allocate(kpath(Npts,3),points_name(Npts))
     kpath(1,:)=[0.5d0,0.5d0,0d0]
     kpath(2,:)=[0.5d0,0.5d0,0.5d0]
     kpath(3,:)=[0d0,0d0,0d0]
     kpath(4,:)=[0.5d0,0d0,0d0]
     kpath(5,:)=[0.5d0,0.5d0,0d0]
     kpath(6,:)=[0d0,0d0,0d0]
     kpath(7,:)=[0d0,0d0,0.5d0]
     kpath(8,:)=[0.5d0,0d0,0.5d0]
     kpath(9,:)=[0.5d0,0.5d0,0.5d0]
     points_name=[character(len=40) ::'M', 'R', '{/Symbol} G', 'X', 'M', '{/Symbol} G', 'Z','A', 'R']
  endif


  !Setup inequivalent sites in the unit cell
  allocate(ineq_sites(Nlat))
  if(bool_ineq)then
     open(free_unit(unit),file=reg(ineqfile))
     do i=1,Nlat
        read(unit,*)ineq_sites(i)
        write(*,"(A,I5,A,I5)")"Site",i,"corresponds to ",ineq_sites(i)
     enddo
     close(unit)
  else
     write(*,"(A)")"Using default Ineq_sites list: all equivalent to 1"
     ineq_sites = 1
  endif


  !Set basis vectors:
  if(bool_lat)then
     open(free_unit(unit),file=reg(latfile))
     read(unit,*)e1
     read(unit,*)e2
     read(unit,*)e3
     close(unit)
  else
     write(*,"(A)")"Using default lattice basis: ortho-normal"
     e1 = [1d0,0d0,0d0]
     e2 = [0d0,1d0,0d0]
     e3 = [0d0,0d0,1d0]
  endif
  call TB_set_ei(e1,e2,e3)
  call TB_build_bk(verbose=.true.)


  !Setup Wannier90 or read H(k) from file:
  call start_timer
  call TB_w90_setup(reg(w90file),nlat=Nlat,norb=Norb,nspin=Nspin,Spinor=spinor,verbose=.true.)
  call stop_timer("TB_w90_setup")
  if(bool_hk)then
     call TB_read_hk(Hk,reg(hkfile),Nlat,Nspin,Norb,Nkvec,kgrid)
     call assert_shape(Hk,[Nlso,Nlso,product(Nkvec)])
  else
     call start_timer  
     call TB_w90_FermiLevel(Nkvec,filling,Ef)
     call stop_timer("TB_w90_FermiLevel")
     !
     allocate(Hk(Nlso,Nlso,Nktot))
     call start_timer
     !Build H(k) and re-order it to the default DMFT_tools order:
     call TB_build_model(Hk,Nlso,Nkvec)
     Hk = TB_reshape_array(Hk,Nin=Nin_w90,OrderIn=OrderIn_w90,&
          OrderOut=[character(len=5)::"Norb","Nspin","Nlat"])
     call TB_write_hk(Hk,reg(hkfile),Nlat,Nspin,Norb,Nkvec)
     call stop_timer("TB_build_model")
  endif


  write(*,*)"Using Nk_total="//str(size(Hk,3))
  allocate(Hloc(Nlso,Nlso))
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  if(master)call TB_write_Hloc(Hloc,"w90Hloc.dat")




  !########################################
  !SOLVE SS: explicitly set User Order, so SS re-organize H(k) according to its needs.
  !note: we could have either i) re-order H(k) directly in SS order (it requires the user
  !to know it) or ii) reorder H(k) externally using again TB_reshape_array.
  !SS order is: Norb,Nlat,Nspin
  call start_timer
  call ss_solve(Hk,ineq_sites=ineq_sites,UserOrder=[character(len=5)::"Norb","Nspin","Nlat"])
  call stop_timer("SS SOLUTION")
  call ss_get_zeta(zeta)
  call ss_get_Self(self)
  call TB_w90_Zeta(zeta)
  call TB_w90_Self(diag(self))
  if(master)call save_array("renorm.save",[zeta,self])
  !########################################




  !Solve for the renormalized bands:
  if(BandsFlag)then
     call start_timer
     if(master)call TB_Solve_model(TB_w90_model,Nlso,kpath,Nkpath,&
          colors_name=[black,red,green,blue,magenta,black,red,green,blue,magenta],&
          points_name=points_name,& 
          file="zBands_ssDFT",iproject=.true.)
     call stop_timer("SS get zBands")
  endif

  !Store the renormalized Hamiltonian:
  if(zHkflag)then
     call start_timer
     call TB_build_model(Hk,Nlso,Nkvec)
     Hk = TB_reshape_array(Hk,Nin=Nin_w90,OrderIn=OrderIn_w90,&
          OrderOut=[character(len=5)::"Norb","Nspin","Nlat"])
     call TB_write_hk("z"//reg(hkfile),Nkvec)
     call stop_timer("SS get zHk")
  endif



  if(FSflag)then
     inquire(file='renorm.save',exist=bool)
     if(bool)then
        allocate(tmp(2*Nlat*Nspin*Norb))
        call read_array("zeta_self.restart",tmp)
        zeta = tmp(:Nlat*Nspin*Norb)
        self = tmp(Nlat*Nspin*Norb+1:)
        call TB_w90_Zeta(zeta)
        call TB_w90_Self(diag(self))
     endif
     call TB_FSurface(Nlso,0d0,Nkvec(1:2),&
          colors_name=[black,red,red,green,blue],&
          file='FS_ssDFT',cutoff=1d-1,Niter=3,Nsize=2)
  endif


  deallocate(zeta,self)
  call TB_w90_delete()


#ifdef _MPI
  call finalize_MPI()
#endif

end program ss_DFT
