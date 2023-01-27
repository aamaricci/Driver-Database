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
  type(rgb_color),allocatable,dimension(:):: colors

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
  if(.not.bool_order)stop "order.conf file does not exist. STOP" 


  !Solve for the renormalized bands:
  if(BandsFlag)then
     allocate(colors(Nlso))
     select case(Nlso)
     case default
        colors = black
     case(2)
        colors = [blue,red]
     case(3)
        colors = [blue,red,green]
     case(4)
        colors = [blue,red,green,black]
     case(5)
        colors = [blue,red,green,black,magenta]
     end select
  end if


  !Get/Set Wannier ordering:
  if(bool_order)then
     open(free_unit(unit),file=reg(orderfile))
     read(unit,*)Nin_w90(1),Nin_w90(2),Nin_w90(3)
     read(unit,*)OrderIn_w90(1),OrderIn_w90(2),OrderIn_w90(3)
     close(unit)
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
     call sleep(1)
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
     Hk = TB_reorder_hk(Hk,Nin=Nin_w90,OrderIn=OrderIn_w90,&
          OrderOut=[character(len=5)::"Norb","Nlat","Nspin"]) !default SS ordering
     call TB_write_hk(Hk,reg(hkfile),Nlat,Nspin,Norb,Nkvec)
     call stop_timer("TB_build_model")
  endif


  write(*,*)"Using Nk_total="//str(size(Hk,3))
  allocate(Hloc(Nlso,Nlso))
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  if(master)call TB_write_Hloc(Hloc,"w90Hloc.dat")


  !Solve for the renormalized bands:
  if(BandsFlag)then
     call start_timer
     if(master)call TB_Solve_model(TB_w90_model,Nlso,kpath,Nkpath,&
          colors_name=[black,red,green,blue,magenta,black,red,green,blue,magenta],&
          points_name=points_name,& 
          file="Bands_DFT",iproject=.true.)
     call stop_timer("SS get zBands")
  endif


  !########################################
  !SOLVE SS
  !SS order is: Norb,Nlat,Nspin
  call start_timer
  call ss_solve(Hk,ineq_sites=ineq_sites)
  call stop_timer("SS SOLUTION")
  !########################################


  !Solve for the renormalized bands:
  if(BandsFlag)then
     call start_timer
     if(master)call TB_Solve_model(ss_Hk_model,Nlso,kpath,Nkpath,&
          colors_name=[black,red,green,blue,magenta,black,red,green,blue,magenta],&
          points_name=points_name,& 
          file="zBands_ssDFT",iproject=.false.)
     call stop_timer("SS get zBands")
  endif



  call TB_w90_delete()


#ifdef _MPI
  call finalize_MPI()
#endif




contains


  function ss_Hk_model(kvec,N) result(Hk)
    real(8),dimension(:)      :: kvec
    integer                   :: N
    complex(8),dimension(N,N) :: Hloc,Hk
    !< Get Hloc from W90:
    call TB_w90_Hloc(Hloc)
    !< Build H(k) from W90:
    Hk = TB_w90_model(kvec,N)
    !< Reorder H(k) according to W90-SS orders
    Hk  = TB_reorder_array(Hk,Nin=Nin_w90,OrderIn=OrderIn_w90,&
         OrderOut=[character(len=5)::"Norb","Nlat","Nspin"])
    Hloc= TB_reorder_array(Hloc,Nin=Nin_w90,OrderIn=OrderIn_w90,&
         OrderOut=[character(len=5)::"Norb","Nlat","Nspin"])
    !< Build effective fermionic H*(k) 
    call ss_get_ssHk(Hk,Hloc)
  end function ss_Hk_model



end program ss_DFT















! if(FSflag)then
!    inquire(file='renorm.save',exist=bool)
!    if(bool)then
!       allocate(tmp(2*Nlat*Nspin*Norb))
!       call read_array("zeta_self.restart",tmp)
!       zeta = tmp(:Nlat*Nspin*Norb)
!       self = tmp(Nlat*Nspin*Norb+1:)
!       call TB_w90_Zeta(zeta)
!       call TB_w90_Self(diag(self))
!    endif
!    call TB_FSurface(Nlso,0d0,Nkvec(1:2),&
!         colors_name=[black,red,red,green,blue],&
!         file='FS_ssDFT',cutoff=1d-1,Niter=3,Nsize=2)
! endif



! contains


!   subroutine write_hk_array(Hk,file,Nlat,Nspin,Norb,Nkvec)
!     character(len=*)                                                     :: file
!     integer                                                              :: Nlat,Nspin,Norb,Nlso
!     integer                                                              :: Nkvec(:)
!     real(8),dimension(product(Nkvec),size(Nkvec))                        :: kgrid ![Nk][Ndim]
!     real(8),dimension(3)                                                 :: kvec
!     integer                                                              :: Nktot,unit
!     integer                                                              :: i,ik,io,jo
!     complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,product(Nkvec)) :: Hk
!     logical :: mpi_master
!     !
!     mpi_master=.true.
!     !
!     call TB_build_kgrid(Nkvec,kgrid)
!     !
!     Nktot  = product(Nkvec)!==size(Hk,3)
!     Nlso   = Nlat*Nspin*Norb
!     !
!     if(mpi_master)then
!        open(free_unit(unit),file=reg(file))
!        write(unit,'(10(A10,1x))')str(Nktot),str(Nlat),str(Nspin),str(Norb)
!        write(unit,'(1A1,3(A12,1x))')"#",(reg(txtfy(Nkvec(ik))),ik=1,size(Nkvec))
!        do ik=1,Nktot
!           kvec=0d0 ; kvec(:size(Nkvec)) = kgrid(ik,:size(Nkvec))
!           write(unit,"(3(F15.9,1x))")(kvec(i),i=1,3) 
!           do io=1,Nlso
!              write(unit,"(1000(F5.2,1x))")(dreal(Hk(io,jo,ik)),jo=1,Nlso)
!           enddo
!        enddo
!        close(unit)
!     endif
!   end subroutine write_hk_array




!   function reorder_hk(Huser,Nin,OrderIn,OrderOut) result(Hss)
!     complex(8),dimension(:,:,:)                                     :: Huser
!     integer,dimension(3)                                            :: Nin   !In sequence of Nlat,Nspin,Norb as integers
!     character(len=*),dimension(3)                                   :: OrderIn  !in  sequence of Nlat,Nspin,Norb as strings
!     !
!     complex(8),dimension(size(Huser,1),size(Huser,2),size(Huser,3)) :: Hss
!     integer,dimension(3)                                            :: Ivec,Jvec
!     integer,dimension(3)                                            :: IndexOut
!     integer,dimension(3)                                            :: Nout
!     integer                                                         :: iss,jss,iuser,juser,i,Nlso,Nk
!     character(len=*),dimension(3),optional                          :: OrderOut !out sequence of Nlat,Nspin,Norb as strings
!     character(len=5),dimension(3)                                   :: OrderOut_ !out sequence of Nlat,Nspin,Norb as strings
!     !
!     OrderOut_=[character(len=5)::"Norb","Nspin","Nlat"];
!     if(present(OrderOut))then
!        do i=1,3
!           OrderOut_(i) = trim(OrderOut(i))
!        enddo
!     endif
!     !   
!     Nlso = size(Huser,1)
!     Nk   = size(Huser,3)
!     call assert_shape(Huser,[Nlso,Nlso,Nk],"tb_reorder_hk_d","Huser")
!     !
!     !Construct an index array InderOut corresponding to the Out ordering.
!     !This is a permutation of the In ordering taken as [1,2,3].
!     !For each entry in OrderIn we look for the position of the
!     !corresponding entry in OrderOut using tb_findloc.
!     !If 0 entries exist, corresponding components are not found. stop. 
!     do i=1,3     
!        IndexOut(i:i)=findloc_char(OrderIn,OrderOut_(i))
!     enddo
!     if(any(IndexOut==0))then
!        print*,"TB_Reorder_vec ERROR: wrong entry in IndexOut at: ",findloc_int(IndexOut,0)
!        stop
!     endif
!     write(*,*)IndexOut
!     !
!     !From IndexOut we can re-order the dimensions array to get the User dimensions array 
!     Nout=indx_reorder(Nin,IndexOut)
!     !    
!     if(any(IndexOut/=[1,2,3]))then
!        do iuser=1,Nlso
!           Ivec  = i2indices(iuser,Nin)           !Map iss to Ivec:(ilat,iorb,ispin) by IN ordering
!           Jvec  = indx_reorder(Ivec,IndexOut)  !Reorder according to Out ordering
!           iss   = indices2i(Jvec,Nout)         !Map back new Jvec to total index iuser by OUT ordering 
!           do juser=1,Nlso
!              Ivec  = i2indices(juser,Nin)
!              Jvec  = indx_reorder(Ivec,IndexOut)
!              jss = indices2i(Jvec,Nout)
!              Hss(iss,jss,:) = Huser(iuser,juser,:)
!           enddo
!        enddo
!     else
!        Hss = Huser
!     endif
!     return
!   end function reorder_hk


!   function indx_reorder(Ain,Index)  result(Aout)
!     integer,dimension(:)         :: Ain
!     integer,dimension(size(Ain)) :: Index
!     integer,dimension(size(Ain)) :: Aout
!     integer                        :: i
!     do i=1,size(Ain)
!        Aout(i) = Ain(Index(i))!Aout(Index(i)) = Ain(i)
!     enddo
!   end function indx_reorder


!   function indices2i(ivec,Nvec) result(istate)
!     integer,dimension(:)          :: ivec
!     integer,dimension(size(ivec)) :: Nvec
!     integer                       :: istate,i
!     istate=ivec(1)
!     do i=2,size(ivec)
!        istate = istate + (ivec(i)-1)*product(Nvec(1:i-1))
!     enddo
!   end function indices2i

!   function i2indices(istate,Nvec) result(ivec)
!     integer                       :: istate
!     integer,dimension(:)          :: Nvec
!     integer,dimension(size(Nvec)) :: Ivec
!     integer                       :: i,count,N
!     count = istate-1
!     N     = size(Nvec)
!     do i=1,N
!        Ivec(i) = mod(count,Nvec(i))+1
!        count   = count/Nvec(i)
!     enddo
!   end function i2indices



!   function findloc_char(array,val) result(pos)
!     character(len=*),dimension(:) :: array
!     character(len=*)              :: val
!     integer                       :: pos,i
!     pos=0
!     do i=1,size(array)
!        if(array(i)==val)then
!           pos = i
!           exit
!        endif
!     enddo
!     return
!   end function findloc_char

!   function findloc_int(array,val) result(pos)
!     integer,dimension(:) :: array
!     integer              :: val
!     integer              :: pos,i
!     pos=0
!     do i=1,size(array)
!        if(array(i)==val)then
!           pos = i
!           exit
!        endif
!     enddo
!     return
!   end function findloc_int
