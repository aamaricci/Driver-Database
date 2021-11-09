program ed_W90
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI
  USE MPI
#endif
  implicit none
  integer                                     :: Nlat,Nktot,Nkpath,Nkvec(3),Npts,Nlso
  integer                                     :: iloop
  integer                                     :: i,j,k,ik,ilat,iorb,jorb,io,ispin
  integer                                     :: unit  
  real(8)                                     :: wmixing
  real(8),dimension(3)                        :: e1,e2,e3
  real(8),dimension(:,:),allocatable          :: kpath
  real(8)                                     :: ef,filling
  complex(8),dimension(:,:,:),allocatable     :: Hk
  complex(8),dimension(:,:),allocatable       :: Hloc
  real(8),dimension(:),allocatable            :: Wtk,dens
  character(len=60)                           :: w90file,InputFile,latfile,kpathfile,hkfile
  character(len=40),allocatable               :: points_name(:)
  logical                                     :: converged
  logical                                     :: EFflag
  !Bath:
  integer                                     :: Nb
  real(8)   ,allocatable,dimension(:)         :: Bath
  !local dmft Weisss:
  complex(8),allocatable,dimension(:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats,Greal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Weiss_prev
  !Mpi:
  integer                                     :: comm,rank,ier
  logical                                     :: master=.true.,bool
  logical                                     :: bool_hk
  logical                                     :: bool_lat
  logical                                     :: bool_kpath



  !
#ifdef _MPI
  call init_MPI
  comm = MPI_COMM_WORLD
  master = get_master_MPI()
  rank = get_Rank_MPI()
#endif
  !
  !
  !#########    VARIABLE PARSING    #########
  !
  call parse_cmd_variable(InputFile,"INPUTFILE",default="inputDFT.conf")
  ! call parse_input_variable(NLAT,"NLAT",InputFile,default=1)
  call parse_input_variable(Nkvec,"NKVEC",InputFile,default=[10,10,10])
  call parse_input_variable(nkpath,"NKPATH",InputFile,default=500)
  call parse_input_variable(wmixing,"WMIXING",InputFile,default=0.5d0)
  call parse_input_variable(w90file,"w90file",InputFile,default="hij.conf")
  call parse_input_variable(hkfile,"hkfile",InputFile,default="hk.conf")
  call parse_input_variable(latfile,"latfile",InputFile,default="lat.conf")
  call parse_input_variable(kpathfile,"kpathfile",InputFile,default="kpath.conf")
  call parse_input_variable(EFflag,"EFflag",InputFile,default=.false.)
  call parse_input_variable(filling,"filling",InputFile,default=1d0)
  call ed_read_input(reg(InputFile),comm)
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  !>DEVEL
  Nlat=1
  if(Nlat/=1)stop "This driver is now only for Nlat=1"
  !<DEVEL
  Nlso = Nlat*Nspin*Norb
  Nktot= product(Nkvec)
  
  !METHOD 1 (setup W90 --> use internal W90 model)
  inquire(file=reg(hkfile),exist=bool_hk)
  inquire(file=reg(latfile),exist=bool_lat)
  inquire(file=reg(kpathfile),exist=bool_kpath)


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


  !Setup Wannier90 or read H(k) from file:                          
  call start_timer
  call TB_w90_setup(reg(w90file),nlat=Nlat,nspin=Nspin,norb=Norb,verbose=.true.)
  call stop_timer("TB_w90_setup")


  if(bool_hk)then
     call TB_read_hk(Hk,reg(hkfile),Nlat,Nspin,Norb,Nkvec)     
     call assert_shape(Hk,[Nlso,Nlso,product(Nkvec)])
  else
     call start_timer
     call TB_w90_FermiLevel(Nkvec,filling,Ef)
     call stop_timer("TB_w90_FermiLevel")
     !
     allocate(Hk(Nlso,Nlso,Nktot));Hk=zero
     call start_timer
     call TB_build_model(Hk,Nlso,Nkvec)
     call TB_write_hk(Hk,reg(hkfile),Nlat,Nspin,Norb,Nkvec)
     call stop_timer("TB_build_model")
  endif

  allocate(Hloc(Nlso,Nlso))
  Hloc= sum(Hk(:,:,:),dim=3)/size(Hk,3)
  where(abs(Hloc)<1d-6)Hloc=zero
  if(master)call TB_write_Hloc(Hloc,"w90Hloc.dat")
  write(*,*)"Using Nk_total="//str(size(Hk,3))

  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));      Smats=zero
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));      Gmats=zero
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal));      Sreal=zero
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal));      Greal=zero
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats));      Weiss=zero
  allocate(Weiss_prev(Nspin,Nspin,Norb,Norb,Lmats)); Weiss_prev=zero
  allocate(dens(Norb))

  !If bath_type /= replica: 1.Hloc should not be passed here
  !If bath_type == replica: 2.if(using Hloc to define each replica) one must pass Hloc here
  !                         3.if(using matrix basis to define each replica) one should not pass Hloc here
  !So far we avoid condition 2.
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nb));    Bath=0.0d0
  call ed_init_solver(comm,Bath) !ed_init_solver DOES NOT accept Hloc as input anymore.
  !  

  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")

     !--------  solve impurity (diagonal basis)  --------
     call ed_solve(comm,Bath,so2nn(Hloc)) !ed_solve REQUIRED Hloc as input by now.
     !
     !--------    get sigmas   (diagonal basis)  --------
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)

     !------  get local Gf's  (Wannier90 basis)  ------
     call dmft_gloc_matsubara(Hk,Gmats,Smats) !dmft_gloc_xyz DOES NOT accept wtk as input anymore.
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)

     !------    get Weiss     (Wannier90 basis)  --------
     call dmft_self_consistency(Gmats,Smats,Weiss,so2nn(Hloc),cg_scheme)
     if(master)call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)


     !------    mix Weiss     ------
     if(iloop>1)then
        Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_prev
     endif
     Weiss_prev=Weiss
     !
     !------    fit Weiss     ------
     call ed_chi2_fitgf(comm,Weiss,Bath,ispin=1)
     !
     !All routines here are MPI aware, so no if(master) is required anymore.
     if(nread/=0d0)then
        call ed_get_dens(dens)
        call ed_search_variable(xmu,sum(dens),converged) !ed_search_variable uses Luca's algorithm, though a bit experimental here
        !call ed_search_chemical_potential(xmu,sum(dens),converged) !this is the other algorithm, slower.
     endif
     !
     !Check convergence (if required change chemical potential)
     converged = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop)

     !
     if(master)call end_loop
     !
  enddo


  !Compute the local gfs:
  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)

  !Compute the Kinetic Energy:
  call dmft_kinetic_energy(Hk,Smats)

  ! !Rebuild required Gimp or Self-energy. Use latest dmft_print_function to print.
  ! !Matsubara example:
  ! beta=beta*10d0 !change Mats freq mesh
  ! Lf  = 4096     !use arbitray freq. number
  ! allocate(wfreq(Lf))
  ! allocate(Greb(Nspin,Nspin,Norb,Norb,Lf)) 
  ! allocate(Sreb(Nspin,Nspin,Norb,Norb,Lf))
  ! wfreq=pi/beta*(2*arange(1,Lf)-1) !Matsubara freq.
  ! !call ed_get_gimp/sigma, specifify Imaginary Freq. == Matsubara
  ! call ed_get_gimp(dcmplx(0d0,wfreq),Hloc,Greb)
  ! call ed_get_sigma(dcmplx(0d0,wfreq),Hloc,Sreb)
  ! call dmft_print_function(Greb,"reGimp",iprint=1,axis="mats",zeta=wfreq)
  ! call dmft_print_function(Sreb,"reSigma",iprint=1,axis="mats",zeta=wfreq)
  ! !
  ! deallocate(wfreq,Greb,Sreb)
  ! !
  ! !Realaxis example  
  ! Lf = 10000                    !arbitrary number of freq.
  ! eps= 0.0099314d0              !arbitary imaginary shift
  ! allocate(wfreq(Lf))
  ! allocate(Greb(Nspin,Nspin,Norb,Norb,Lf)) 
  ! allocate(Sreb(Nspin,Nspin,Norb,Norb,Lf))
  ! wfreq=linspace(-10d0,pi2,Lf)  !arbitrary domain
  ! !call ed_get_gimp/sigma, specifify Complex Freq. == realaxis+xi*eps
  ! call ed_get_gimp(dcmplx(wfreq,eps),Hloc,Greb)
  ! call ed_get_sigma(dcmplx(wfreq,eps),Hloc,Sreb)
  ! call dmft_print_function(Greb,"reGimp",iprint=1,axis="real",zeta=wfreq)
  ! call dmft_print_function(Sreb,"reSigma",iprint=1,axis="real",zeta=wfreq)

  call finalize_MPI()

contains


  function so2nn(Hlso) result(Hnnn)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hlso
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnnn
    integer                                     :: iorb,jorb
    integer                                     :: ispin,jspin
    integer                                     :: is,js
    Hnnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb
                js = jorb + (jspin-1)*Norb
                Hnnn(ispin,jspin,iorb,jorb) = Hlso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function so2nn


  function nn2so(Hnnn) result(Hlso)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnnn
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hlso
    integer                                     :: iorb,jorb
    integer                                     :: ispin,jspin
    integer                                     :: is,js
    Hlso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb
                js = jorb + (jspin-1)*Norb
                Hlso(is,js) = Hnnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function nn2so

end program ed_W90
