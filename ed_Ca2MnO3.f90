!48x48 H_ij
! 2 Mn t_2g: [xz, yz, xy] - site 2 x orb 3 x spin 2: 12
! 6 Os p: x,y,z - site 6 x orb 3 x spin 2: 36
!
! NOTE ABOUT THE ORDER OF T_2g orbitals:
! In order to adhere to L*S representation we shall re-order
! the Hamiltonian H(k) obtained from W90 file to that the
! actual order will be
! [yz, xz, xy].
! This will be done just after the creation or reading of H(k).
program ed_CMO
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI
  USE MPI
#endif
  implicit none
  integer                                     :: Nlats(2),Nlat
  integer                                     :: Norbs(2),Nspins(2)
  integer                                     :: Nktot,Nkpath,Nkvec(3),Npts,Nlso,Ntot,Nso,Nbasis
  integer                                     :: iloop
  integer                                     :: i,j,k,ik,ilat,iorb,jorb,io,ispin
  integer                                     :: unit  
  real(8)                                     :: wmixing
  real(8),dimension(3)                        :: e1,e2,e3
  real(8),dimension(:,:),allocatable          :: kpath
  real(8)                                     :: ef,filling
  complex(8),dimension(:,:,:),allocatable     :: Hk,H0
  complex(8),dimension(:,:),allocatable       :: Hloc
  real(8),dimension(:),allocatable            :: Wtk,dens
  character(len=60)                           :: w90file,InputFile,latfile,kpathfile,hkfile
  character(len=40),allocatable               :: points_name(:)
  logical                                     :: converged
  logical                                     :: EFflag
  !Bath:
  integer                                     :: Nb
  real(8)   ,allocatable,dimension(:,:)       :: Bath
  !local dmft Weisss:
  complex(8),allocatable,dimension(:,:,:)     :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:)     :: Gmats,Greal
  complex(8),allocatable,dimension(:,:,:)     :: Weiss,Weiss_prev
  !Mpi:
  integer                                     :: comm,rank,ier
  logical                                     :: master=.true.,bool
  !logicals
  logical                                     :: bool_hk
  logical                                     :: bool_lat
  logical                                     :: bool_kpath
  !
  logical                                     :: full_solve
  !
  !Replica Hamiltonian
  real(8),dimension(:,:,:),allocatable        :: lambda_basis ![Nlat,Nbath,Nsym]
  complex(8),dimension(:,:,:,:,:),allocatable :: H_basis      ![size(Hloc),Nsym]
  !
  type(rgb_color),allocatable,dimension(:)    :: colors
  !
  complex(8),dimension(2,2)                   :: Sx,Sy,Sz
  complex(8),dimension(3,3)                   :: Lx,Ly,Lz,Lcf
  complex(8),dimension(6,6)                   :: Ix,Iy,Iz,Icf,Izf


  !SETUP MPI (under pre-processing instruction) 
#ifdef _MPI
  call init_MPI
  comm = MPI_COMM_WORLD
  master = get_master_MPI()
  rank = get_Rank_MPI()
#endif
  !
  !
  ! INTERNAL/LOCAL VARIABLE PARSINGe
  call parse_cmd_variable(InputFile,"INPUTFILE",default="inputDFT.conf")
  call parse_input_variable(full_solve,"FULL_SOLVE",InputFile,default=.false.,comment="Fix the type of solution/approximation to use. TRUE: full Hloc=H_cf_H_soc (mode=nonsu2 + bath=replica/general). FALSE  Hloc=H_cf, soc treated in the self-consistency, (mode=normal + bath=normal)")
  call parse_input_variable(Nkvec,"NKVEC",InputFile,default=[10,10,10])
  call parse_input_variable(nkpath,"NKPATH",InputFile,default=500)
  call parse_input_variable(wmixing,"WMIXING",InputFile,default=0.5d0)
  call parse_input_variable(w90file,"w90file",InputFile,default="hij.conf")
  call parse_input_variable(hkfile,"hkfile",InputFile,default="hk.conf")
  call parse_input_variable(latfile,"latfile",InputFile,default="lat.conf")
  call parse_input_variable(kpathfile,"kpathfile",InputFile,default="kpath.conf")
  call parse_input_variable(EFflag,"EFflag",InputFile,default=.false.)
  call parse_input_variable(filling,"filling",InputFile,default=1d0)
  call ed_read_input(reg(InputFile))
  !
  !PASS CONTROL VARIABLES TO DMFT_TOOLS USING +CTRL_VAR MEMORY POOL:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  !SETUP INTERNAL VARIABLES
  Nlats = [2,6]
  Norbs = [3,3]
  Nspins= [2,2]

  if(Norb/=Norbs(1))stop "ERROR: Norb != Norbs(1). Set the input file correctly" !Safety measure
  if(Nspin/=Nspins(1))stop "ERROR: Nspin != 2. Set the input file correctly"     !Safety measure
  Nlat = Nlats(1)
  Nso  = Nspin*Norb             !6 spin-orbit problem per site of d-electrons
  Nlso = Nlat*Nspin*Norb       !Total dimension of the correlated problem Nlat=2,Norb=3,Nspin=2  
  Ntot = sum(Nlats*Nspins*Norbs) !Total dimension of the whole problem (p+d_t2g)
  Nktot= product(Nkvec)         !Total dimension of the k-grid

  !DEFINE SPIN, ORBITAL AND TOTAL MOMENT MATRICES   
  Sx = pauli_sigma_x
  Sy = pauli_sigma_y
  Sz = pauli_sigma_z
  !
  Lx = reshape([zero,zero,zero,  zero,zero, -xi,  zero,  xi,zero],shape(Lx))  !Fortran column major
  Ly = reshape([zero,zero,  xi,  zero,zero,zero,   -xi,zero,zero],shape(Ly))  !Fortran column major
  Lz = reshape([zero, -xi,zero,    xi,zero,zero,  zero,zero,zero],shape(Lz))  !Fortran column major
  !
  Ix = kron(Sx,Lx)
  Iy = kron(Sy,Ly)
  Iz = kron(Sz,Lz)
  !
  !+ Crystal field basis diag(0,0,1) [(yz,xz) - xy + Delta]
  Lcf= reshape([zero,zero,zero,  zero,zero,zero,  zero,zero,one],shape(Lcf))
  Icf= kron(zeye(2),Icf)
  !
  !+ Zeeman separation:
  Izf= kron(Sz,zeye(3))



  !CHECK IF GIVEN FILES EXISTS
  inquire(file=reg(hkfile),exist=bool_hk)       !H(k) from previous calculation (save time)
  inquire(file=reg(latfile),exist=bool_lat)     !R-lattice basis coordinates (for bands plots)
  inquire(file=reg(kpathfile),exist=bool_kpath) !A list of high symmetry points in K-space



  !SETUP THE PATH IN THE BZ (for later plots)
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
     write(*,"(A)")"Using default path for 3d BZ: [M,R,\G,X,M,\G,Z,A,R]"
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




  !SET BASIS VECTORS:
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



  !SETUP WANNIER90 OR READ H(K) FROM FILE:                          
  call start_timer
  call TB_w90_setup(w90file,nlat=Nlats,norb=Norbs,nspin=Nspin,spinor=.true.,verbose=.true.)
  call stop_timer("TB_w90_setup")
  if(bool_hk)then
     call TB_read_hk(Hk,reg(hkfile),Nkvec)     
     call assert_shape(Hk,[Ntot,Ntot,product(Nkvec)])
  else
     call start_timer
     call TB_w90_FermiLevel(Nkvec,filling,Ef)
     call stop_timer("TB_w90_FermiLevel")
     !
     allocate(Hk(Ntot,Ntot,Nktot));Hk=zero
     call start_timer
     call TB_build_model(Hk,Ntot,Nkvec)
     !Reorder into [Norbs,Nspins,Nlats]->[d_orbs,d_spins,d_sites][p_orbs,p_spins,p_sites]->[Nlso][Ntot-Nlso]
     call TB_reshuffle_hk(Hk,[Nspins,Norbs,Nlats],[2,1,3])
     call TB_write_hk(Hk,reg(hkfile),Nkvec)
     call stop_timer("TB_build_model")
  endif

  !PERFORM ORDER EXCHANGE IN H(k) MATRIX
  ![xz, yz, xy] --> [yz, xz, xy]
  do ik=1,size(Hk,3)!==Nktot
     Hk(:,:,ik) = soc_reshape(Hk(:,:,ik)) !this should be safe in Fortran 
  enddo



  !GET THE LOCAL HAMILTONIAN 
  !Hloc_d1d2 = Sum_k H(k)[1:Nd1,1:Nd2] restricted to the correlated d-block
  allocate(Hloc(Ntot,Ntot))
  Hloc= sum(Hk,dim=3)/size(Hk,3)
  where(abs(Hloc)<1d-4)Hloc=zero
  call TB_write_Hloc(Hloc(1:6,1:6),"w90Hloc_d1.dat")
  call TB_write_Hloc(Hloc(7:12,7:12),"w90Hloc_d2.dat")


  !GET THE NON-INTERACTING BAND STRUCTURE
  call start_timer
  allocate(colors(Ntot))
  colors = gray
  colors(1:Nlso) = [[black,red,green,black,red,green],[black,red,green,black,red,green]]
  if(master)call TB_Solve_model(ed_Hk_model,Ntot,kpath,Nkpath,&
       colors_name=colors,&
       points_name=points_name,& 
       file="BandStructure",iproject=.false.)
  call stop_timer("TB get Bands")



  !CREATE THE REQUIRED DMFT ARRAYS FOR LOCAL FUNCTIONS:
  allocate(Smats(Ntot,Ntot,Lmats));      Smats=zero
  allocate(Gmats(Ntot,Ntot,Lmats));      Gmats=zero
  allocate(Sreal(Ntot,Ntot,Lreal));      Sreal=zero
  allocate(Greal(Ntot,Ntot,Lreal));      Greal=zero
  allocate(Weiss(Ntot,Ntot,Lmats));      Weiss=zero
  allocate(Weiss_prev(Ntot,Ntot,Lmats)); Weiss_prev=zero



  !SETUP THE LOCAL HAMILTONIAN AND THE BATH
  !ACCORDING TO THE VALUE OF FULL_SOLVE
  allocate(H0(Nlat,Nspin*Norb,Nspin*Norb))
  !
  select case(full_solve)
  case(.false.)
     ! use the H_cf in the local Hamiltonian
     ! the local SOC H_soc is treated within the self-consistency.
     ! In this case we can solve using ed_mode=normal.
     if(ed_mode/="normal")stop "ed_Ca2MnO3 ERROR: .not.full_solve but ed_mode!=normal"
     !
     !Get the diagonal of Hloc on both Mn sites in the unit cell
     H0(1,:,:) = diag(diagonal(Hloc(1:6,1:6)))
     H0(2,:,:) = diag(diagonal(Hloc(7:12,7:12)))
     !Set H0 as Hlocal in the impurity problem.  
     call ed_set_Hloc(H0,Nlat)
     !
  case (.true.)
     ! use the full local Hamiltonian H_loc = H_cf + H_soc.
     ! This case requires to have ed_mode=nonsu2 because H_soc breaks spin conservation and
     ! we require to use a replica bath, which is set here.
     if(ed_mode/="replica")stop "ed_Ca2MnO3 ERROR: full_solve but ed_mode!=replica"
     !
     Nbasis = 5 !3 S*L terms + Crystal field + Zeeman + do we need identity?
     allocate(H_basis(Nspin,Nspin,Norb,Norb,Nbasis))
     H_basis(:,:,:,:,1) = so2nn_reshape(Ix, Nspin,Norb)
     H_basis(:,:,:,:,2) = so2nn_reshape(Iy, Nspin,Norb)
     H_basis(:,:,:,:,3) = so2nn_reshape(Iz, Nspin,Norb)
     H_basis(:,:,:,:,4) = so2nn_reshape(Icf,Nspin,Norb)
     H_basis(:,:,:,:,5) = so2nn_reshape(Izf,Nspin,Norb)
     allocate(lambda_basis(Nlat,Nbath,Nbasis))
     !Mn_1
     lambda_basis(1,:,1) = 0.001d0  !tiny breaking along Ix
     lambda_basis(1,:,2) = 0.001d0  !tiny breaking along Iy
     lambda_basis(1,:,3) = 0.001d0  !tiny breaking along Iz
     lambda_basis(1,:,4) = 0.01d0   !tiny breaking along Icf
     lambda_basis(1,:,5) = 0.001d0  !tiny breaking along Izf
     !Mn_2
     lambda_basis(2,:,1) = 0.001d0  !tiny breaking along Ix
     lambda_basis(2,:,2) = 0.001d0  !tiny breaking along Iy
     lambda_basis(2,:,3) = 0.001d0  !tiny breaking along Iz
     lambda_basis(2,:,4) = 0.01d0   !tiny breaking along Icf
     lambda_basis(2,:,5) = -0.001d0 !tiny breaking along Izf
     !
     !SETUP H_replica = sum_{i=1,Nbasis} \lambda_i * H_i
     call ed_set_Hreplica(H_basis,lambda_basis)
     !
     !Set local Hamiltonian full symmetry:
     H0(1,:,:) = Hloc(1:6,1:6)
     H0(2,:,:) = Hloc(7:12,7:12)
     !Set H0 as Hlocal in the impurity problem.  
     call ed_set_Hloc(H0,Nlat)
  end select



  !INITIALIZE THE BATH AND THE SOLVER:
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nlat,Nb)); Bath=0d0
  call ed_init_solver(Bath)
  !

  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")

     !--------  solve impurity for the correlated orbitals/sites
     call ed_solve(Bath)


     !--------  retrieve sigmas and embed 
     Smats = zero
     call ed_get_sigma(Smats(1:Nlso,1:Nlso,:), Nlat, axis='matsubara')
     call dmft_write_gf(Smats,"Sigma",axis='mats',iprint=1)


     !------  get local Gfs
     call dmft_get_gloc(Hk,Gmats,Smats,axis='mats')
     call dmft_write_gf(Gmats,"Gloc",axis='mats',iprint=1)


     !------    get Weiss
     call dmft_self_consistency(Gmats,Smats,Weiss) 
     call dmft_write_gf(Weiss,"Weiss",axis='mats',iprint=1)


     !------    mix Weiss
     if(iloop>1)then
        Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_prev
     endif
     Weiss_prev=Weiss


     !------    fit Weiss     ------
     call ed_chi2_fitgf(Weiss(1:Nlso,1:Nlso,:),Bath)

     !Turn it off for the time being
     ! if(nread/=0d0)then
     !    call ed_get_dens(dens)
     !    call ed_search_variable(xmu,sum(dens),converged)
     ! endif
     !
     !Check convergence
     converged = check_convergence(Weiss(1,1,:),dmft_error,nsuccess,nloop)

     if(master)call end_loop
  enddo


  !Compute the local gfs:
  Sreal = zero
  call ed_get_sigma(Sreal(1:Nlso,1:Nlso,:), Nlat, axis='realaxis')
  call dmft_get_gloc(Hk,Greal,Sreal,axis='real')
  call dmft_write_gf(Greal,"Gloc",axis='real',iprint=1)

  ! !Compute the Kinetic Energy:
  ! call dmft_kinetic_energy(Hk,Smats)

  call finalize_MPI()


contains


  function ed_Hk_model(kvec,N) result(Hk)
    real(8),dimension(:)      :: kvec
    integer                   :: N
    complex(8),dimension(N,N) :: Hloc,Hk
    !< Build H(k) from W90:
    Hk = TB_w90_model(kvec,N)
    !< Reorder H(k) according to W90-SS orders
    call TB_reshuffle(Hk,[Nspins,Norbs,Nlats],[2,1,3])
    !< Build effective fermionic H*(k) 
    !>TODO
  end function ed_Hk_model




  !AUX ROUTINES FOR RESHAPE:
  function so2nn_reshape(Fin,Nspin,Norb) result(Fout)
    integer                                               :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb)           :: Fin
    complex(8),dimension(Nspin,Nspin,Norb,Norb)           :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Fout(ispin,jspin,iorb,jorb) = Fin(is,js)
             enddo
          enddo
       enddo
    enddo
  end function so2nn_reshape

  function nn2so_reshape(Fin,Nspin,Norb) result(Fout)
    integer                                               :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb)           :: Fin
    complex(8),dimension(Nspin*Norb,Nspin*Norb)           :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Fout(is,js) = Fin(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function nn2so_reshape


  function lso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: ilat
    integer                                               :: iorb,jorb
    integer                                               :: ispin,jspin
    integer                                               :: is,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function lso2nnn


  function nnn2lso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: ilat
    integer                                               :: iorb,jorb
    integer                                               :: ispin,jspin
    integer                                               :: is,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function nnn2lso




  function soc_reshape(Mat) result(Mtmp)
    complex(8),dimension(Nso,Nso) :: Mat
    complex(8),dimension(Nso,Nso) :: Mtmp
    integer                       :: i,j
    integer                       :: ispin,jspin
    integer                       :: iorb,jorb
    integer                       :: io,jo
    do iorb=1,Norb
       do ispin=1,Nspin
          i  = iorb + (ispin-1)*Norb
          io = soc_order(iorb) + (ispin-1)*Norb
          do jorb=1,Norb
             do jspin=1,Nspin
                j  = jorb + (jspin-1)*Norb
                jo = soc_order(jorb) + (jspin-1)*Norb
                Mtmp(i,j) = Mat(io,jo)
             enddo
          enddo
       enddo
    enddo
  end function soc_reshape

  function soc_order(i) result(j)
    integer :: i
    integer :: j
    j=i
    if(i<3)j=3-i
  end function soc_order



end program ed_CMO
