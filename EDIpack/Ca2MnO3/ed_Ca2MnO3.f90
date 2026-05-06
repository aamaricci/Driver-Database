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
  integer                                  :: Nlats(2),Nlat
  integer                                  :: Norbs(2),Nspins(2)
  integer                                  :: Nktot,Nkpath,Nkvec(3),Npts,Nlso,Ntot,Nso,Nbasis
  integer                                  :: iloop
  integer                                  :: i,j,k,ik,ilat,iorb,jorb,io,ispin
  integer                                  :: unit  
  real(8)                                  :: wmixing
  real(8),dimension(3)                     :: e1,e2,e3
  real(8),dimension(:,:),allocatable       :: kpath
  real(8)                                  :: ef,filling
  complex(8),dimension(:,:,:),allocatable  :: Hk,H0
  complex(8),dimension(:,:),allocatable    :: Hloc,sigH,zetaH
  real(8),dimension(:),allocatable         :: Wtk,dens
  character(len=60)                        :: w90file,InputFile,latfile,kpathfile,hkfile
  character(len=40),allocatable            :: points_name(:)
  logical                                  :: converged
  logical                                  :: EFflag
  !Bath:
  integer                                  :: Nb
  real(8)   ,allocatable,dimension(:,:)    :: Bath
  !local dmft Weisss:
  complex(8),allocatable,dimension(:,:,:)  :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:)  :: Gmats,Greal
  complex(8),allocatable,dimension(:,:,:)  :: Weiss,Weiss_prev
  !Mpi:
  integer                                  :: comm,rank,ier
  logical                                  :: master=.true.,bool
  !logicals
  logical                                  :: bool_hk
  logical                                  :: bool_lat
  logical                                  :: bool_kpath
  !
  !Replica Hamiltonian
  real(8),dimension(:,:,:),allocatable     :: lambda_basis ![Nlat,Nbath,Nsym]
  complex(8),dimension(:,:,:),allocatable  :: H_basis      ![size(Hloc),Nsym]
  !
  type(rgb_color),allocatable,dimension(:) :: colors
  !
  complex(8),dimension(2,2)                :: Sx,Sy,Sz
  complex(8),dimension(3,3)                :: Lx,Ly,Lz,Lcf
  complex(8),dimension(6,6)                :: Ix,Iy,Iz,Icf,Izf


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
  Lcf= reshape([zero,zero,zero,  zero,zero,zero,  zero,zero,one],[3,3])
  Icf= kron(zeye(2),Lcf)
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



  !SETUP WANNIER90 
  call start_timer
  call TB_w90_setup(w90file,nlat=Nlats,norb=Norbs,nspin=Nspin,spinor=.true.,verbose=.true.)
  !PERFORM ORDER EXCHANGE: [xz, yz, xy] --> [yz, xz, xy]
  call TB_w90_Transform(t2g_order_exchange) 
  call stop_timer("TB_w90_setup")

  !GENERATE OR READ H(K):
  if(bool_hk)then
     call TB_read_hk(Hk,reg(hkfile),Nkvec)     
     call assert_shape(Hk,[Ntot,Ntot,product(Nkvec)])
  else
     if(EFflag)then
        call start_timer
        call TB_w90_FermiLevel(Nkvec,filling,Ef)
        call stop_timer("TB_w90_FermiLevel")
     endif
     !
     allocate(Hk(Ntot,Ntot,Nktot));Hk=zero
     call start_timer
     call TB_build_model(Hk,Ntot,Nkvec)
     !Reorder degrees of freedom into to suit DMFT ordering 
     ![Norbs,Nspins,Nlats]->[d_orbs,d_spins,d_sites][p_orbs,p_spins,p_sites]->[Nlso][Ntot-Nlso]
     call TB_reshuffle_hk(Hk,[Nspins,Norbs,Nlats],[2,1,3])
     call TB_write_hk(Hk,reg(hkfile),Nkvec)
     call stop_timer("TB_build_model")
  endif



  !GET THE LOCAL HAMILTONIAN 
  !Hloc_d1d2 = Sum_k H(k)[1:Nd1,1:Nd2] restricted to the correlated d-block
  allocate(sigH(Ntot,Ntot)) ;sigH=zero
  allocate(zetaH(Ntot,Ntot));zetaH=zero
  allocate(Hloc(Ntot,Ntot)) ;Hloc=zero
  Hloc= sum(Hk,dim=3)/size(Hk,3)
  where(abs(Hloc)<1d-3)Hloc=zero
  call TB_write_Hloc(Hloc(1:6,1:6),"w90Hloc_d1.dat")
  call TB_write_Hloc(Hloc(7:12,7:12),"w90Hloc_d2.dat")


  !GET THE NON-INTERACTING BAND STRUCTURE
  call  get_band_structure()


  !CREATE THE REQUIRED DMFT ARRAYS FOR LOCAL FUNCTIONS:
  allocate(Smats(Ntot,Ntot,Lmats));      Smats=zero
  allocate(Gmats(Ntot,Ntot,Lmats));      Gmats=zero
  allocate(Sreal(Ntot,Ntot,Lreal));      Sreal=zero
  allocate(Greal(Ntot,Ntot,Lreal));      Greal=zero
  allocate(Weiss(Ntot,Ntot,Lmats));      Weiss=zero
  allocate(Weiss_prev(Ntot,Ntot,Lmats)); Weiss_prev=zero



  !SETUP THE LOCAL HAMILTONIAN AND THE BATH
  allocate(H0(Nlat,Nso,Nso))
  !
  select case(ed_mode)
  case default
     stop "ed_Ca2MnO3 error: ed_mode is neither normal nor nonsu2."
  case("normal")
     ! use the H_cf in the local Hamiltonian
     ! the local SOC H_soc is treated within the self-consistency.
     ! In this case we can solve using ed_mode=normal.
     if(bath_type/="normal")stop "ed_Ca2MnO3 ERROR: bath_type!=normal"
     !
     !Get the diagonal of Hloc on both Mn sites in the unit cell
     H0(1,:,:) = diag(diagonal(Hloc(1:6,1:6)))
     H0(2,:,:) = diag(diagonal(Hloc(7:12,7:12)))
     !Set H0 as Hlocal in the impurity problem.  
     call ed_set_Hloc(H0,Nlat)
     Nb=ed_get_bath_dimension()
     !
  case ("nonsu2")
     ! use the full local Hamiltonian H_loc = H_cf + H_soc.
     ! This case requires to have ed_mode=nonsu2 because H_soc breaks spin conservation and
     ! we require to use a replica bath, which is set here.
     if(bath_type/="replica")stop "ed_Ca2MnO3 ERROR: bath_type!=replica"
     !
     !SETUP H_replica = sum_{i=1,Nbasis} \lambda_i * H_i
     Nbasis = 5 !3 S*L terms + Crystal field + Zeeman + do we need identity?
     allocate(H_basis(Nso,Nso,Nbasis))
     H_basis(:,:,1) = Ix
     H_basis(:,:,2) = Iy
     H_basis(:,:,3) = Iz
     H_basis(:,:,4) = Icf
     H_basis(:,:,5) = Izf
     allocate(lambda_basis(Nlat,Nbath,Nbasis))
     !Mn_1
     lambda_basis(1,:,1) = 0.001d0  !tiny breaking along Ix
     lambda_basis(1,:,2) = 0.001d0  !tiny breaking along Iy
     lambda_basis(1,:,3) = 0.001d0  !tiny breaking along Iz
     lambda_basis(1,:,4) = 0.001d0   !tiny breaking along Icf
     lambda_basis(1,:,5) = 0.01d0  !tiny breaking along Izf
     !Mn_2
     lambda_basis(2,:,1) = 0.001d0  !tiny breaking along Ix
     lambda_basis(2,:,2) = 0.001d0  !tiny breaking along Iy
     lambda_basis(2,:,3) = 0.001d0  !tiny breaking along Iz
     lambda_basis(2,:,4) = 0.001d0   !tiny breaking along Icf
     lambda_basis(2,:,5) = -0.01d0 !tiny breaking along Izf
     !     
     call ed_set_Hreplica(H_basis,lambda_basis)
     !
     !Set local Hamiltonian full symmetry:
     H0(1,:,:) = Hloc(1:6,1:6)
     H0(2,:,:) = Hloc(7:12,7:12)
     !
     !Set H0 as Hlocal in the impurity problem.  
     call ed_set_Hloc(H0,Nlat)
     Nb=ed_get_bath_dimension(5)
  end select

  !INITIALIZE THE BATH AND THE SOLVER:
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
     Sreal = zero
     call ed_get_sigma(Smats(1:Nlso,1:Nlso,:), Nlat, axis='matsubara')
     call ed_get_sigma(Sreal(1:Nlso,1:Nlso,:), Nlat, axis='realaxis')

     !------  get local Gfs
     call dmft_get_gloc(Hk,Gmats,Smats,axis='mats')

     !------    get Weiss
     call dmft_self_consistency(Gmats,Smats,Weiss) !LOCAL: calG0_ii = {[Gloc^-1]_ii + Sigma_ii}^-1 ii=Mn1,Mn2


     !------ write them all
     call dmft_write_gf(Smats(1:Nlso,1:Nlso,:),"Sigma",axis='mats',iprint=1,&
          Nvec=[1,1,1,1,1,1,1,1,1,1,1,1],&
          labels=[&
          "Mn1_dyz_up","Mn1_dxz_up","Mn1_dxy_up",&
          "Mn1_dyz_dw","Mn1_dxz_dw","Mn1_dxy_dw",&
          "Mn2_dyz_up","Mn2_dxz_up","Mn2_dxy_up",&
          "Mn2_dyz_dw","Mn2_dxz_dw","Mn2_dxy_dw"])

     call dmft_write_gf(Gmats,"Gloc",axis='mats',iprint=1,&
          Nvec=[1,1,1,1,1,1,1,1,1,1,1,1,Ntot-Nlso],&
          labels=[&
          "Mn1_dyz_up","Mn1_dxz_up","Mn1_dxy_up",&
          "Mn1_dyz_dw","Mn1_dxz_dw","Mn1_dxy_dw",&
          "Mn2_dyz_up","Mn2_dxz_up","Mn2_dxy_up",&
          "Mn2_dyz_dw","Mn2_dxz_dw","Mn2_dxy_dw",&
          "p         "])

     call dmft_write_gf(Weiss,"Weiss",axis='mats',iprint=1,&
          Nvec=[1,1,1,1,1,1,1,1,1,1,1,1,Ntot-Nlso],&
          labels=[&
          "Mn1_dyz_up","Mn1_dxz_up","Mn1_dxy_up",&
          "Mn1_dyz_dw","Mn1_dxz_dw","Mn1_dxy_dw",&
          "Mn2_dyz_up","Mn2_dxz_up","Mn2_dxy_up",&
          "Mn2_dyz_dw","Mn2_dxz_dw","Mn2_dxy_dw",&
          "p         "])



     !------    mix Weiss
     if(iloop>1)then
        Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_prev
     endif
     Weiss_prev=Weiss


     !------    fit Weiss     ------
     select case(ed_mode)
     case default
        stop "ed_mode!=Normal/Nonsu2"
     case("normal")
        call ed_chi2_fitgf(Weiss(1:Nlso,1:Nlso,:),Bath,ispin=1)
        call ed_chi2_fitgf(Weiss(1:Nlso,1:Nlso,:),Bath,ispin=2)
     case("nonsu2")
        call ed_chi2_fitgf(Weiss(1:Nlso,1:Nlso,:),Bath)
     end select

     !if target density is given (nread) chemical potential is adjusted
     !This needs to be tuned better to comply with P-D structure.
     !SO far it fixes total density of the D bands on Mn's
     if(nread/=0d0)then
        call ed_get_dens(dens)
        !get P-density sum to D-dens
        call ed_search_variable(xmu,sum(dens),converged)
     endif
     !
     !Check convergence
     converged = check_convergence(Weiss(1,1,:),dmft_error,nsuccess,nloop)

     if(master)call end_loop
  enddo


  !COMPUTE THE LOCAL REAL-AXIS GF::
  Sreal = zero
  call ed_get_sigma(Sreal(1:Nlso,1:Nlso,:), Nlat, axis='realaxis')
  call dmft_get_gloc(Hk,Greal,Sreal,axis='real')
  call dmft_write_gf(Sreal(1:Nlso,1:Nlso,:),"Sigma",axis='real',iprint=1,&
       Nvec=[1,1,1,1,1,1,1,1,1,1,1,1],&
       labels=[&
       "Mn1_dyz_up","Mn1_dxz_up","Mn1_dxy_up",&
       "Mn1_dyz_dw","Mn1_dxz_dw","Mn1_dxy_dw",&
       "Mn2_dyz_up","Mn2_dxz_up","Mn2_dxy_up",&
       "Mn2_dyz_dw","Mn2_dxz_dw","Mn2_dxy_dw"])

  call dmft_write_gf(Greal,"Gloc",axis='real',iprint=1,&
       Nvec=[1,1,1,1,1,1,1,1,1,1,1,1,Ntot-Nlso],&
       labels=[&
       "Mn1_dyz_up","Mn1_dxz_up","Mn1_dxy_up",&
       "Mn1_dyz_dw","Mn1_dxz_dw","Mn1_dxy_dw",&
       "Mn2_dyz_up","Mn2_dxz_up","Mn2_dxy_up",&
       "Mn2_dyz_dw","Mn2_dxz_dw","Mn2_dxy_dw",&
       "p         "])


  !COMPUTE THE KINETIC ENERGY:
  call dmft_kinetic_energy(Hk,Smats)

  !GET THE INTERACTING BAND STRUCTURE
  call  get_band_structure(Smats(:,:,1))


  call save_array("Smats",Smats)
  call save_array("Sreal",Sreal)


  call finalize_MPI()


contains



  !get the system band structure, renormalized or bare
  subroutine get_band_structure(sigma)
    complex(8),dimension(Ntot,Ntot),optional :: sigma !Sigma(iw_1)
    character(len=32) :: file_name
    call start_timer
    call set_sigmaH();if(present(sigma))call set_sigmaH(sigma)
    file_name="BandStructure";if(present(sigma))file_name="Htop_BandStructure"
    if(.not.allocated(colors))allocate(colors(Ntot))
    colors = gray
    colors(1:Nlso) = [[black,red,green,black,red,green],[black,red,green,black,red,green]]
    if(master)call TB_Solve_model(ed_Hk_model,Ntot,kpath,Nkpath,&
         colors_name=colors,&
         points_name=points_name,& 
         file=str(file_name),iproject=.false.)
    call stop_timer("TB get Bands")
  end subroutine get_band_structure



  !Set the Self-Energy for renormalized band structure calculation
  subroutine set_sigmaH(sigma)
    complex(8),dimension(Ntot,Ntot),optional :: sigma
    integer                                  :: ii
    sigH  = zero
    zetaH = zeye(Ntot) 
    if(present(sigma))then
       sigH = dreal(sigma)
       zetaH=zero
       do ii=1,Ntot
          zetaH(ii,ii)  = one/abs( one +  abs(dimag(sigma(ii,ii))/(pi/beta)) )
       end do
    endif
  end subroutine set_sigmaH


  !Renormalized H(k) model: G^-1(k,w=0) = Z*(H(k)+ReSigma(w=0))
  function ed_Hk_model(kvec,N) result(Hk)
    real(8),dimension(:)      :: kvec
    integer                   :: N
    complex(8),dimension(N,N) :: Hk
    !< Build H(k) from W90:
    Hk = TB_w90_model(kvec,N)
    !< Reorder H(k) according to W90-SS orders
    call TB_reshuffle(Hk,[Nspins,Norbs,Nlats],[2,1,3])
    !< Build effective fermionic H*(k)
    ! !add the SigmaBHZ term to get Topologial Hamiltonian if required:
    Hk = Hk + dreal(sigH)
    Hk = matmul(zetaH,Hk)
  end function ed_Hk_model





  !Function to reorder [xz,yz,xy]-->[yz,xz,xy]
  function t2g_order_exchange(Hin,N) result(Hout)
    integer                   :: N
    complex(8),dimension(N,N) :: Hin
    complex(8),dimension(N,N) :: Hout
    integer                   :: i,j,ii,jj
    do i=1,N
       ii = swap_indices(i)
       do j=1,N
          jj = swap_indices(j)
          Hout(ii,jj) = Hin(i,j)
       enddo
    enddo
  end function t2g_order_exchange
  !
  function swap_indices(i) result(j)
    integer :: i
    integer :: j
    select case (i)
    case default ; j=i          !
    case(1)      ; j=3          !xz_up_Mn1:1-->yz_up_Mn1:3 
    case(2)      ; j=4          !xz_dw_Mn1:2-->yz_dw_Mn1:4
    case(3)      ; j=1          !vice-versa
    case(4)      ; j=2          !
    case(7)      ; j=9          !xz_up_Mn2:7-->yz_up_Mn2:9 
    case(8)      ; j=10         !xz_dw_Mn2:8-->yz_dw_Mn2:10
    case(9)      ; j=7          !vice-versa
    case(10)     ; j=8          !
    end select
  end function swap_indices




  !THIS IS NOT USED ANYMORE HERE, LEFT FOR LEGACY...
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

  function lso2nnn_iw(Hlso,Nlat,Nspin,Norb,L) result(Hnnn)
    integer                                                 :: Nlat,Nspin,Norb,L
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,L) :: Hlso
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,L)      :: Hnnn
    integer                                                 :: iorb,ispin,ilat,is
    integer                                                 :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb,:) = Hlso(is,js,:)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function lso2nnn_iw

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











  subroutine write_matrix_i(a)
    integer,dimension(:,:) :: a
    write(*,*)
    do i = lbound(a,1), ubound(a,1)
       write(*,*) (a(i,j), j = lbound(a,2), ubound(a,2))
    end do
  end subroutine write_matrix_i

  subroutine write_matrix_d(a)
    real(8),dimension(:,:) :: a
    write(*,*)
    do i = lbound(a,1), ubound(a,1)
       write(*,*) (str(a(i,j)), j = lbound(a,2), ubound(a,2))
    end do
  end subroutine write_matrix_d

  subroutine write_matrix_c(a)
    complex(8),dimension(:,:) :: a
    write(*,*)
    do i = lbound(a,1), ubound(a,1)
       write(*,"(100(A1,F4.1,1x,F4.1,A1,2X))") ('(',a(i,j),')', j = lbound(a,2), ubound(a,2))
    end do
  end subroutine write_matrix_c


end program ed_CMO
