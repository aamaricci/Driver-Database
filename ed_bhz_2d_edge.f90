!!! ### ACTHUNG: post-processing has been partly done using the
!!! TightBinding code, updated to read and process SelfEnergy.
program ed_bhz_2d_edge
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none

  integer                                       :: iloop
  integer                                       :: Nlso
  integer                                       :: Nso
  integer                                       :: Nineq,Nlat
  integer                                       :: ilat,iy,iorb,ispin,ineq,i,print_mode
  logical                                       :: converged
  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:,:)            :: Bath_ineq
  real(8),allocatable,dimension(:,:)            :: Bath_prev
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal,Sreal_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal
  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hkr
  complex(8),allocatable,dimension(:,:)         :: bhzHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc,Hloc_ineq,S0

  !gamma matrices:
  complex(8),dimension(4,4)                     :: Gamma1,Gamma2,Gamma5,GammaN
  complex(8),dimension(4,4)                     :: GammaE0,GammaEx,GammaEy,GammaEz
  complex(8),dimension(4,4)                     :: GammaR0,GammaRx,GammaRy,GammaRz
  real(8),allocatable,dimension(:)              :: kxgrid
  real(8),dimension(:,:),allocatable            :: kpath
  integer                                       :: Nk,Ly,Nkpath
  real(8)                                       :: e0,mh,lambda,wmixing
  character(len=30)                             :: Params
  logical                                       :: spinsym,tridiag,lrsym,rebuild_sigma
  character(len=60)                             :: finput
  character(len=32)                             :: hkfile
  real(8),dimension(:,:),allocatable            :: Zmats
  complex(8),dimension(:,:,:),allocatable       :: Zfoo
  real(8),dimension(:,:),allocatable            :: lambdasym_vector
  complex(8),dimension(:,:,:,:,:),allocatable   :: Hsym_basis

  integer                                       :: comm,rank
  logical                                       :: master


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ_EDGE.conf')
  call parse_input_variable(Params,"Params",finput,default="E0EzEx",&
       comment="Ex; EzEx; E0Ex; E0Ez; E0EzEx")
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(Ly,"Ly",finput,default=20)
  call parse_input_variable(Nkpath,"NKPATH",finput,default=501)
  call parse_input_variable(tridiag,"TRIDIAG",finput,default=.true.)
  call parse_input_variable(mh,"MH",finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_input_variable(e0,"e0",finput,default=1d0)
  call parse_input_variable(lrsym,"LRSYM",finput,default=.true.)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  !
  call ed_read_input(trim(finput),comm)


  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  !set the global number of lattice sites equal to the number of layers along the y-axis
  Nlat = Ly
  Nineq= Ly
  if(lrsym)then
     if(mod(Ly,2)/=0)stop "Wrong setup from input file: Ly%2 > 0 (odd number of sites)"
     Nineq=Ly/2
     print*,"Using L-R Symmetry. Solve",Nineq," of",Nlat," sites."
     call sleep(2)
  endif

  !set the local number of total spin-orbitals (4)
  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso  = Nspin*Norb

  !set the total lattice-spin-orbit dimension:
  Nlso=Nlat*Nspin*Norb



  !Allocate Functions:
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb));Hloc=zero
  allocate(S0(Nlat,Nspin,Nspin,Norb,Norb));S0=zero
  allocate(Zmats(Nlso,Nlso));Zmats=eye(Nlso)
  allocate(Zfoo(Nlat,Nso,Nso));Zfoo=0d0
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Weiss_ineq=zero
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Smats_ineq=zero
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal));Sreal_ineq=zero
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Gmats_ineq=zero
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb));Hloc_ineq=zero

  !SETUP THE GAMMA MATRICES:
  gamma1=kron_pauli( pauli_sigma_z, pauli_tau_x)
  gamma2=kron_pauli( pauli_sigma_0,-pauli_tau_y)
  gamma5=kron_pauli( pauli_sigma_0, pauli_tau_z)
  gammaN=kron_pauli( pauli_sigma_0, pauli_tau_0)
  !
  gammaE0=kron_pauli( pauli_sigma_0, pauli_tau_x )
  gammaEx=kron_pauli( pauli_sigma_x, pauli_tau_x )
  gammaEy=kron_pauli( pauli_sigma_y, pauli_tau_x )
  gammaEz=kron_pauli( pauli_sigma_z, pauli_tau_x )

  print_mode=1
  if(ed_mode=="nonsu2")print_mode=4


  !Buil the Hamiltonian on a grid or on  path
  call build_hkr(trim(hkfile))
  Hloc = lso2nnn_reshape(bhzHloc,Nlat,Nspin,Norb)
  do ineq=1,Nineq
     ilat = ineq2ilat(ineq)
     Hloc_ineq(ineq,:,:,:,:) = Hloc(ilat,:,:,:,:)
  enddo


  !Setup solver
  !Setup solver
  if(bath_type=="replica")then
     select case(trim(Params))
     case default
        stop "Params not in [Ex; EzEx; E0Ex; ExEy; E0Ez; E0EzEx; E0EzExEy]"
     case("Ex")
        allocate(lambdasym_vector(Nineq,2))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,2))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(:,1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaEx) ;lambdasym_vector(:,2)=-sb_field

     case("EzEx")
        allocate(lambdasym_vector(Nineq,3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(:,1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaEz) ;lambdasym_vector(:,2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEx) ;lambdasym_vector(:,3)=-sb_field

     case("E0Ex")
        allocate(lambdasym_vector(Nineq,3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(:,1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0) ;lambdasym_vector(:,2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEx) ;lambdasym_vector(:,3)=-sb_field

     case("E0Ez")
        allocate(lambdasym_vector(Nineq,3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(:,1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0) ;lambdasym_vector(:,2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEz) ;lambdasym_vector(:,3)=sb_field

     case("E0EzEx")
        allocate(lambdasym_vector(Nineq,4))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,4))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(:,1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0) ;lambdasym_vector(:,2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEz) ;lambdasym_vector(:,3)=sb_field
        Hsym_basis(:,:,:,:,4)=j2so(GammaEx) ;lambdasym_vector(:,4)=-sb_field
     end select
     call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
     Nb=ed_get_bath_dimension(Hsym_basis)
  else     
     Nb=ed_get_bath_dimension()
  endif
  allocate(Bath_ineq(Nineq,Nb) )
  allocate(Bath_prev(Nineq,Nb) )
  call ed_init_solver(comm,Bath_ineq)




  !DMFT loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     ! solve the impurities on each inequivalent y-layer
     call ed_solve(comm,Bath_ineq,Hloc_ineq)
     ! store the 1st Matsubara freq. into S0, used to get H_topological = Hk + S0
     call ed_get_sigma_matsubara(Smats_ineq,Nineq)
     do ilat=1,Nlat
        ineq = ilat2ineq(ilat)
        Smats(ilat,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:)
        S0(ilat,:,:,:,:)      = Smats_ineq(ineq,:,:,:,:,1)
     enddo
     do ilat=1,Nlat
        Zfoo(ilat,:,:)        = select_block(ilat,S0)
        do iorb=1,Nso
           i = iorb + (ilat-1)*Nso
           Zmats(i,i)  = 1.d0/( 1.d0 + abs( dimag(Zfoo(ilat,iorb,iorb))/(pi/beta) ))
        enddo
     enddo
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(Hkr,Gmats,Smats,tridiag=tridiag)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=print_mode)


     !Update Weiss Field
     do ineq=1,Nineq
        ilat = ineq2ilat(ineq)
        Gmats_ineq(ineq,:,:,:,:,:) = Gmats(ilat,:,:,:,:,:)
     enddo
     !
     ! compute the Weiss field (only the Nineq ones)
     call dmft_self_consistency(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,cg_scheme)
     call dmft_print_gf_matsubara(Weiss_ineq,"Weiss",iprint=print_mode)
     !
     ! fit baths and mix result with old baths
     select case(ed_mode)
     case default
        stop "ed_mode!=Normal/Nonsu2"
     case("normal")
        call ed_chi2_fitgf(Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=1)
        call ed_spin_symmetrize_bath(Bath_ineq,save=.true.)
     case("nonsu2")
        call ed_chi2_fitgf(Bath_ineq,Weiss_ineq,Hloc_ineq)
     end select
     !

     if(iloop>1)Bath_ineq=wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath_ineq
     converged = check_convergence(Weiss_ineq(:,1,1,1,1,:),dmft_error,nsuccess,nloop)


     call end_loop
  enddo


  call ed_get_sigma_realaxis(Sreal_ineq,Nineq)
  do ilat=1,Nlat
     ineq = ilat2ineq(ilat)
     Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
  enddo
  call dmft_gloc_realaxis(Hkr,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=print_mode)

  call dmft_kinetic_energy(Hkr,Smats)


  call solve_hk_topological()


  call save_array("Smats",Smats)
  call save_array("Sreal",Sreal)

  call finalize_MPI()


contains



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: build the BHZ Hamiltonian H(k_x,R_y) on the STRIPE along Y
  !+-----------------------------------------------------------------------------+!
  subroutine build_hkr(file)
    character(len=*),optional          :: file
    integer :: i,ik
    !
    !
    !SETUP THE H(kx,Ry):
    if(master)then
       write(LOGfile,*)"Build H(kx,y) for BHZ-stripe:"
       write(*,*)"# of kx-points     :",Nk
       write(*,*)"# of y-layers      :",Nlat
    endif
    !
    if(allocated(Hkr))deallocate(Hkr)
    allocate(Hkr(Nlso,Nlso,Nk))
    !
    call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])
    call TB_build_model(Hkr,bhz_edge_model,Ly,Nso,[Nk,1,1],pbc=.false.)
    !
    !SETUP THE LOCAL PART Hloc(Ry)
    allocate(bhzHloc(Nlso,Nlso))
    bhzHloc = sum(Hkr,dim=3)/Nk
    where(abs(bhzHloc)<1d-6)bhzHloc=zero
  end subroutine build_hkr











  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve H_BHZ(k_x,R_y) along the 1d -pi:pi path in the BZ.
  !+-----------------------------------------------------------------------------+!  
  subroutine solve_hk_topological(kpath_)
    real(8),dimension(:,:),optional    :: kpath_
    real(8),dimension(:,:),allocatable :: kpath
    integer                            :: Npts
    character(len=64)                  :: file
    type(rgb_color),dimension(:,:),allocatable :: colors
    !
    if(present(kpath_))then
       if(master)write(LOGfile,*)"Solve H(kx,y) along a given path:"
       Npts = size(kpath_,1)
       allocate(kpath(Npts,size(kpath_,2)))
       kpath=kpath_
       file="Eigenbands_path.nint"
    else
       !PRINT H(kx,Ry) ALONG A -pi:pi PATH
       if(master)write(LOGfile,*)"Solve H(kx,y) along [-pi:pi]:"
       Npts=3
       allocate(Kpath(Npts,1))
       kpath(1,:)=[-1]*pi
       kpath(2,:)=[ 0]*pi
       kpath(3,:)=[ 1]*pi
       file="Eigenbands.nint"
    endif

    allocate(colors(Ly,Nso))
    colors = gray88
    colors(1,:) = [red1,gray88,blue1,gray88,blue1,gray88,red1,gray88]
    colors(Ly,:) =[blue1,gray88,red1,gray88,red1,gray88,blue1,gray88]

    call TB_solve_model(bhz_edge_model,Ly,Nso,kpath,Nkpath,&
         colors_name=colors,&
         points_name=[character(len=10) :: "-pi","0","pi"],&
         file="Eigenbands.nint",pbc=.false.)
  end subroutine solve_hk_topological







  !+-----------------------------------------------------------------------------+!
  !PURPOSE: the BHZ-edge model hamiltonian
  !+-----------------------------------------------------------------------------+!
  !BHZ on a stripe geometry;
  function bhz_edge_model(kpoint,Nlat,N,pbc) result(Hrk)
    real(8),dimension(:)                :: kpoint
    real(8)                             :: kx
    integer                             :: Nlat,N
    complex(8),dimension(N,N)           :: Hmat,Tmat,TmatH
    complex(8),dimension(Nlat*N,Nlat*N) :: Hrk
    integer                             :: i,Idmin,Idmax,Itmin,Itmax
    logical                             :: pbc
    kx=kpoint(1)
    Hrk=zero
    Hmat=h0_rk_bhz(kx,N)
    Tmat=t0_rk_bhz(N)
    TmatH=conjg(transpose(Tmat))
    do i=1,Nlat
       Idmin=1+(i-1)*N
       Idmax=      i*N
       Hrk(Idmin:Idmax,Idmin:Idmax)=Hmat + dreal(select_block(i,S0)) !< H(k) + Re(Sigma_iy(:Nso,:Nso;omega=0))
    enddo
    do i=1,Nlat-1
       Idmin=1 + (i-1)*N
       Idmax=        i*N
       Itmin=1 +     i*N
       Itmax=    (i+1)*N
       Hrk(Idmin:Idmax,Itmin:Itmax)=Tmat
       Hrk(Itmin:Itmax,Idmin:Idmax)=TmatH
    enddo
    if(pbc)then
       Itmin=1+(Nlat-1)*N
       Itmax=0+Nlat*N
       Hrk(1:N,Itmin:Itmax)=TmatH
       Hrk(Itmin:Itmax,1:N)=Tmat
    endif
    Hrk = matmul(Zmats,Hrk)
  end function bhz_edge_model

  function h0_rk_bhz(kx,N) result(H)
    real(8)                    :: kx
    integer                    :: N
    complex(8),dimension(N,N)  :: H
    H = (mh-e0*cos(kx))*gamma5 + lambda*sin(kx)*gamma1
  end function h0_rk_bhz

  function t0_rk_bhz(N) result(H)
    integer                    :: N
    complex(8),dimension(N,N)  :: H
    H = -0.5d0*e0*gamma5 + xi*0.5d0*lambda*gamma2
  end function T0_rk_bhz





  function ilat2ineq(ilat) result(ineq)
    integer,intent(in) :: ilat
    integer            :: ineq
    ineq=ilat
    if( lrsym .AND. (ilat>Nineq) )ineq=Nlat-ilat+1
  end function ilat2ineq

  function ineq2ilat(ineq) result(ilat)
    integer,intent(in) :: ineq
    integer            :: ilat
    ilat=ineq
    if(ineq>Nineq)stop "ineq2ilat error: called with ineq > Nineq"
  end function ineq2ilat




  function select_block(ip,Matrix) result(Vblock)
    integer                                          :: ip
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Matrix
    complex(8),dimension(Nspin*Norb,Nspin*Norb)      :: Vblock
    integer                                          :: is,js,ispin,jspin,iorb,jorb
    Vblock=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Vblock(is,js) = Matrix(ip,ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function select_block



  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop "error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index


  function so2j(fg) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function so2j

  function j2so(fg) result(g)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2so

  function lso2nnn_reshape(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function lso2nnn_reshape

end program ed_bhz_2d_edge
