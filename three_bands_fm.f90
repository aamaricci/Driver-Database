program ed_ti_slab
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none

  integer                                       :: iloop
  integer                                       :: Nineq,Nlat,Nlso,Nso
  integer                                       :: ilat,iorb,ispin,ineq
  logical                                       :: converged,converged_xmu,PBC
  
  !Bath:
  integer                                       :: Nb
  real(8)                                       :: wmixing
  real(8),allocatable,dimension(:,:)            :: Bath_ineq
  real(8),allocatable,dimension(:,:)            :: Bath_prev
  
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal,Sreal_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: S0
  
  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hkr
  complex(8),allocatable,dimension(:,:)         :: tiHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc,Hloc_ineq
  real(8),dimension(:,:),allocatable            :: kpath
  integer                                       :: Nk,Lk,Ly,Nkpath
  real(8)                                       :: dens

  !gamma matrices and fields:
  complex(8),dimension(4,4)                     :: emat,soxmat,soymat,sozmat
  real(8)                                       :: e0,mh,lambda
  complex(8),dimension(4,4)                     :: NHmat,GapOpeningMat
  real(8)                                       :: GapOpeningField,NHfield_up,NHfield_dw

  !misc
  real                                          :: XMU_OLD, XMU_REAL
  logical                                       :: tridiag,lrsym
  character(len=60)                             :: finput
  character(len=32)                             :: hkfile
  complex(8),dimension(:,:,:),allocatable       :: toconverge

  !mpi
  integer                                       :: comm,rank
  logical                                       :: master,getbands,Refine_xmu


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  ! 
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED_ti_SLAB.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(Ly,"Ly",finput,default=20)
  call parse_input_variable(Nkpath,"NKPATH",finput,default=501)
  call parse_input_variable(tridiag,"TRIDIAG",finput,default=.true.)
  call parse_input_variable(mh,"MH",finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_input_variable(GapOpeningField,"GAPOPENINGFIELD",finput,default=0.0d0)
  call parse_input_variable(NHfield_up,"NHFIELD_UP",finput,default=0.0d0)
  call parse_input_variable(NHfield_dw,"NHFIELD_DW",finput,default=0.0d0)
  call parse_input_variable(e0,"e0",finput,default=1d0)
  call parse_input_variable(PBC,"PBC",finput,default=.false.)
  call parse_input_variable(lrsym,"LRSYM",finput,default=.true.)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(getbands,"GETBANDS",finput,default=.true.)
  call parse_input_variable(Refine_xmu,"Refine_xmu",finput,default=.false.)
  !
  call ed_read_input(trim(finput),comm)

  !SETUP THE GAMMA MATRICES:
  emat = kron_pauli( pauli_sigma_0, pauli_tau_z)
  !
  soxmat=kron_pauli( pauli_sigma_z, pauli_tau_x)
  soymat=kron_pauli( pauli_sigma_0, pauli_tau_y)
  sozmat=kron_pauli( pauli_sigma_x, pauli_tau_x)
  !
  GapOpeningMat=kron_pauli( pauli_sigma_y, pauli_tau_X)
  NHmat =zero
  NHmat(1:2,1:2) = xi * NHFIELD_UP * pauli_tau_0
  NHmat(3:4,3:4) = xi * NHFIELD_DW * pauli_tau_0
  !

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
  XMU_OLD=XMU
  XMU_REAL=XMU
  converged_xmu=.false.
  if(lrsym)then
     if(mod(Ly,2)/=0)stop "Wrong setup from input file: Ly%2 > 0 (odd number of sites)"
     Nineq=Ly/2
     print*,"Using L-R Symmetry. Solve",Nineq," of",Nlat," sites."
     call sleep(2)
  endif

  !set the local number of total spin-orbitals (4)
  if(Nspin/=2.OR.Norb/=3)stop "Wrong setup from input file: Nspin=2, Norb=3 -> 6Spin-Orbitals"
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
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Weiss_ineq=zero
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Smats_ineq=zero
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal));Sreal_ineq=zero
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Gmats_ineq=zero
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb));Hloc_ineq=zero
  allocate(toconverge(Nineq,6,Lmats));toconverge=zero

  !Buil the Hamiltonian on a grid or on  path
  call build_hkr(trim(hkfile))
  Hloc = lso2nnn(tiHloc,Nlat,Nspin,Norb)
  do ineq=1,Nineq
     ilat = ineq2ilat(ineq)
     Hloc_ineq(ineq,:,:,:,:) = Hloc(ilat,:,:,:,:)
  enddo

  !If postprocessing only
  if(getbands)then
    call read_sigma_matsubara(Smats)
    if(master)then
      do ilat=1,Nlat
        S0(ilat,:,:,:,:) = Smats(ilat,:,:,:,:,1)
      enddo
      call build_eigenbands()
    endif
    call dmft_gloc_matsubara(Hkr,Gmats,Smats,tridiag=tridiag)
    dens=0d0
     do ilat=1,Nlat
       dens = dens+(fft_get_density(Gmats(ilat,1,1,1,1,:),beta)+fft_get_density(Gmats(ilat,2,2,1,1,:),beta))/Nlso
       dens = dens+(fft_get_density(Gmats(ilat,1,1,2,2,:),beta)+fft_get_density(Gmats(ilat,2,2,2,2,:),beta))/Nlso
       dens = dens+(fft_get_density(Gmats(ilat,1,1,3,3,:),beta)+fft_get_density(Gmats(ilat,2,2,3,3,:),beta))/Nlso
     enddo
     print*,"Dens = ",dens
    call finalize_MPI()
    STOP
  endif
  
  !Setup solver
  Nb=ed_get_bath_dimension()
  allocate(Bath_ineq(Nineq,Nb) )
  allocate(Bath_prev(Nineq,Nb) )
  call ed_init_solver(comm,Bath_ineq)

  !DMFT loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")   
     !
     call ed_solve(comm,Bath_ineq,Hloc_ineq)
     ! 
     call ed_get_sigma_matsubara(Smats_ineq,Nineq)
     !
     do ilat=1,Nlat
        ineq = ilat2ineq(ilat)
        Smats(ilat,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:)
        S0(ilat,:,:,:,:)      = Smats_ineq(ineq,:,:,:,:,1)
     enddo
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(Hkr,Gmats,Smats,tridiag=tridiag)
     do ineq=1,Nineq
        ilat = ineq2ilat(ineq)
        Gmats_ineq(ineq,:,:,:,:,:) = Gmats(ilat,:,:,:,:,:)
     enddo
     !
     ! compute the Weiss field (only the Nineq ones)
     call dmft_self_consistency(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,cg_scheme)
     !
     ! fit baths
     call ed_chi2_fitgf(Comm,Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=1)
     call ed_chi2_fitgf(Comm,Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=2)
     !
     call ed_orb_symmetrize_bath(Bath_ineq,2,3)
     !
     !MIXING the current bath with the previous:
     if(iloop>1)Bath_ineq=wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath_ineq
     !
     !Components to converge
     do ilat=1,Nineq
       toconverge(ilat,1,:) = Weiss_ineq(ilat,1,1,2,2,:)
       toconverge(ilat,2,:) = Weiss_ineq(ilat,2,2,2,2,:)
       toconverge(ilat,3,:) = Weiss_ineq(ilat,1,1,3,3,:)
       toconverge(ilat,4,:) = Weiss_ineq(ilat,2,2,3,3,:)
       toconverge(ilat,5,:) = Weiss_ineq(ilat,1,1,1,1,:)
       toconverge(ilat,6,:) = Weiss_ineq(ilat,2,2,1,1,:)
     enddo
     !
     converged = check_convergence(toconverge,dmft_error,nsuccess,nloop)
     !
     !Density evaluation and xmu search
     if(.not. Refine_xmu) then
       dens = 0d0
       do ilat=1,Nlat
         dens = dens+(fft_get_density(Gmats(ilat,1,1,1,1,:),beta)+fft_get_density(Gmats(ilat,2,2,1,1,:),beta))/Nlso
         dens = dens+(fft_get_density(Gmats(ilat,1,1,2,2,:),beta)+fft_get_density(Gmats(ilat,2,2,2,2,:),beta))/Nlso
         dens = dens+(fft_get_density(Gmats(ilat,1,1,3,3,:),beta)+fft_get_density(Gmats(ilat,2,2,3,3,:),beta))/Nlso
       enddo
       if(nread/=0.d0)then
          call ed_search_chemical_potential(xmu,dens,converged)
       endif
       if(converged)refine_xmu=.true.
     else
       call tune_xmu(converged_xmu)
     endif
     !
     converged=converged.and.converged_xmu
     !
     call end_loop
  enddo


  !call read_sigma_real(Sreal)
  !call dmft_gloc_realaxis(Hkr,Greal,Sreal)
  !call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)
  
  do ilat=1,Nlat
     ineq = ilat2ineq(ilat)
     Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
     Smats(ilat,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:)
     S0(ilat,:,:,:,:)      = Smats_ineq(ineq,:,:,:,:,1)
  enddo
  !
  if(master)call build_eigenbands()


  call finalize_MPI()


contains



  !+---------------------------------------------------------------------------+!
  !PURPOSE: build the ti Hamiltonian H(k_x,kz,R_y) on the STRIPE along Y
  !+---------------------------------------------------------------------------+!
  subroutine build_hkr(file)
    character(len=*),optional          :: file
    integer                            :: i,ik
    !
    Lk=Nk**2
    !
    !SETUP THE H(kx,Ry,kz):
    if(master)then
       write(LOGfile,*)"Build H(kx,y,kz) for ti-stripe:"
       write(*,*)"# of kx and kz points     :",Nk
       write(*,*)"# of y-layers      :",Nlat
    endif
    !
    if(allocated(Hkr))deallocate(Hkr)
    allocate(Hkr(Nlso,Nlso,Lk))
    !
    call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])
    call TB_build_model(Hkr,ti_edge_model,Ly,Nso,[Nk,1,Nk],pbc=PBC,wdos=.false.)
    !   
    !
    !SETUP THE LOCAL PART Hloc(Ry)
    if(allocated(tiHloc))deallocate(tiHloc)
    allocate(tiHloc(Nlso,Nlso))
    tiHloc = sum(Hkr(:,:,:),dim=3)/Lk
    where(abs((tiHloc))<1.d-9)tiHloc=0d0
    !call TB_write_Hloc(tiHloc)
    !
  end subroutine build_hkr


  !-----------------------------------------------------------------------------!
  ! purpose: read the local self-energy from disk
  !-----------------------------------------------------------------------------!

  subroutine read_sigma_matsubara(Self)
    complex(8),allocatable,dimension(:,:,:,:,:,:) :: Self,Self_ineq
    character(len=30)                             :: suffix
    integer                                       :: ilat,ispin,iorb,ineq
    real(8),dimension(:),allocatable              :: wm
    call assert_shape(Self,[Nlat,Nspin,Nspin,Norb,Norb,Lmats],"read_sigma_matsubara","Self")
    allocate(wm(Lmats))
    allocate(Self_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    if(master)then
    do ilat=1,Nineq
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix="_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw_ineq"//reg(txtfy(ilat,Npad=4))//".ed"
             call sread("impSigma"//trim(suffix),wm,Self_ineq(ilat,ispin,ispin,iorb,iorb,:))
          enddo
       enddo
    enddo
    do ilat=1,Nlat
     ineq = ilat2ineq(ilat)
     Self(ilat,:,:,:,:,:) = Self_ineq(ineq,:,:,:,:,:)
    enddo
    endif
  end subroutine read_sigma_matsubara
  !
  !
  !
  subroutine read_sigma_real(Self)
    complex(8),allocatable,dimension(:,:,:,:,:,:) :: Self,Self_ineq
    character(len=30)                             :: suffix
    integer                                       :: ilat,ispin,iorb,ineq
    real(8),dimension(:),allocatable              :: wr
    call assert_shape(Self,[Nlat,Nspin,Nspin,Norb,Norb,Lreal],"read_sigma_real","Self_ineq")
    allocate(wr(Lreal))
    allocate(Self_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    wr = linspace(wini,wfin,Lreal)
    if(master)then
    do ilat=1,Nineq
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix="_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw_ineq"//reg(txtfy(ilat,Npad=4))//".ed"
             call sread("impSigma"//trim(suffix),wr,Self_ineq(ilat,ispin,ispin,iorb,iorb,:))
          enddo
       enddo
    enddo
    do ilat=1,Nlat
     ineq = ilat2ineq(ilat)
     Self(ilat,:,:,:,:,:) = Self_ineq(ineq,:,:,:,:,:)
    enddo
    endif
  end subroutine read_sigma_real






  !+---------------------------------------------------------------------------+!
  !PURPOSE: solve H_ti(k_x,R_y,kz) along the 1d -pi:pi path in the BZ.
  !+---------------------------------------------------------------------------+!
  subroutine build_eigenbands(kpath_)
    real(8),dimension(:,:),optional    :: kpath_
    real(8),dimension(:,:),allocatable :: kpath
    type(rgb_color),dimension(:,:),allocatable :: colors
    integer                            :: Npts
    real(8)                            :: offset
    character(len=64)                  :: file
    !
    !PRINT H(kx,Ry) ALONG A -pi:pi PATH
    if(master)write(LOGfile,*)"Solve H(kx,y,kz) along [-Z:Z]:"
    Npts=3
    allocate(Kpath(Npts,3))
    offset=find_x_coordinate(50)/pi
    !
    kpath(1,:)=[offset,0d0,-1d0]*pi
    kpath(2,:)=[offset,0d0,0d0]*pi
    kpath(3,:)=[offset,0d0,1d0]*pi
    file="Eigenbands.nint"
    allocate(colors(Ly,Nso))
    colors = black
    call solve_nh_model(ti_edge_model,Ly,Nso,kpath,Nkpath,&
         colors_name=colors,&
         points_name=[character(len=10) :: "-Z+dx","G+dx","Z+dx","G","X"],&
         file="Eigenbands.nint",&
         pbc=PBC)
  end subroutine build_eigenbands


  function find_x_coordinate(Nintervals) result(xcoord)
    integer                           :: interval,Nintervals
    real(8)                           :: xcoord,xmin,xmax,xcoord_tmp,yval_tmp,yval_tmp_old,yval
    !
    write(LOGfile,*)"Searching for the kx offset:"
    write(LOGfile,*)"Using ",Nintervals," intervals"
    !
    yval_tmp_old=1000
    do interval = 1,Nintervals
      xmin= -pi + 2*pi/Nintervals*interval
      xmax= -pi + 2*pi/Nintervals*(interval+1)
      call  brent_strict(gap_minimum,xcoord_tmp,[xmin,xmax])
      yval_tmp=abs(gap_minimum(xcoord_tmp))
      if(yval_tmp.lt.yval_tmp_old)then 
        xcoord=xcoord_tmp
        yval=yval_tmp
      endif
    enddo
    write(LOGfile,*)"The EP stand at kx = ",xcoord
    write(LOGfile,*)"The real energy gap is ",yval
  end function find_x_coordinate

  function gap_minimum(xcoord) result (gap)
    real(8),dimension(3)                                        :: kpoint
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: Hrk,evec
    real(8)                                                     :: gap,xcoord
    complex(8),dimension(Nlat*Nso)                              :: Eval,eval_sorted
    !
    kpoint=[xcoord,0d0,0d0]
    Hrk=ti_edge_model(kpoint,Nlat,Nso,.false.)
    call eig(hrk,Eval,Evec)
    Eval_sorted=lazy_sort(REAL(Eval))
    gap=abs(Eval_sorted(2*Nlat+1)-Eval_sorted(2*Nlat))
  end function gap_minimum

  function get_e(xcoord) result (gap)
    real(8),dimension(3)                                        :: kpoint
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: Hrk,evec
    real(8)                                                     :: gap,xcoord
    complex(8),dimension(Nlat*Nso)                              :: Eval,eval_sorted
    !
    kpoint=[xcoord,0d0,0d0]
    Hrk=ti_edge_model(kpoint,Nlat,Nso,.false.)
    call eig(hrk,Eval,Evec)
    Eval_sorted=lazy_sort(REAL(Eval))
    gap=(Eval_sorted(2*Nlat+1)+Eval_sorted(2*Nlat))/2d0
  end function get_e

  !+---------------------------------------------------------------------------+!
  !PURPOSE: the ti-edge model hamiltonian
  !+---------------------------------------------------------------------------+!
  !ti on a stripe geometry;
  function ti_edge_model(kpoint,Nlat,N,pbc) result(Hrk)
    real(8),dimension(:)                :: kpoint
    real(8)                             :: kx,kz
    integer                             :: Nlat,N
    complex(8),dimension(N,N)           :: Hmat,Tmat,TmatH
    complex(8),dimension(Nlat*N,Nlat*N) :: Hrk
    integer                             :: i,Idmin,Idmax,Itmin,Itmax
    logical                             :: pbc
    kx=kpoint(1)
    kz=kpoint(3)
    Hrk=zero
    Hmat=h0_rk_ti(kx,kz,N)
    Tmat=t0_rk_ti(N)
    TmatH=conjg(transpose(Tmat))
    do i=1,Nlat
       Idmin=1+(i-1)*N
       Idmax=      i*N
       Hrk(Idmin:Idmax,Idmin:Idmax)=Hmat + select_block(i,S0)
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
  end function ti_edge_model

  function h0_rk_ti(kx,kz,N) result(H)
    real(8)                    :: kx,kz
    integer                    :: N
    complex(8),dimension(4,4)  :: Hblock
    complex(8),dimension(N,N)  :: H
    !
    Hblock = (Mh - e0*(cos(kx) + cos(kz)))*emat+&
         lambda*(sin(kx)*soxmat + sin(kz)*sozmat) + NHmat + GapOpeningField*GapOpeningMat
    !
    H=zero
    !
    H(1,1)=Mh+3*e0
    H(4,4)=Mh+3*e0
    H(2:3,2:3)=Hblock(1:2,1:2)
    H(2:3,5:6)=Hblock(1:2,3:4)
    H(5:6,2:3)=Hblock(3:4,1:2)
    H(5:6,5:6)=Hblock(3:4,3:4)
  end function h0_rk_ti

  function t0_rk_ti(N) result(H)
    integer                    :: N
    complex(8),dimension(4,4)  :: Hblock
    complex(8),dimension(N,N)  :: H
    !
    Hblock = -0.5d0*e0*emat - xi*0.5d0*lambda*soymat
    H=0
    !
    H(2:3,2:3)=Hblock(1:2,1:2)
    H(2:3,5:6)=Hblock(1:2,3:4)
    H(5:6,2:3)=Hblock(3:4,1:2)
    H(5:6,5:6)=Hblock(3:4,3:4)
  end function T0_rk_ti





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

  
  function lso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
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
  end function lso2nnn
  
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


  
  subroutine solve_nh_model(hkr_model,Nlat,Nso,kpath,Nkpath,colors_name,points_name,file,pbc,iproject)
    interface 
       function hkr_model(kpoint,Nlat,Nso,pbc)
         real(8),dimension(:)                    :: kpoint
         integer                                 :: Nlat,Nso
         logical                                 :: pbc
         complex(8),dimension(Nlat*Nso,Nlat*Nso) :: hkr_model
       end function hkr_model
    end interface
    integer                                   :: Nlat,Nso,Nlso
    real(8),dimension(:,:)                    :: kpath
    integer                                   :: Nkpath,Nktot
    type(rgb_color),dimension(Nlat,Nso)       :: colors_name
    character(len=*),dimension(size(kpath,1)) :: points_name
    character(len=*),optional                 :: file
    logical,optional                          :: pbc,iproject
    character(len=256)                        :: file_,file_real,file_imag
    logical                                   :: pbc_,iproject_
    character(len=256)                        :: xtics
    integer                                   :: Npts,Ndim
    integer                                   :: ipts,ik,ic,unit,iorb,ilat,io,nrot,u1,u2
    real(8)                                   :: coeff(Nlat*Nso),klen,ktics(size(Kpath,1))
    type(rgb_color)                           :: corb(Nlat*Nso),c(Nlat*Nso)
    real(8),dimension(size(kpath,2))          :: kstart,kstop,kpoint,kdiff,bk_x,bk_y,bk_z
    complex(8),dimension(Nlat*Nso,Nlat*Nso)   :: h,evec
    complex(8),dimension(Nlat*Nso)            :: Eval
    real(8),allocatable                       :: kseg(:)
    complex(8),allocatable                    :: Ekval(:,:)
    real(8),allocatable                       :: Ekval_sorted(:,:)
    integer,allocatable                       :: Ekcol(:,:)
    !
    !
    file_    = "Eigenbands.tb";if(present(file))file_=file
    file_real=reg(file_)//".real"
    file_imag=reg(file_)//".imag"
    iproject_= .false.        ;if(present(iproject))iproject_=iproject
    pbc_     = .true.         ;if(present(pbc))pbc_=pbc
    !
    Nlso  = Nlat*Nso
    Npts  = size(kpath,1)
    Ndim  = size(kpath,2)
    Nktot = (Npts-1)*Nkpath
    !
    do ilat=1,Nlat
       do io=1,Nso
          corb(io + (ilat-1)*Nso) = colors_name(ilat,io)
       enddo
    enddo
    !
    bk_x=[pi2,0d0,0d0]
    bk_y=[0d0,pi2,0d0]
    bk_z=[0d0,0d0,pi2]
    !
    if(iproject_)then
       select case(Ndim)
       case (1)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x
       case(2)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y
       case (3)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y + kpath(ipts,3)*bk_z
       end select
    endif
    !
    if(master)then
       write(*,*)"Solving model along the path:"
       do ipts=1,Npts
          write(*,"(A,10(A,1x),A1)")"Point"//str(ipts)//": [",(str(kpath(ipts,ic)),ic=1,size(kpath,2)),"]"
       enddo
    endif
    !
    ic=0
    allocate(kseg(Nktot))
    allocate(ekval(Nktot,Nlso))
    allocate(ekval_sorted(Nktot,Nlso))
    allocate(ekcol(Nktot,Nlso))
    klen = 0d0
    if(master)call start_timer()
    do ipts=1,Npts-1
       kstart = kpath(ipts,:)
       kstop  = kpath(ipts+1,:)
       kdiff  = (kstop-kstart)/Nkpath
       ktics(ipts)  = klen
       do ik=1,Nkpath
          ic=ic+1
          kpoint = kstart + (ik-1)*kdiff
          h = hkr_model(kpoint,Nlat,Nso,pbc) - XMU_REAL*Eye(Nlat*Nso)
          call eig(h,Eval,Evec)
          if(master)call eta(ic,Nktot)
          do io=1,Nlso
             coeff(:)=h(:,io)*conjg(h(:,io))
             c(io) = coeff.dot.corb
             Ekval(ic,io) = Eval(io)
             Ekcol(ic,io) = rgb(c(io))
          enddo
          kseg(ic) = klen
          klen = klen + sqrt(dot_product(kdiff,kdiff))
       enddo
    enddo
    ktics(Npts) = Kseg(ic-1)
    if(master)call stop_timer()
    !
    if(master)then
       open(free_unit(unit),file=str(file_real))
       do ic=1,Nktot
        Ekval_sorted(ic,:)=lazy_sort(REAL(Ekval(ic,:)))
       enddo
       do io=1,Nlso
          do ic=1,Nktot
             write(unit,*)kseg(ic),Ekval_sorted(ic,io),Ekcol(ic,io)
          enddo
          write(unit,*)""
       enddo
       close(unit)
       !
       !
       xtics=""
       xtics="'"//reg(points_name(1))//"'"//str(ktics(1))//","
       do ipts=2,Npts-1
          xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//str(ktics(ipts))//","
       enddo
       xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//str(ktics(Npts))//""
       !
       open(unit,file=reg(file_real)//".gp")
       write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
       write(unit,*)"#set out '"//reg(file_real)//".png'"
       write(unit,*)""
       write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
       write(unit,*)"#set out '"//reg(file_real)//".svg'"
       write(unit,*)""
       write(unit,*)"#set term postscript eps enhanced color 'Times'"
       write(unit,*)"#set output '|ps2pdf  -dEPSCrop - "//reg(file_real)//".pdf'"
       write(unit,*)"unset key"
       write(unit,*)"set xtics ("//reg(xtics)//")"
       write(unit,*)"set grid ytics xtics"
       !
       write(unit,*)"plot '"//reg(file_real)//"' every :::0 u 1:2:3 w l lw 3 lc rgb 'black'"
       write(unit,*)"# to print from the i-th to the j-th block use every :::i::j"
       !
       close(unit)
       !
       call system("chmod +x "//reg(file_real)//".gp")
    endif
    !
    if(master)then
       open(free_unit(unit),file=str(file_imag))
       do ic=1,Nktot
        Ekval_sorted(ic,:)=lazy_sort(IMAG(Ekval(ic,:)))
       enddo
       do io=1,Nlso
          do ic=1,Nktot
             write(unit,*)kseg(ic),Ekval_sorted(ic,io),Ekcol(ic,io)
          enddo
          write(unit,*)""
       enddo
       close(unit)
       !
       !
       xtics=""
       xtics="'"//reg(points_name(1))//"'"//str(ktics(1))//","
       do ipts=2,Npts-1
          xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//str(ktics(ipts))//","
       enddo
       xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//str(ktics(Npts))//""
       !
       open(unit,file=reg(file_imag)//".gp")
       write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
       write(unit,*)"#set out '"//reg(file_imag)//".png'"
       write(unit,*)""
       write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
       write(unit,*)"#set out '"//reg(file_imag)//".svg'"
       write(unit,*)""
       write(unit,*)"#set term postscript eps enhanced color 'Times'"
       write(unit,*)"#set output '|ps2pdf  -dEPSCrop - "//reg(file_imag)//".pdf'"
       write(unit,*)"unset key"
       write(unit,*)"set xtics ("//reg(xtics)//")"
       write(unit,*)"set grid ytics xtics"
       !
       write(unit,*)"plot '"//reg(file_imag)//"' every :::0 u 1:2:3 w l lw 3 lc rgb 'black'"
       write(unit,*)"# to print from the i-th to the j-th block use every :::i::j"
       !
       close(unit)
       !
       call system("chmod +x "//reg(file_imag)//".gp")
    endif
    
  end subroutine solve_nh_model
  
  
  subroutine tune_xmu(xmu_converged)
    real(8)          :: ep_energy,ep_location,ep_egap
    logical          :: xmu_converged
    !
    !
    if(master)then
      ep_location=find_x_coordinate(50)
      ep_egap=gap_minimum(ep_location)
      ep_energy=get_e(ep_location)
      !
      if(abs(ep_energy-XMU)<0.01d0)then
        xmu_converged=.true.
      else
        xmu_converged=.false.
        XMU_REAL=ep_energy
        XMU=(1d0-wmixing)*XMU_REAL+wmixing*XMU_OLD
      endif
      XMU_OLD=XMU
    endif
    call bcast_MPI(comm,XMU)
    call bcast_MPI(comm,xmu_converged)
    if(master)then
      print*,"real XMU is ",XMU_REAL
      print*,"mixd XMU is ",XMU
      print*,"Converged = ",xmu_converged
      print*,""
    endif
  end subroutine tune_xmu
  
  
  function lazy_sort(arr_in) result (arr)
    REAL(8), DIMENSION(:)                :: arr_in
    REAL(8), DIMENSION(:), allocatable   :: arr
    INTEGER                              :: i,j,inc,n
    REAL(8)                              :: v
    n=size(arr_in)
    allocate(arr(n))
    arr=arr_in
    inc=1
    do
      inc=3*inc+1
      if (inc > n) exit
    end do
    do
      inc=inc/3
      do i=inc+1,n
        v=arr(i)
        j=i
        do
          if (arr(j-inc) <= v) exit
          arr(j)=arr(j-inc)
          j=j-inc
          if (j <= inc) exit
        end do
        arr(j)=v
      end do
      if (inc <= 1) exit
    end do
  end function lazy_sort
  
    subroutine end_loop_local(unit,id)
    integer,optional :: unit,id
    integer          :: unit_,id_
    logical          :: mpi_master
    unit_=6 ; if(present(unit))unit_=unit
    id_  =0 ; if(present(id))id_=id
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    if(mpi_master)write(unit_,"(A)")"====================================="
    if(mpi_master)write(unit_,"(A)") get_rank_MPI()
    if(mpi_master)write(unit_,*)
    if(mpi_master)write(unit_,*)
  end subroutine end_loop_local
  
  !BRENT MINIMIZING FUNCTION

  subroutine brent_strict(func,xmin,brack,tol,niter)
    interface
       function func(x)
         real(8) :: x
         real(8) :: func
       end function func
    end interface
    real(8),intent(inout)         :: xmin
    real(8),dimension(:),optional :: brack
    real(8),optional              :: tol
    integer,optional              :: niter
    real(8)                       :: tol_
    integer                       :: niter_
    integer                       :: iter
    real(8)                       :: ax,xx,bx,fa,fx,fb,fret
    !
    tol_=1d-6;if(present(tol))tol_=tol
    Niter_=200;if(present(Niter))Niter_=Niter
    !
    if(present(brack))then
       select case(size(brack))
       case(1)
          stop "Brent error: calling brent with size(brack)==1. None or two points are necessary."
       case(2)
          ax = brack(1)
          xx = brack(2)
       case (3)
          ax = brack(1)
          xx = brack(2)
          bx = brack(3)
       end select
    else
       ax=0d0
       xx=1d0
    endif
    fret=brent_optimize(ax,xx,bx,func,tol_,niter_,xmin)
  end subroutine brent_strict
  !



  function brent_optimize(ax,bx,cx,func,tol,itmax,xmin)
    real(8), intent(in)  :: ax,bx,cx,tol
    real(8), intent(out) :: xmin
    real(8)              :: brent_optimize
    integer              :: itmax
    real(8), parameter   :: cgold=0.3819660d0,zeps=1.0d-3*epsilon(ax)
    integer              :: iter
    real(8)              :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    interface
       function func(x)
         real(8) :: x
         real(8) :: func
       end function func
    end interface
    a=min(ax,cx)
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.d0
    fx=func(x)
    fv=fx
    fw=fx
    do iter=1,itmax
       xm=0.5d0*(a+b)
       tol1=tol*abs(x)+zeps
       tol2=2.0*tol1
       if (abs(x-xm) <= (tol2-0.5d0*(b-a))) then
          xmin=x
          brent_optimize=fx
          return
       end if
       if (abs(e) > tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.d0*(q-r)
          if (q > 0.d0) p=-p
          q=abs(q)
          etemp=e
          e=d
          if (abs(p) >= abs(0.5d0*q*etemp) .or. &
               p <= q*(a-x) .or. p >= q*(b-x)) then
             e=merge(a-x,b-x, x >= xm )
             d=cgold*e
          else
             d=p/q
             u=x+d
             if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
          end if
       else
          e=merge(a-x,b-x, x >= xm )
          d=cgold*e
       end if
       u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
       fu=func(u)
       if (fu <= fx) then
          if (u >= x) then
             a=x
          else
             b=x
          end if
          call shft(v,w,x,u)
          call shft(fv,fw,fx,fu)
       else
          if (u < x) then
             a=u
          else
             b=u
          end if
          if (fu <= fw .or. w == x) then
             v=w
             fv=fw
             w=u
             fw=fu
          else if (fu <= fv .or. v == x .or. v == w) then
             v=u
             fv=fu
          end if
       end if
    end do
    !pause 'brent: exceed maximum iterations'

  end function brent_optimize

    subroutine shft(a,b,c,d)
      real(8), intent(out) :: a
      real(8), intent(inout) :: b,c
      real(8), intent(in) :: d
      a=b
      b=c
      c=d
    end subroutine shft

end program ed_ti_slab
