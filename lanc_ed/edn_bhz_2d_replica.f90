program cdn_bhz_2d
  USE DMFT_ED
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE MPI
  !
  implicit none
  integer                                       :: Nx,Ny,Nlso,iloop,Nb,Nkx,Nky,iw,iii,jjj,kkk,Nlat
  integer,dimension(2)                          :: recover
  logical                                       :: converged
  real(8)                                       :: ts,Mh,lambda,wmixing,observable
  !Bath:
  real(8),allocatable                           :: Bath(:),Bath_prev(:)
  !The local hybridization function:
  complex(8),allocatable                        :: Hloc(:,:)
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Gmats,Greal,Smats,Sreal,Weiss,Weiss_old
  character(len=16)                             :: finput
  real(8),allocatable                           :: wt(:)
  complex(8),allocatable                        :: wm(:),wr(:)
  complex(8),allocatable                        :: Hk(:,:,:),Smats_lso(:,:,:)
  !SYMMETRIES TEST
  real(8),dimension(:),allocatable              :: lambdasym_vector
  complex(8),dimension(:,:,:,:,:),allocatable   :: Hsym_basis
  !MPI VARIABLES (local use -> ED code has its own set of MPI variables)
  integer                                       :: comm
  integer                                       :: rank
  integer                                       :: mpi_size
  logical                                       :: master

  !Init MPI: use of MPI overloaded functions in SciFor
  call init_MPI(comm,.true.)
  rank   = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  !

  !Parse input variables
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
  call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="hopping parameter")
  call parse_input_variable(Mh,"Mh",finput,default=1.d0,comment="crystal field splitting")
  call parse_input_variable(lambda,"lambda",finput,default=0.3d0,comment="spin-orbit coupling")
  call parse_input_variable(Nx,"Nx",finput,default=1,comment="Number of cluster sites in x direction")
  call parse_input_variable(Ny,"Ny",finput,default=1,comment="Number of cluster sites in y direction")
  call parse_input_variable(Nkx,"Nkx",finput,default=10,comment="Number of kx point for BZ integration")
  call parse_input_variable(Nky,"Nky",finput,default=10,comment="Number of ku point for BZ integration")

  !
  call ed_read_input(trim(finput),comm)
  !
  !Add dmft control variables
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  !set global variables
  !if (Nspin/=1.or.Norb/=1) stop "You are using too many spin-orbitals"
  Nlat=Nx*Ny
  Nlso=Nlat*Nspin*Norb
  if(.not.allocated(wm))allocate(wm(Lmats))
  wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)

  if(ED_VERBOSE > 0)call naming_convention()
  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats),Weiss_old(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats),Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats_lso(Nlso,Nlso,Lmats))

  !Build Hk and Hloc
  call generate_hk_hloc()
  allocate(lambdasym_vector(3))
  allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))

  !SETUP SYMMETRIES (EXPERIMENTAL)
  lambdasym_vector(1)=Mh
  Hsym_basis(:,:,:,:,1)=so2nn(hloc_model(Nlso,1.d0,0.d0,0.d0))
  !
  lambdasym_vector(2)=ts
  Hsym_basis(:,:,:,:,2)=so2nn(hloc_model(Nlso,0.d0,1.d0,0.d0))
  !
  lambdasym_vector(3)=lambda
  Hsym_basis(:,:,:,:,3)=so2nn(hloc_model(Nlso,0.d0,0.d0,1.d0))
  !
  !setup solver
  call set_Hloc(Hsym_basis,lambdasym_vector)
  Nb=get_bath_dimension(dreal(Hsym_basis));print*,Nb
  allocate(bath(Nb))
  allocate(bath_prev(Nb))
  call ed_init_solver(comm,bath)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath) 
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)

     !Compute the local gfs:
     call dmft_gloc_matsubara(comm,Hk,Wt,Gmats,Smats)
     if(master)call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
     !
     !Get the Weiss field/Delta function to be fitted
     call dmft_self_consistency(comm,Gmats,Smats,Weiss,so2nn(Hloc),cg_scheme)
     if(master)call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=4)
     ! call Bcast_MPI(comm,Weiss)
     !
     !
     !Perform the SELF-CONSISTENCY by fitting the new bath
     if(master)then
        call ed_chi2_fitgf(Weiss,bath)
        !
        !MIXING:
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_Prev
        Bath_Prev=Bath
        !
        !Check convergence (if required change chemical potential)
        converged = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop)
     endif
     !
     call Bcast_MPI(comm,bath)
     call Bcast_MPI(comm,converged)
     call Bcast_MPI(comm,xmu)
     !
     if(master)call end_loop
  enddo

  !Compute the local gfs:
  call dmft_gloc_realaxis(comm,Hk,Wt,Greal,Sreal)
  if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)

  !Compute the Kinetic Energy:
  do iw=1,Lmats
     Smats_lso(:,:,iw)=nn2so(Smats(:,:,:,:,iw))
  enddo
  call dmft_kinetic_energy(comm,Hk(:,:,:),Wt,Smats_lso)


  call finalize_MPI()


contains

  !+------------------------------------------------------------------+
  !PURPOSE  : Hloc for the 2d BHZ model
  !+------------------------------------------------------------------+


  function Hloc_model(N,Mh_,ts_,lambda_) result (H0)
    integer                                               :: N,ilat,jlat,ispin,iorb,jorb,ind1,ind2
    real(8)                                               :: Mh_,ts_,lambda_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hopping_matrix
    complex(8),dimension(N,N)                             :: H0
    character(len=64)                                     :: file_
    integer                                               :: unit
    file_ = "tcluster_matrix.dat"
    !
    hopping_matrix=zero
    !
    do ispin=1,Nspin
       do ilat=1,Nx
          do jlat=1,Ny
             ind1=indices2N([ilat,jlat])
             hopping_matrix(ind1,ind1,ispin,ispin,:,:)= t_m(mh_)
             if(ilat<Nx)then
                ind2=indices2N([ilat+1,jlat])
                hopping_matrix(ind1,ind2,ispin,ispin,:,:)= t_x(ts_,lambda_,ispin)
             endif
             if(ilat>1)then
                ind2=indices2N([ilat-1,jlat])
                hopping_matrix(ind1,ind2,ispin,ispin,:,:)= dconjg(transpose(t_x(ts_,lambda_,ispin)))
             endif
             if(jlat<Ny)then
                ind2=indices2N([ilat,jlat+1])
                hopping_matrix(ind1,ind2,ispin,ispin,:,:)= t_y(ts_,lambda_)
             endif
             if(jlat>1)then
                ind2=indices2N([ilat,jlat-1])
                hopping_matrix(ind1,ind2,ispin,ispin,:,:)= transpose(t_y(ts_,lambda_))
             endif
          enddo
       enddo
    enddo
    !
    H0=nn2so(hopping_matrix(1,1,:,:,:,:))
    !
  end function hloc_model


  function hk_model(kpoint,N) result(Hk)
    integer                                                      :: N,ilat,jlat,ispin,iorb,jorb,i,j,ind1,ind2
    real(8),dimension(:)                                         :: kpoint
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: hopping_matrix
    complex(8),dimension(N,N)                                    :: hk
    !
    hopping_matrix=zero
    !
    do ispin=1,Nspin
       do ilat=1,Ny
          ind1=indices2N([1,ilat])
          ind2=indices2N([Nx,ilat])
          hopping_matrix(ind1,ind2,ispin,ispin,:,:)=hopping_matrix(ind1,ind2,ispin,ispin,:,:) + dconjg(transpose(t_x(ts,lambda,ispin)))*exp(xi*kpoint(1)*Nx)
          hopping_matrix(ind2,ind1,ispin,ispin,:,:)=hopping_matrix(ind2,ind1,ispin,ispin,:,:) + t_x(ts,lambda,ispin)*exp(-xi*kpoint(1)*Nx)
       enddo
       do ilat =1,Nx
          ind1=indices2N([ilat,1])
          ind2=indices2N([ilat,Ny])
          hopping_matrix(ind1,ind2,ispin,ispin,:,:)=hopping_matrix(ind1,ind2,ispin,ispin,:,:) + transpose(t_y(ts,lambda))*exp(xi*kpoint(2)*Ny)
          hopping_matrix(ind2,ind1,ispin,ispin,:,:)=hopping_matrix(ind2,ind1,ispin,ispin,:,:) + t_y(ts,lambda)*exp(-xi*kpoint(2)*Ny)
       enddo
    enddo
    !
    Hk=nn2so(hopping_matrix(1,1,:,:,:,:))+hloc_model(N,Mh,ts,lambda)
    ! 
  end function hk_model


  !AUXILLIARY HOPPING MATRIX CONSTRUCTORS

  function t_m(mass) result(tmpmat)
    complex(8),dimension(Norb,Norb) :: tmpmat
    real(8)                         :: mass
    !
    tmpmat=zero
    tmpmat=mass*pauli_sigma_z
    !
  end function t_m

  function t_x(hop1,hop2,spinsign) result(tmpmat)
    complex(8),dimension(Norb,Norb) :: tmpmat
    real(8)                         :: hop1,hop2,sz
    integer                         :: spinsign
    !
    tmpmat=zero
    sz=(-1.d0)**(spinsign+1)
    tmpmat=-hop1*pauli_sigma_z+0.5d0*sz*xi*hop2*pauli_sigma_x
    !
  end function t_x

  function t_y(hop1,hop2) result(tmpmat)
    complex(8),dimension(Norb,Norb) :: tmpmat
    real(8)                         :: hop1,hop2
    !
    tmpmat=zero
    tmpmat=-hop1*pauli_sigma_z
    tmpmat(1,2)=-hop2*0.5d0
    tmpmat(2,1)=hop2*0.5d0
    !
  end function t_y

  !-------------------------------------------------------------------------------------------
  !PURPOSE: generate Hloc and Hk
  !-------------------------------------------------------------------------------------------

  subroutine generate_hk_hloc()
    integer                                     :: ik
    real(8),dimension(Nkx*Nky,2)                :: kgrid
    !
    call TB_build_kgrid([Nkx,Nky],kgrid)
    kgrid(:,1)=kgrid(:,1)/Nx
    kgrid(:,2)=kgrid(:,2)/Ny
    !
    if(allocated(hk))deallocate(hk)
    if(allocated(Wt))deallocate(Wt)
    if(allocated(hloc))deallocate(Hloc)
    !
    allocate(Hk(Nlso,Nlso,Nkx*Nky),Wt(Nkx*Nky),Hloc(Nlso,Nlso))
    hk=zero
    wt=zero
    hloc=zero
    !
    call TB_build_model(Hk,hk_model,Nlso,kgrid)
    Wt = 1d0/(Nkx*Nky)
    Hloc=hloc_model(Nlso,Mh,ts,lambda)
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
    !
  end subroutine generate_hk_hloc

  !+------------------------------------------------------------------+
  !PURPOSE  : Auxilliary reshape functions
  !+------------------------------------------------------------------+

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


  function indices2N(indices) result(N)
    integer,dimension(2)         :: indices
    integer                      :: N,i
    !
    !
    N=Nx*(indices(2)-1)+indices(1)
  end function indices2N

  function N2indices(N) result(indices) 
    integer,dimension(2)         :: indices
    integer                      :: N,i,N_
    !
    indices(1)=mod(N,Nx)
    if(indices(1)==0)indices(1)=3
    indices(2)=N/Nx
  end function N2indices

  subroutine naming_convention()
    integer                       :: i,j
    integer,dimension(Nx,Ny)      :: matrix
    !
    do j=1,Ny
       do i=1,Nx
          matrix(i,j)=indices2N([i,j])
       enddo
    enddo
    !
    write(LOGfile,"(A)")"The unique index of each site (on the cartesian plane) is as follows:"
    write(LOGfile,"(A)")" "
    do j=1,Ny
       write(LOGfile,"(20(I2,2x))")(matrix(i,Ny+1-j),i =1,Nx)
    enddo
    write(LOGfile,"(A)")" "
  end subroutine naming_convention
  !
end program cdn_bhz_2d








