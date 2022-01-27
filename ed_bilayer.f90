program ed_bilayer
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                             :: iloop,Lk,Nso
  logical                             :: converged
  !Bath:
  integer                             :: Nb
  real(8),allocatable                 :: Bath(:),Bath_prev(:)
  !The local hybridization function:
  complex(8),allocatable              :: Weiss(:,:,:,:,:)
  complex(8),allocatable              :: Smats(:,:,:,:,:)
  complex(8),allocatable              :: Sreal(:,:,:,:,:)
  complex(8),allocatable              :: Gmats(:,:,:,:,:)
  complex(8),allocatable              :: Greal(:,:,:,:,:)
  complex(8),allocatable,dimension(:) :: Gtest
  real(8),allocatable,dimension(:)    :: dens
  !hamiltonian input:
  !variables for the model:
  complex(8),allocatable              :: Hk(:,:,:)
  complex(8),allocatable              :: Hloc(:,:)
  complex(8),allocatable              :: SigmaHk(:,:),Zmats(:,:)
  integer                             :: Nk,Nkpath
  real(8)                             :: ts,tperp,Vel,lambda,wmixing
  character(len=16)                   :: finput
  !MPI Vars:
  integer                             :: comm,rank
  logical                             :: master


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BL.in')  
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(ts,"TS",finput,default=1d0)
  call parse_input_variable(tperp,"TPERP",finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
  call parse_input_variable(Vel,"Vel",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)

  !
  call ed_read_input(trim(finput),comm)
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(Nspin/=1.OR.Norb/=2)stop "Wrong setup from input file: Nspin=1, Norb=2 -> 2Spin-Orbitals"
  Nso=Nspin*Norb

  !Allocate Weiss Field:

  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gtest(Lmats))
  allocate(dens(Norb))
  allocate(SigmaHk(Nso,Nso))
  allocate(Zmats(Nso,Nso))

  !Buil the Hamiltonian on a grid or on  path
  call set_SigmaHk()
  call build_hk()

  
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_prev(Nb))
  call ed_init_solver(comm,bath)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath,Hloc=j2so(Hloc))
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)
     call ed_get_dens(dens)


     !Get GLOC:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)


     !Update WeissField:
     call dmft_self_consistency(Gmats,Smats,Weiss,j2so(Hloc),cg_scheme)
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)


     !Fit the new bath, starting from the old bath + the supplied Weiss
     call ed_chi2_fitgf(comm,Weiss,bath,ispin=1)
     call ed_spin_symmetrize_bath(bath,save=.true.)

     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath
     !

     Gtest=Weiss(1,1,1,1,:)
     converged = check_convergence(Gtest,dmft_error,nsuccess,nloop)
     if(nread/=0d0)call ed_search_variable(xmu,sum(dens),converged)

     call end_loop
  enddo


  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)

  call dmft_kinetic_energy(Hk,Smats)
  call solve_hk_bands(so2j(Smats(:,:,:,:,1)))

  call save_array("Smats",Smats)
  call save_array("Sreal",Sreal)

  call finalize_MPI()



contains


  !---------------------------------------------------------------------
  !PURPOSE: GET HAMILTONIAN (from the NonInteracting code)
  !---------------------------------------------------------------------
  subroutine build_hk()
    !
    if(allocated(Hk))deallocate(Hk)
    if(allocated(Hloc))deallocate(Hloc)
    !
    call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
    !
    if(master)write(LOGfile,*)"Build H(k):"
    Lk=Nk*Nk
    if(master)write(*,*)"# of k-points     :",Lk
    if(master)write(*,*)"# of SO-bands     :",Nso
    !
    allocate(Hk(Nso,Nso,Lk)) ;Hk=zero
    call TB_build_model(Hk,hk_model,Nso,[Nk,Nk])
    !
    allocate(Hloc(Nso,Nso))
    Hloc = zero
    Hloc = sum(Hk,dim=3)/Lk
    where(abs(Hloc)<1d-6)Hloc=zero
    if(master) call TB_write_Hloc(Hloc)
  end subroutine build_hk






  !--------------------------------------------------------------------!
  !PURPOSE: Set the Self-Energy
  !--------------------------------------------------------------------!
  subroutine set_SigmaHk(sigma)
    complex(8),dimension(Nso,Nso),optional :: sigma(Nso,Nso)
    SigmaHk = zero;if(present(sigma))SigmaHk=sigma
  end subroutine set_SigmaHk


  !--------------------------------------------------------------------!
  !PURPOSE: Solve the topological Hamiltonian
  !--------------------------------------------------------------------!
  subroutine solve_hk_bands(sigma)
    integer                                :: Npts
    complex(8),dimension(Nso,Nso)          :: sigma(Nso,Nso)
    real(8),dimension(:,:),allocatable     :: kpath
    !
    if(master)then
       write(LOGfile,*)"Build G^-1(k,0) along path:"
       !
       call set_SigmaHk()
       !
       Npts = 4
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_gamma
       kpath(2,:)=kpoint_x1
       kpath(3,:)=kpoint_m1
       kpath(4,:)=kpoint_gamma
       call set_SigmaHk(sigma)
       call TB_solve_model(hk_model_renorm,Nso,kpath,Nkpath,&
            colors_name=[red,blue],&
            points_name=[character(len=20) :: "{\Symbol G}","X","M","{\Symbol G}"],&
            file="Eig_Htop.ed")
    endif
  end subroutine solve_hk_bands




  !--------------------------------------------------------------------!
  ! HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_model(kvec,N) result(hk)
    integer                   :: N
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: kx,ky
    real(8)                   :: ek,vk
    if(N/=Nso)stop "hk_model error: N != Nspin*Norb == 2"
    kx=kvec(1)
    ky=kvec(2)
    ek = -2*ts*(cos(kx)+cos(ky))
    vk = tperp + lambda*(cos(kx)-cos(ky))   
    Hk = Vel*pauli_tau_z + ek*pauli_tau_0 + vk*pauli_tau_x
  end function hk_model


  function hk_model_renorm(kvec,N) result(hk)
    integer                   :: N
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: kx,ky
    real(8)                   :: ek,vk
    integer                   :: ii
    if(N/=Nso)stop "hk_model error: N != Nspin*Norb == 2"
    kx=kvec(1)
    ky=kvec(2)
    ek = -2*ts*(cos(kx)+cos(ky))
    vk = tperp + lambda*(cos(kx)-cos(ky))   
    Hk = Vel*pauli_tau_z + ek*pauli_tau_0 + vk*pauli_tau_x
    !
    Hk = Hk + dreal(SigmaHk)
    Zmats=zero
    do ii=1,Nso
       Zmats(ii,ii)  = 1d0/abs( 1d0 +  abs(dimag(SigmaHk(ii,ii))/(pi/beta)) )
    end do
    Hk = matmul(Zmats,Hk)
  end function hk_model_renorm







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


end program
