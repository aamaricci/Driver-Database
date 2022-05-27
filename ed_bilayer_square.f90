program ed_bilayer
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                     :: iloop,Lk,Nso,Nlso,Nlat
  logical                                     :: converged,conv_dens
  integer :: ilat
  !Bath:
  integer                                     :: Nb
  real(8),allocatable                         :: Bath(:,:),Bath_prev(:,:)
  !The local hybridization function:
  complex(8),allocatable                      :: Weiss(:,:,:,:,:,:)
  complex(8),allocatable                      :: Smats(:,:,:,:,:,:)
  complex(8),allocatable                      :: Sreal(:,:,:,:,:,:)
  complex(8),allocatable                      :: Gmats(:,:,:,:,:,:)
  complex(8),allocatable                      :: Greal(:,:,:,:,:,:)
  complex(8),allocatable,dimension(:,:)         :: Gtest
  real(8),allocatable,dimension(:)            :: dens,dens_prev,docc
  !hamiltonian input:
  !variables for the model:
  complex(8),allocatable                      :: Hk(:,:,:),Hk_loop(:,:,:)
  complex(8),allocatable                      :: Hloc(:,:,:,:,:)
  complex(8),allocatable                      :: SigmaHk(:,:),Zmats(:,:)
  integer                                     :: Nk,Nkpath
  real(8)                                     :: ts,tperp,Vel,lambda,wmixing,ntarget,wmixing_dens
  character(len=16)                           :: finput
  !MPI Vars:
  integer                                     :: comm,rank,mpierr
  logical                                     :: master
  !
  !
  real(8),dimension(2)                        :: e1,e2   !real-space lattice basis
  real(8),dimension(2)                        :: bk1,bk2 !reciprocal space lattice basis
  real(8),dimension(2)                        :: bklen
  !
  real(8)                                     :: Vnn  !+- nearest neighbourg interaction
  real(8)                                     :: ntop,nbot
  integer                                     :: uio,uio_loop,imu
  real(8)                                     :: ntest1,ntest2,xmu_imp,xmu1,xmu2,alpha_mu,dens_error
  logical                                     :: fix_mu,flag_mpi
  logical                                     :: fixing_newton

  character(len=6)     :: fix_mu_scheme           !newton/f_zero

  !
  !
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

  call parse_input_variable(Vnn,"Vnn",finput,default=0.d0)
  call parse_input_variable(ntop,"ntop",finput,default=1.d0)
  call parse_input_variable(nbot,"nbot",finput,default=1.d0)
  
  call parse_input_variable(alpha_mu,"alpha_mu",finput,default=0.5d0)

  call parse_input_variable(fix_mu,"fix_mu",finput,default=.false.)

  call parse_input_variable(flag_mpi,"flag_mpi",finput,default=.false.)
  call parse_input_variable(conv_dens,"conv_dens",finput,default=.false.)
  
  call parse_input_variable(dens_error,"dens_error",finput,default=1.d-5)

  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)

  call parse_input_variable(wmixing_dens,"WMIXING_DENS",finput,default=0.5d0)

  call parse_input_variable(ntarget,"NTARGET",finput,default=2d0)

  call parse_input_variable(fix_mu_scheme,"fix_mu_scheme",finput,default='f_zero')

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


  Nlat=2
  if(Nspin/=1.OR.Nlat/=2)stop "Wrong setup from input file: Nspin=1, Norb=2 -> 2Spin-Orbitals"
  Nso=Nspin*Norb
  Nlso=Nlat*Nso
  
  !Allocate Weiss Field:

  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gtest(Nlat,Lmats))
  allocate(dens(Norb)); allocate(dens_prev(Norb))
  allocate(docc(Nlat)); 
  allocate(SigmaHk(Nso,Nso))
  allocate(Zmats(Nso,Nso))

  !Buil the Hamiltonian on a grid or on  path
  call set_SigmaHk()
  
  !+- a=1 is the hexagon edge
  !+- a' = sqrt(3.0)*a is the modulus of the primitive vector of the triangular lattice
  !+- e_1 = sqrt(3.0)/2.0 * a [-1, sqrt(3.0)]
  !+- e_2 = sqrt(3.0)/2.0 * a [ 1, sqrt(3.0)]
  !e1 = sqrt(3d0)/2d0*[-1d0, sqrt(3d0)]
  !e2 = sqrt(3d0)/2d0*[ 1d0, sqrt(3d0)]
  e1 = [1d0,0d0]
  e2 = [0d0,1d0]
  
  !RECIPROCAL LATTICE VECTORS:
  bk1=2*pi*[ 1d0, 0d0]
  bk2=2*pi*[ 0d0, 1d0]
  call TB_set_bk(bkx=bk1,bky=bk2) 
  call build_hk_square('Hk_square.out')
  !
  !
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nlat,Nb))
  allocate(Bath_prev(Nlat,Nb))
  call ed_init_solver(comm,bath)
  Bath_prev=Bath
  
  !
  !DMFT loop
  iloop=0;converged=.false.
  dens=[ntop,nbot]
  dens_prev=dens
  uio=free_unit()
  if(master) then
     open(uio,file='init_ndens_hf.out')
     write(uio,'(3F18.10)') dens
     close(uio)
  end if
  !
  uio=free_unit()
  if(master) then
     open(uio,file='iloop_observables.out')
     close(uio)
  end if
  call set_Hloc(dens)    
  !
  !
  xmu1=0.d0
  xmu2=xmu1


  xmu_imp=xmu
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call set_Hloc(dens,'Hloc_set'//str(iloop,Npad=3)//'.dat')     
     
     !+- set fix_mu=F and skip all this part
     if(fix_mu) then
        !   bracketing the xmu  !

        select case(fix_mu_scheme)
        case('newton')
           call newton(get_delta_dens_imp_,xmu,eps=1d-4)
        case('f_zero')
           xmu1= -50
           xmu2=  50
           call fzero(get_delta_dens_imp,xmu1,xmu2,imu,tol_rel=1.d-8,tol_abs=1.d-6)
           xmu=xmu1
        case default           
           xmu1= -50
           xmu2=  50
           call fzero(get_delta_dens_imp,xmu1,xmu2,imu,tol_rel=1.d-8,tol_abs=1.d-6)
           xmu=xmu1
        end select
        !
     end if
     !
     call ed_solve(comm,bath,Hloc,mpi_lanc=flag_mpi)     
     call ed_get_sigma_matsubara(Smats,Nlat)
     call ed_get_sigma_realaxis(Sreal,Nlat)

     !
     !
     call ed_get_dens(dens,Nlat,iorb=1)
     call ed_get_docc(docc,Nlat,iorb=1)
     !
     uio=free_unit()
     if(master) then
        open(uio,file='iloop_observables.out',status='old',position='append')
        write(uio,'(10F18.10)') dens,docc,xmu
        close(uio)
     end if
     
     ! !Get GLOC:
     
     !+- here I should use Hk_loop !!!!
     if(.not.allocated(hk_loop)) then
        call mpi_barrier(comm,mpiERR)
        stop
     end if
     call dmft_gloc_matsubara(Hk_loop,Gmats,Smats)
     call dmft_gloc_realaxis(Hk_loop,Greal,Sreal)
     
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4,ineq_pad=2)
     call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4,ineq_pad=2)
     
     ! !Update WeissField:
     call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,cg_scheme)
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1,ineq_pad=2)
               
     !
     call ed_chi2_fitgf(bath,Weiss,Hloc,ispin=1)  !+- perche qui si mangia anche Hloc?
     !
     !+- mix things
     Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     dens = wmixing_dens*dens + (1.d0-wmixing_dens)*dens_prev
     
     Bath_prev = Bath
     dens_prev = dens
     xmu_imp   = xmu
     !
     Gtest=Weiss(:,1,1,1,1,:)
     if(.not.conv_dens) then
        converged = check_convergence(Gtest,dmft_error,nsuccess,nloop)
     else
        converged = check_convergence(dens,dens_error,nsuccess,nloop)
     end if


     !if(nread/=0d0)call ed_search_variable(xmu,sum(dens),converged)
     
     call end_loop


  enddo


  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)

  call dmft_kinetic_energy(Hk,Smats)

  call save_array("Smats",Smats)
  call save_array("Sreal",Sreal)

  call finalize_MPI()

contains
  
  function get_delta_dens_imp(xmu_in) result(dens_out)
    real(8),intent(in) :: xmu_in
    real(8) :: dens_out
    character(len=20) :: suffix
    !
    xmu=xmu_in
    !
    call ed_solve(comm,bath,Hloc,mpi_lanc=flag_mpi,iflag=.false.)
    call ed_get_dens(dens,Nlat,iorb=1)
    !
    call ed_reset_suffix    
    dens_out=dens(1)+dens(2)-ntarget        
    !
  end function get_delta_dens_imp


  function get_delta_dens_imp_(xmu_in) result(dens_out)
    real(8) :: xmu_in
    real(8) :: dens_out
    character(len=20) :: suffix
    !
    xmu=xmu_in
    !
    call ed_solve(comm,bath,Hloc,iflag=.false.)
    call ed_get_dens(dens,Nlat,iorb=1)
    !
    call ed_reset_suffix    
    dens_out=(dens(1)+dens(2)-ntarget)**2.d0        
    !
  end function get_delta_dens_imp_
    
  

  subroutine build_hk_honeycomb(file)
    character(len=*),optional          :: file
    real(8),dimension(2)               :: pointK,pointKp,pointM
    real(8),dimension(:,:),allocatable :: KPath
    real(8),dimension(:,:),allocatable :: kgrid
    real(8),dimension(:),allocatable   :: gridx,gridy
    integer                            :: i,j,ik
    !
    Lk= Nk*Nk
    if(master)write(*,*)"Build H(k) for the honeycomb lattice:",Lk
    if(master)write(*,*)"# of SO-bands     :",Nlso
    !
    if(allocated(Hk))deallocate(Hk)
    !
    allocate(Hk(Nlso,Nlso,Lk));Hk=zero
    !
    call TB_build_model(Hk,hk_honeycomb,Nlso,[Nk,Nk])
    if(present(file).and.master) call TB_write_Hk(Hk,file,1,Nspin,Norb,[Nk,Nk])
    !
    pointK  = 1d0/3d0*bk1 + 2d0/3d0*bk2
    pointKp = 2d0/3d0*bk1 + 1d0/3d0*bk2
    pointM= 0.5d0*pointK+0.5d0*pointKp
    write(*,*) pointK
    write(*,*) pointKp
    allocate(Kpath(4,2))
    KPath(1,:)=[0,0]
    KPath(2,:)=pointK
    Kpath(3,:)=pointM
    KPath(4,:)=[0d0,0d0]
    !
    call TB_Solve_model(hk_honeycomb,Nlso,KPath,Nkpath,&
         colors_name=[red1,blue1,red1,blue1],&
         points_name=[character(len=10) :: "G","K","K`","G"],&
         file="Eigenbands.nint",iproject=.false.)
    !
  end subroutine build_hk_honeycomb

  subroutine build_hk_square(file)
    character(len=*),optional          :: file
    real(8),dimension(2)               :: pointK,pointKp,pointM
    real(8),dimension(:,:),allocatable :: KPath
    real(8),dimension(:,:),allocatable :: kgrid
    real(8),dimension(:),allocatable   :: gridx,gridy
    integer                            :: i,j,ik
    !
    Lk= Nk*Nk
    if(master)write(*,*)"Build H(k) for the honeycomb lattice:",Lk
    if(master)write(*,*)"# of SO-bands     :",Nlso
    !
    if(allocated(Hk))deallocate(Hk)
    !
    allocate(Hk(Nlso,Nlso,Lk));Hk=zero
    !
    call TB_build_model(Hk,hk_square,Nlso,[Nk,Nk])
    if(present(file).and.master) call TB_write_Hk(Hk,file,1,Nspin,Norb,[Nk,Nk])
    !
    !+- nomi di fantasia -+!
    pointK  = 0.d0
    pointKp = bk1*0.5d0 
    pointM= 0.5d0*(bk1+bk2)
    allocate(Kpath(4,2))
    KPath(1,:)=[0d0,0d0]
    KPath(2,:)=pointKp
    Kpath(3,:)=pointM
    KPath(4,:)=[0d0,0d0]
    !
    call TB_Solve_model(hk_square,Nlso,KPath,Nkpath,&
         colors_name=[red1,blue1,red1,blue1],&
         points_name=[character(len=10) :: "G","K","K`","G"],&
         file="Eigenbands.nint",iproject=.false.)
    !
  end subroutine build_hk_square

  subroutine set_hloc(nave,file_out)
    real(8),dimension(2) :: nave
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hloc_tmp
    character(len=*),optional :: file_out
    integer :: unit,iso,ik
    !
    Hloc_tmp = zero
    !+- sum the k Hk -+!
    Hloc_tmp = sum(Hk,dim=3)/Lk    
    where(abs(Hloc_tmp)<1d-6) Hloc_tmp=zero
    !
    if(allocated(Hloc)) deallocate(Hloc)    
    allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))    
    !
    if(allocated(Hk_loop)) deallocate(Hk_loop)
    allocate(Hk_loop(Nlso,Nlso,Lk)) ; Hk_loop=Hk
    !
    iso=lso2j_index(1,1,1)
    Hloc_tmp(iso,iso) = Hloc_tmp(iso,iso) + 3.d0*Vnn*nave(2) !- 3.d0*Vnn*(nave(2)+nave(1))*0.5d0
    !
    iso=lso2j_index(2,1,1)
    Hloc_tmp(iso,iso) = Hloc_tmp(iso,iso) + 3.d0*Vnn*nave(1) !- 3.d0*Vnn*(nave(2)+nave(1))*0.5d0
    !
    forall(ik=1:Lk) Hk_loop(:,:,ik) = Hk_loop(:,:,ik) + Hloc_tmp
    !
    Hloc = lso2nnn_reshape(Hloc_tmp,Nlat,Nspin,Norb)
    !
    if(present(file_out).and.master) then
       unit=free_unit()
       open(unit,file=trim(file_out))
       do ilat=1,Nlat
          write(unit,'(10F18.10)') Hloc(ilat,:,:,:,:)
       end do
       close(unit)
    end if
    !
  end subroutine set_hloc

  !--------------------------------------------------------------------!
  !PURPOSE: Set the Self-Energy
  !--------------------------------------------------------------------!
  subroutine set_SigmaHk(sigma)
    complex(8),dimension(Nso,Nso),optional :: sigma(Nso,Nso)
    SigmaHk = zero;if(present(sigma))SigmaHk=sigma
  end subroutine set_SigmaHk
  !
  function hk_honeycomb(kpoint,Nlso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nlso
    complex(8),dimension(2,2)       :: hk11,hk22
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                         :: h0,hx,hy,hz, kdote1, kdote2
    !
    kdote1 = dot_product(kpoint,e1)
    kdote2 = dot_product(kpoint,e2)
    !
    h0 = -2*ts*( cos(kdote1) + cos(kdote2) + cos(kdote1+kdote2) )
    hx = -tperp*( cos(kdote1) + cos(kdote2) + 1)
    hy = -tperp*( sin(kdote1) + sin(kdote2) )
    hz = Vel
    !
    hk = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z 
    !
  end function hk_honeycomb


  
  function hk_square(kpoint,Nlso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nlso
    complex(8),dimension(2,2)       :: hk11,hk22
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                         :: h0,hx,hy,hz, kdote1, kdote2
    !
    kdote1 = dot_product(kpoint,e1)
    kdote2 = dot_product(kpoint,e2)
    !
    h0 = -2*ts*( cos(kdote1) + cos(kdote2) )
    hx = -tperp*( cos(kdote1) + cos(kdote2) + 1 + cos(kdote1+kdote2))
    hy = -tperp*( sin(kdote1) + sin(kdote2) + sin(kdote1+kdote2) )
    hz = Vel
    !
    hk = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z 
    !
  end function hk_square
  



  function lso2nnn_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Fin
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Fout(ilat,ispin,jspin,iorb,jorb) = Fin(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function lso2nnn_reshape



  function lso2j_index(ilat,ispin,iorb) result(isporb)
    integer :: ispin,iorb,ilat
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop "error so2j_index: ispin>Nspin"
    isporb=(ilat-1)*Norb*Nspin + (ispin-1)*Nspin + iorb
  end function lso2j_index


  function lso2j(fg,Nlat) result(g)
    integer :: Nlat
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin,ilat
    g=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   i=lso2j_index(ilat,ispin,iorb)
                   j=lso2j_index(ilat,jspin,jorb)
                   g(ilat,i,j) = fg(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    end do
  end function lso2j

  function j2lso(fg,Nlat) result(g)
    integer :: Nlat
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb) :: fg
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: g
    integer                                          :: i,j,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   i=lso2j_index(ilat,ispin,iorb)
                   j=lso2j_index(ilat,jspin,jorb)
                   g(ilat,ispin,jspin,iorb,jorb)  = fg(ilat,i,j)
                enddo
             enddo
          enddo
       enddo
    end do
  end function j2lso


end program
        ! uio=free_unit()
        ! if(master) open(uio,file='bracket_xmu'//str(iloop,Npad=3)//'.dat')
        ! if(iloop>1) then
        !    xmu1=xmu_imp
        !    xmu2=xmu1
        ! end if
        ! !+- check where to move - +!
        ! ntest1=get_delta_dens_imp(xmu1)       
        ! ntest2=ntest1
        ! alpha_mu= abs(alpha_mu)  !+- assuming starting from a system short of particles
        ! if(ntest1.gt.0) alpha_mu=-1.d0*abs(alpha_mu)  !+- if not adjust
        ! imu=0
        ! do while(ntest2*ntest1.ge.0.d0.and.imu.le.100)
        !    xmu1=xmu2
        !    ntest1=ntest2
        !    !
        !    xmu2=xmu2+alpha_mu
        !    ntest2=get_delta_dens_imp(xmu2)  
        !    if(master) write(uio,'(5F18.10)') xmu2,ntest2,ntest1
        !    imu=imu+1
        ! end do
        ! if(master) close(uio)
        ! xmu_imp = brentq(get_delta_dens_imp,xmu1,xmu2)            
        ! xmu=xmu_imp
        ! xmu_imp=50
        ! do imu=1,100
        !    ntest1=get_delta_dens_imp_(xmu_imp)
        !    if(master) write(450+iloop,*) xmu_imp,ntest1
        !    xmu_imp=xmu_imp-1.0
        ! end do
        ! if(iloop.eq.2) then
        !    call mpi_barrier(comm,mpiERR)
        !    stop
        ! end if
        ! xmu=1
