program ed_bilayer
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer :: ip,im,ipp,iw
  real(8) :: wx
  real(8) :: Vex
  real(8) :: Ef
  real(8) :: h0
  !
  real(8) :: ecut,pcut
  integer :: lp,lphi  
  integer :: unit_io
  !
  real(8),dimension(:),allocatable :: wr,modp,phip  
  complex(8),dimension(:,:),allocatable :: Pi_irreducible,Pi_irreducible_tmp
  !
  complex(8),dimension(:,:),allocatable :: Sigma,Sigma_tmp
  !
  complex(8),dimension(:,:),allocatable :: Lambda_vertex
  !
  complex(8),dimension(:),allocatable :: int_tmp,int_tmp_inn
  real(8),dimension(:),allocatable :: epp_tmp
  complex(8) :: wcmplx
  real(8) :: wp,mx,mel,bwp
  !  
  character(len=16)                           :: finput
  !MPI Vars:
  integer                                     :: comm,rank,mpierr,mpiSize
  logical                                     :: master
  !
  !
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  mpiSize = get_Size_MPI(comm)
  !
  !

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='input_Xpol.in')  
  call parse_input_variable(wx,"wx",finput,default=1.d0,comment='q=0  exciton frequency ')
  call parse_input_variable(Vex,"Vex",finput,default=1.d0,comment='e-X interaction')
  call parse_input_variable(ef,"ef",finput,default=0.d0,comment='fermi energy of the 2d electron gas')
  
  call parse_input_variable(ecut,"ecut",finput,default=1.d0,comment='cut off energy for momentum integration [eV]')  
  call parse_input_variable(lp,"lp",finput,default=100,comment='number of points for momenutm integration (modulus)')
  call parse_input_variable(lphi,"lphi",finput,default=100,comment='number of points for momenutm integration (phase)')


  call parse_input_variable(mX,"mX",finput,default=1.d0,comment='mass enhancement exciton')
  call parse_input_variable(mel,"mel",finput,default=1.d0,comment='mass enhancement electron')

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
  !
  !
  h0=7.62d-02 ! hbar^2/m [eV x nm^2]
  pcut=sqrt(2.d0*ecut/h0)
  !
  ! 
  allocate(wr(lreal));  wr=linspace(wini,wfin,lreal)
  allocate(modp(lp));   modp=linspace(0.d0,pcut,lp)
  allocate(phip(lphi)); phip=linspace(0.d0,2.d0*pi,lphi)
  !
  allocate(Pi_irreducible(lreal,lp));     Pi_irreducible = 0.d0
  allocate(Pi_irreducible_tmp(lreal,lp)); Pi_irreducible_tmp = 0.d0
  !  
  allocate(int_tmp(lp)); int_tmp=0.d0
  allocate(int_tmp_inn(lphi)); int_tmp_inn=0.d0

  allocate(epp_tmp(lphi)); epp_tmp=0.d0

  do iw=1+rank,lreal,mpiSize     
     wcmplx=wr(iw) + xi*eps
     
     if(master) write(*,*) iw
     
     do ip=1,lp
        !
        Pi_irreducible_tmp(iw,ip) = 0.d0
        !
        int_tmp=0d0
        do ipp=1,lp
           !
           int_tmp_inn = 0.d0
           epp_tmp(1:lphi) = h0*0.5*(modp(ip)**2.d0+modp(ipp)**2.d0-2.d0*modp(ipp)*modp(ip)*dcos(phip(1:lphi)))/mel-Ef
           wp = wx + h0*0.5*0.5*modp(ipp)**2.d0/mx        
           bwp = 1.d0/(exp(-beta*wp)-1.d0)
           int_tmp_inn(1:lphi) = (fermi(epp_tmp,beta) + bwp)/(wcmplx - epp_tmp - wp)*modp(ipp)/(2d0*pi)/(2d0*pi)
           !
           int_tmp(ipp) = trapz(int_tmp_inn,phip(1),phip(lphi))           
           !
        end do
        Pi_irreducible_tmp(iw,ip) = -1.d0*trapz(int_tmp,modp(1),modp(lp))
        !
     end do
  end do
  !
  call mpi_allreduce(Pi_irreducible_tmp,Pi_irreducible,lreal*lp,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  
  allocate(Lambda_vertex(lreal,lp));     Lambda_vertex = 0.d0

  Lambda_vertex=Vex*Vex*Pi_irreducible/(1.0-Vex*Pi_irreducible)

  if(master) then     
     unit_io=free_unit()
     open(unit=unit_io,file="pi_irreducible.out")
     do iw=1,lreal        
        do ip=1,lp
           write(unit_io,'(5F18.10)') wr(iw),modp(ip),Pi_irreducible(iw,ip)
        end do
        write(unit_io,*)
     end do
     close(unit_io)
     !
     open(unit=unit_io,file="Lambda_vertex.out")
     do iw=1,lreal        
        do ip=1,lp
           write(unit_io,'(5F18.10)') wr(iw),modp(ip),lambda_vertex(iw,ip)
        end do
        write(unit_io,*)
     end do
     close(unit_io)

     open(unit=unit_io,file="q0_pi_irreducible.out")
     do iw=1,lreal        
        write(unit_io,'(5F18.10)') wr(iw),Pi_irreducible(iw,1)
     end do     
     close(unit_io)

     open(unit=unit_io,file="q0_lambda_vertex.out")
     do iw=1,lreal        
        write(unit_io,'(5F18.10)') wr(iw),Lambda_vertex(iw,1)
     end do     
     close(unit_io)

  end if
  !
  allocate(Sigma(lreal,lp));     Sigma = 0.d0
  allocate(Sigma_tmp(lreal,lp)); Sigma_tmp = 0.d0
  !
  

  

  ! Nlat=2
  ! if(Nspin/=1.OR.Nlat/=2)stop "Wrong setup from input file: Nspin=1, Norb=2 -> 2Spin-Orbitals"
  ! Nso=Nspin*Norb
  ! Nlso=Nlat*Nso
  
  ! !Allocate Weiss Field:

  ! allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  ! allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  ! allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  ! allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  ! allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  ! allocate(Gtest(Nlat,Lmats))
  ! allocate(dens(Norb)); allocate(dens_prev(Norb))
  ! allocate(docc(Nlat)); 
  ! allocate(SigmaHk(Nso,Nso))
  ! allocate(Zmats(Nso,Nso))

  ! !Buil the Hamiltonian on a grid or on  path
  ! call set_SigmaHk()
  
  ! !+- a=1 is the hexagon edge
  ! !+- a' = sqrt(3.0)*a is the modulus of the primitive vector of the triangular lattice
  ! !+- e_1 = sqrt(3.0)/2.0 * a [-1, sqrt(3.0)]
  ! !+- e_2 = sqrt(3.0)/2.0 * a [ 1, sqrt(3.0)]
  ! e1 = sqrt(3d0)/2d0*[-1d0, sqrt(3d0)]
  ! e2 = sqrt(3d0)/2d0*[ 1d0, sqrt(3d0)]
  
  ! !RECIPROCAL LATTICE VECTORS:
  ! bklen=2d0*pi/3d0
  ! bk1=bklen*[ -sqrt(3d0), 1d0]
  ! bk2=bklen*[  sqrt(3d0), 1d0]
  ! call TB_set_bk(bkx=bk1,bky=bk2) 
  ! call build_hk_honeycomb()
  ! !
  ! !
  ! Nb=ed_get_bath_dimension()
  ! allocate(Bath(Nlat,Nb))
  ! allocate(Bath_prev(Nlat,Nb))
  ! call ed_init_solver(comm,bath)
  ! Bath_prev=Bath
  
  ! !
  ! !DMFT loop
  ! iloop=0;converged=.false.;converged0=.false.
  ! dens=[ntop,nbot]
  ! dens_prev=dens
  ! uio=free_unit()
  ! if(master) then
  !    open(uio,file='init_ndens_hf.out')
  !    write(uio,'(3F18.10)') dens
  !    close(uio)
  ! end if
  ! !
  ! uio=free_unit()
  ! if(master) then
  !    open(uio,file='iloop_observables.out')
  !    close(uio)
  ! end if
  ! call set_Hloc(dens)    
  ! !
  ! !
  ! xmu1=0.d0
  ! xmu2=xmu1


  ! xmu_imp=xmu
  ! do while(.not.converged.AND.iloop<nloop)
  !    iloop=iloop+1
  !    call start_loop(iloop,nloop,"DMFT-loop")
     
  !    !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
  !    call set_Hloc(dens,'Hloc_set'//str(iloop,Npad=3)//'.dat')     
     
  !    !+- set fix_mu=F and skip all this part
  !    if(fix_mu) then
  !       !   bracketing the xmu  !

  !       select case(fix_mu_scheme)
  !       case('newton')
  !          call newton(get_delta_dens_imp_,xmu,eps=1d-4)
  !       case('f_zero')
  !          xmu1= -50
  !          xmu2=  50
  !          call fzero(get_delta_dens_imp,xmu1,xmu2,imu,tol_rel=1.d-8,tol_abs=1.d-6)
  !          xmu=xmu1
  !       case default           
  !          xmu1= -50
  !          xmu2=  50
  !          call fzero(get_delta_dens_imp,xmu1,xmu2,imu,tol_rel=1.d-8,tol_abs=1.d-6)
  !          xmu=xmu1
  !       end select
  !       !
  !    end if
  !    !
  !    call ed_solve(comm,bath,Hloc,mpi_lanc=flag_mpi)     
  !    call ed_get_sigma_matsubara(Smats,Nlat)
  !    call ed_get_sigma_realaxis(Sreal,Nlat)

  !    !
  !    !
  !    call ed_get_dens(dens,Nlat,iorb=1)
  !    call ed_get_docc(docc,Nlat,iorb=1)
  !    !
  !    uio=free_unit()
  !    if(master) then
  !       open(uio,file='iloop_observables.out',status='old',position='append')
  !       write(uio,'(10F18.10)') dens,docc,xmu
  !       close(uio)
  !    end if
     
  !    ! !Get GLOC:
     
  !    !+- here I should use Hk_loop !!!!
  !    if(.not.allocated(hk_loop)) then
  !       call mpi_barrier(comm,mpiERR)
  !       stop
  !    end if
  !    call dmft_gloc_matsubara(Hk_loop,Gmats,Smats)
  !    call dmft_gloc_realaxis(Hk_loop,Greal,Sreal)
     
  !    ! !Update WeissField:
  !    call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,cg_scheme)
  !    !
  !    if(printG) then
  !       call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4,ineq_pad=2)
  !       call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4,ineq_pad=2)
  !       call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1,ineq_pad=2)
  !    end if               
  !    !
  !    call ed_chi2_fitgf(bath,Weiss,Hloc,ispin=1)  !+- perche qui si mangia anche Hloc?
  !    !
  !    !+- mix things
  !    Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
  !    dens = wmixing_dens*dens + (1.d0-wmixing_dens)*dens_prev

  !    Gtest=Weiss(:,1,1,1,1,:)
  !    if(.not.conv_dens) then
  !       converged = check_convergence(Gtest,dmft_error,nsuccess,nloop)
  !    else
  !       converged = check_convergence_local(dens,dens_error,nsuccess,nloop)
  !       ! dens_check=abs(dens(1)-dens_prev(1))**2.d0+abs(dens(2)-dens_prev(2))**2.d0
  !       ! if(dens_check.lt.dens_err) converged0=.true.
  !       ! if(dens_check.lt.dens_err.and.converged0) converged=.true.
  !    end if
  !    !
  !    Bath_prev = Bath
  !    dens_prev = dens
  !    xmu_imp   = xmu
  !    !
  !    !if(nread/=0d0)call ed_search_variable(xmu,sum(dens),converged)
     
  !    call end_loop

  ! enddo
  ! !
  ! call set_Hloc(dens,'Hloc_last.dat')     
  ! !
  ! call dmft_gloc_realaxis(Hk_loop,Greal,Sreal)
  ! call dmft_kinetic_energy(Hk_loop,Smats)
  ! !
  ! if(printG) then
  !    call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)
  !    call save_array("Smats",Smats)
  !    call save_array("Sreal",Sreal)
  ! end if
  !
  call finalize_MPI()

contains
  


end program
