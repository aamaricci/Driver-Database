<<<<<<< HEAD
program ed_hm_2b_cubic
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none  
  integer                                     :: iloop,Nb,Nk,Nso,Nktot
  logical                                     :: converged
  real(8)                                     :: wband,ts,de
  real(8),allocatable                         :: dens(:)
  !
  !Bath:
  real(8),allocatable,dimension(:)            :: Bath,Bath_Prev
  !
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:) :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Weiss_
  real(8)                                     :: crystal_field,var,wmixing
  !
  complex(8),allocatable,dimension(:,:,:)     :: Hk
  complex(8),dimension(:,:,:,:),allocatable   :: Hloc




  
  !Local, user defined variables. Parse_input (in Scifor) add them to the list
  call parse_input_variable(wband,"WBAND","inputED.conf",default=1.d0)
  call parse_input_variable(Nk,"Nk","inputED.conf",default=10)
  call parse_input_variable(crystal_field,"CRYSTAL_FIELD","inputED.conf",default=0.d0)
  call parse_input_variable(wmixing,"WMIXING","inputED.conf",default=1.d0)
  !ED read input
  call ed_read_input("inputED.conf")
  !Add ctrl variables to DMFT_TOOLS memory pool
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  !                               
  !
  Nso = Nspin*Norb
  !
  !
  ts=wband/6.d0

  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));Gmats=0.d0
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));Smats=0.d0
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats));Weiss=0.d0
  allocate(Weiss_(Nspin,Nspin,Norb,Norb,Lmats));Weiss_=0.d0
  !
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal));Greal=0.d0
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal));Sreal=0.d0
  !
  allocate(dens(Norb));

  !+- the non-interacting hamiltonian
  !Here Hk and Hloc are allocated and built
  call build_Hk
  
  !+- Setup the Solver -+!
  !Get the dimension of the bath, user side.
  !this is evaluated internally by code given all the inputs
  Nb=ed_get_bath_dimension()
  !
  allocate(Bath(Nb))            !actual bath
  allocate(Bath_prev(Nb))       !prev bath, for mixing
  call ed_init_solver(bath)


  !+- DMFT loop -+!
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath,Hloc)
     !

     !Retrieve impurity self-energies 
     call ed_get_Sigma_matsubara(Smats)
     call ed_get_Sigma_realaxis(Sreal)

     !
     !Compute the local greens function of the lattice
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     !

     !+- get the new Weiss field by imposing that G_{lattice} = G_{impurity}
     call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,trim(cg_scheme))
     !
     !+- printout the local GF and the Weiss field
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)
     !
     !

     !+- get the new bath parameters
     call ed_chi2_fitgf(Weiss,bath,ispin=1)

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath

     converged = check_convergence(Weiss(1,1,1,1,:)+Weiss(1,1,2,2,:),dmft_error,nsuccess,nloop)

     !
     !+- retrieve the density and check for chemical potential adjustments
     call ed_get_dens(dens)
     if(nread/=0.d0) call ed_search_variable(xmu,dens(1)+dens(2),converged)

  enddo
  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)
  !
  call dmft_kinetic_energy(Hk,Sreal)
  !

contains



  !---------------------------------------------------------------------
  !GET 2bands cubic lattice HAMILTONIAN (from the NonInteracting code)
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j
    !
    !Set the reciprocal lattice basis vector
    call TB_set_bk(bkx=[pi2,0d0,0d0],bky=[0d0,pi2,0d0],bkz=[0d0,0d0,pi2])
    !
    write(LOGfile,*)"Build H(k):"
    Nktot=Nk**3
    write(*,*)"# of k-points     :",Nktot
    write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    if(allocated(Hloc))deallocate(Hloc)
    allocate(Hk(Nso,Nso,Nktot)) ;Hk=zero
    !
    call TB_build_model(Hk,hk_2bands,Nso,[Nk,Nk,Nk])
    ! if(present(file))then
    call TB_write_hk(Hk,"Hk.dat",&
         Nlat=1,&
         Nspin=1,&
         Norb=Norb,&
         Nkvec=[Nk,Nk,Nk])
    ! endif
    allocate(Hloc(Nspin,Nspin,Norb,Norb));  Hloc = 0d0
    Hloc = j2so(sum(Hk,dim=3))/Nktot !send a 4dim array into a 2dim one
    where(abs(dreal(Hloc))<1d-6)Hloc=zero
    ! call TB_write_Hloc(so2j(Hloc))
  end subroutine build_hk


  !--------------------------------------------------------------------!
  !2bands cubic lattice HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_2bands(kvec,N) result(hk)
    integer                   :: N
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: ek,kx,ky,kz
    if(N/=2)stop "hk_2bands error: N != Nspin*Norb == 2"
    kx=kvec(1)
    ky=kvec(2)
    kz=kvec(3)
    ek = -2.d0*ts*(cos(kx)+cos(ky)+cos(kz))
    Hk = ek*eye(N)              !build the non-local part of the H(k)
    Hk(1,1) = Hk(1,1) + crystal_field*0.5d0
    Hk(2,2) = Hk(2,2) - crystal_field*0.5d0
  end function hk_2bands





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


end program ed_hm_2b_cubic




=======
program ed_hm_2b_cubic
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none  
  integer                :: iloop,Nb,Le,ie,mu_loop,Nk
  logical                :: converged,dos_3d
  real(8),allocatable    :: wm(:),wr(:)
  real(8)                :: wband,ts,de
  real(8),allocatable    :: dens(:),docc(:) 
  !Bath:
  real(8),allocatable,dimension(:)              :: Bath,Bath_Prev
  !The local hybridization function:

  complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:) :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss

  complex(8),allocatable :: Delta(:,:,:)
  real(8),allocatable    :: epsik(:),wt(:)
  real(8)                :: crystal_field,var,wmixing
  integer                :: unit,iorb
  
  complex(8),allocatable,dimension(:,:,:) :: Hk
  complex(8),dimension(:,:,:,:),allocatable     :: Hloc
  !
  call parse_input_variable(wband,"WBAND","inputED.in",default=1.d0)
  !call parse_input_variable(Le,"NE","inputED.in",default=2000)
  call parse_input_variable(Nk,"Nk","inputED.in",default=10)
  !call parse_input_variable(dos_3d,"DOS","inputED.in",default=.false.)
  call parse_input_variable(crystal_field,"CRYSTAL_FIELD","inputED.in",default=0.d0)
  call parse_input_variable(wmixing,"WMIXING","inputED.in",default=1.d0)
  call ed_read_input("inputED.in")

  
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  !                               
  !
  !
  !
  !+- allocate k-grids -+! 
  ts=wband/6.d0
  ! if(dos_3d) then
  !    call get_cubic_dos(100000000)
  ! else
  !    call get_cubic_k(Nk)
  ! end if
  !
  !

  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));Gmats=0.d0
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));Smats=0.d0
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats));Weiss=0.d0
  !
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal));Greal=0.d0
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal));Sreal=0.d0

  
  !
  !+- the local hamiltonian
  allocate(Hloc(Nspin,Nspin,Norb,Norb));  Hloc = 0d0
  Hloc(1,1,1,1) =  crystal_field*0.5d0
  Hloc(1,1,2,2) = -crystal_field*0.5d0
  !+- the bare bands
  call get_cubic_k(Nk)  
  call build_Hk

  !
  allocate(dens(Norb));
  !
  
  !
  !+- Setup the Solver -+!
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_prev(Nb))
  call ed_init_solver(bath)
  
  
  !+- DMFT loop -+!
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath,Hloc)
     !

     !Retrieve impurity self-energies 
     call ed_get_Sigma_matsubara(Smats)
     !
     !Compute the local greens function of the lattice
     call dmft_gloc_matsubara(Hk,Gmats,Smats)

     !+- get the new Weiss field by imposing that G_{lattice} = G_{impurity}
     call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,trim(cg_scheme))
     !
     !+- printout the local GF and the Weiss field
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)
     
     !+- get the new bath parameters
     call ed_chi2_fitgf(Weiss,bath,ispin=1)

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath

     converged = check_convergence(Weiss(1,1,1,1,:)+Weiss(1,1,2,2,:),dmft_error,nsuccess,nloop,reset=.false.)
     
     !
     !+- retrieve the density and check for chemical potential adjustments
     call ed_get_dens(dens)
     if(nread/=0.d0) call ed_search_variable(xmu,dens(1)+dens(2),converged)
     
  enddo
  !Get Kinetic Energy too
!   if(ED_MPI_ID==0) call ed_kinetic_energy(impSmats(1,1,:,:,:),Hk,wt)

! #ifdef _MPI
!   call MPI_FINALIZE(ED_MPI_ERR)
! #endif

contains

  ! subroutine get_delta
  !   integer                     :: i,j,ie,iorb,jorb
  !   complex(8)                  :: iw,zita,g0and,g0loc,gg
  !   complex(8),dimension(Lmats) :: self
  !   complex(8),dimension(Lreal) :: selfr,grloc
  !   !
  !   complex(8),allocatable      :: GLoc_mats(:,:,:),GLoc_real(:,:,:)
  !   complex(8),allocatable      :: tmp_gloc(:,:)
  !   !
  !   allocate(GLoc_mats(Norb,Norb,Lmats),GLoc_real(Norb,Norb,Lreal))
  !   allocate(tmp_gloc(Norb,Norb))    
  !   !+- MATSUBARA FREQ -+!
  !    if(ED_MPI_ID==0) print*,"Get Gloc_iw:"
  !   delta = zero
  !   do i=1,Lmats
  !      iw = xi*wm(i)+xmu
  !      !+- compute loacl greens function -+!
  !      GLoc_mats(:,:,i)=zero
  !      do ie=1,Le
  !         GLoc_mats(:,:,i) = GLoc_mats(:,:,i) +  gk(iw,impSmats(1,1,:,:,i),Hk(:,:,ie))*wt(ie)
  !      end do
  !      if(cg_scheme=='weiss') then
  !         tmp_gloc=GLoc_mats(:,:,i)
  !         call matrix_inverse(tmp_gloc)
  !         delta(:,:,i) = tmp_gloc + impSmats(1,1,:,:,i)
  !         call matrix_inverse(delta(:,:,i))
  !      else
  !         tmp_gloc=GLoc_mats(:,:,i)
  !         call matrix_inverse(tmp_gloc)
  !         delta(:,:,i) = inverse_g0imp(iw) - tmp_gloc - impSmats(1,1,:,:,i)
  !      end if
  !   end do
  !   !+- REAL FREQ -+!
  !   if(ED_MPI_ID==0) print*,"Get Gloc_real:"
  !   do i=1,Lreal
  !      iw=cmplx(wr(i),eps)+xmu
  !      GLoc_real(:,:,i) = zero       
  !      do ie=1,Le
  !         GLoc_real(:,:,i) = GLoc_real(:,:,i) + gk(iw,impSreal(1,1,:,:,i),Hk(:,:,ie))*wt(ie)
  !      enddo
  !   enddo
  !   !
  !   if(ED_MPI_ID==0) then
  !      if(bath_type.eq.'normal') then
  !         do iorb=1,Norb
  !            call splot("Gloc_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_iw.ed",wm,GLoc_mats(iorb,iorb,:))
  !            call splot("Gloc_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_real.ed",wr,GLoc_real(iorb,iorb,:))
  !            call splot("DOS_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_real.ed",wr,-dimag(GLoc_real(iorb,iorb,:))/pi)
  !            call splot("Delta_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_iw.ed",wm,delta(iorb,iorb,:))
  !         end do
  !      else
  !         do iorb=1,Norb
  !            do jorb=1,Norb          
  !               call splot("Gloc_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_iw.ed",wm,GLoc_mats(iorb,jorb,:))
  !               call splot("Gloc_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_real.ed",wr,GLoc_real(iorb,jorb,:))
  !               call splot("Delta_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_iw.ed",wm,delta(iorb,jorb,:))
  !            end do
  !         end do
  !      end if
  !   end if
  !   deallocate(GLoc_mats,GLoc_real)

 
  ! end subroutine get_delta

  ! !+- 2BAND GREENS FUNCTION -+!  
  ! function g0k(iw,hk) 
  !   integer                     :: i
  !   complex(8),dimension(2,2)   :: hk
  !   complex(8)                  :: iw
  !   complex(8),dimension(2,2)   :: g0k,g0k_
  !   complex(8)                  :: det,ginv11,ginv22
  !   complex(8)                  :: vmix12,vmix21
  !   g0k =zero
  !   g0k_=zero
  !   !+- build G0^(-1) -+!
  !   g0k_(1,1) = iw 
  !   g0k_(2,2) = iw 
  !   g0k_=g0k_-hk
  !   det = g0k_(1,1)*g0k_(2,2) - g0k_(1,2)*g0k_(2,1)
  !   !+- compute inverse -+!
  !   g0k(1,1) =  g0k_(2,2)/det
  !   g0k(2,2) =  g0k_(1,1)/det
  !   g0k(1,2) = -g0k_(1,2)/det
  !   g0k(2,1) = -g0k_(2,1)/det
  ! end function g0k
  ! !
  ! function gk(iw,sigma,hk) 
  !   complex(8)                  :: iw
  !   complex(8),dimension(2,2)   :: sigma
  !   complex(8),dimension(2,2)   :: hk
  !   complex(8),dimension(2,2)   :: gk,gk_
  !   complex(8),dimension(2,2)   :: g0k
  !   complex(8)                  :: det
  !   gk  = zero
  !   gk_ = zero
  !   g0k = zero
  !   g0k = inverse_g0k(iw,hk)
  !   !+- build G^(-1) -+!
  !   gk_ = g0k-sigma
  !   !+- compute inverse -+!
  !   det = gk_(1,1)*gk_(2,2) - gk_(1,2)*gk_(2,1)
  !   gk(1,1) =  gk_(2,2)/det
  !   gk(2,2) =  gk_(1,1)/det
  !   gk(1,2) = -gk_(1,2)/det
  !   gk(2,1) = -gk_(2,1)/det
  ! end function gk

  ! !+- INVERSE 2BAND GREENS FUNCTION -+!
  ! function inverse_g0k(iw,hk) result(g0k)
  !   integer                     :: i
  !   complex(8),dimension(2,2)   :: hk
  !   complex(8)                  :: iw
  !   complex(8),dimension(2,2)   :: g0k
  !   complex(8)                  :: ginv11,ginv22
  !   complex(8)                  :: vmix12,vmix21
  !   g0k=zero
  !   !+- build G0^(-1) -+!
  !   g0k(1,1) = iw
  !   g0k(2,2) = iw
  !   g0k = g0k - hk
  ! end function inverse_g0k
  ! !
  ! function inverse_g0imp(iw) result(g0k)
  !   integer                     :: i
  !   complex(8),dimension(2,2)   :: hk
  !   complex(8)                  :: iw
  !   complex(8),dimension(2,2)   :: g0k
  !   complex(8)                  :: ginv11,ginv22
  !   complex(8)                  :: vmix12,vmix21
  !   g0k=zero
  !   !+- build G0_imp^(-1) -+!
  !   g0k(1,1) = iw
  !   g0k(2,2) = iw    
  !   g0k = g0k - Hloc(1,1,:,:)
  ! end function inverse_g0imp
  ! !
  ! function inverse_gk(iw,sigma,hk) result(gk)
  !   complex(8)                  :: iw
  !   complex(8),dimension(2,2)   :: sigma
  !   complex(8),dimension(2,2)   :: hk
  !   complex(8),dimension(2,2)   :: gk,gk_
  !   complex(8),dimension(2,2)   :: g0k
  !   complex(8)                  :: det
  !   gk  = zero
  !   g0k = zero
  !   g0k = inverse_g0k(iw,hk)
  !   !+- build G^(-1) -+!
  !   gk = g0k-sigma
  ! end function inverse_gk
  ! !

  !+- HK routines -+!
  subroutine build_hk
    integer  :: ik
    allocate(Hk(Norb,Norb,Le))
    Hk=0.d0
    do ik=1,Le
       call get_hk(Hk(:,:,ik),ik)
    end do
  end subroutine build_hk

  subroutine get_hk(hk,ik)
    complex(8),dimension(Norb,Norb),intent(inout) :: hk
    integer  :: ik
    hk=zero
    do iorb=1,Norb
       hk(iorb,iorb) = epsik(ik) 
    end do
    hk = hk + Hloc(1,1,:,:)
  end subroutine get_hk
  
  
  
  ! subroutine get_cubic_dos(Nrnd)
  !   integer :: Nrnd
  !   integer :: i,ik,unit
  !   integer(8),allocatable,dimension(:) :: count_
  !   real(8) :: kx,ky,kz,ek
  !   allocate(wt(Le),epsik(Le))
  !   allocate(count_(Le))
  !   epsik=linspace(-wband,wband,Le,mesh=de)
  !   count_=0
  !   call random_seed(put=[1234567])

  !   do i=1,Nrnd
  !      call random_number(kx)
  !      kx=(kx-1.d0)*pi
  !      call random_number(ky)
  !      ky=(ky-1.d0)*pi
  !      call random_number(kz)    
  !      kz=(kz-1.d0)*pi              
  !      ek=-2.d0*ts*(cos(kx)+cos(ky)+cos(kz))
  !      ik=1+ek/de + wband/de
  !      count_(ik)=count_(ik)+1   
  !   end do
  !   !+- symmetryze -+!
  !   do i=1,Le/2
  !      count_(i) = count_(Le-(i-1))
  !   end do
  !   wt=dble(count_)/dble(Nrnd) 
  !   wt=wt/trapz(de,wt)
  !   wt=wt*de
  !   !if(ED_MPI_ID==0) then
  !   unit=free_unit()
  !   open(unit,file='DOS.3d')
  !   do i=1,Le
  !      write(unit,*) epsik(i),wt(i)
  !   end do
  !   close(unit)
  !   !end if
  ! end subroutine get_cubic_dos

  subroutine get_cubic_k(Nk)
    integer :: Nk
    integer :: i,ik
    integer :: ix,iy,iz
    integer(8),allocatable,dimension(:) :: count_
    real(8),allocatable,dimension(:) :: kx
    real(8) :: tmp
    
    allocate(kx(Nk))
    kx=linspace(-pi,pi,Nk)
    Le=Nk*Nk*Nk
    allocate(wt(Le),epsik(Le))    
    ik=0
    do ix=1,Nk
       do iy=1,Nk
          do iz=1,Nk
             ik=ik+1
             epsik(ik)=-2.d0*ts*(cos(kx(ix))+cos(kx(iy))+cos(kx(iz)))
             wt(ik) = 1.d0/dble(Le)
          end do
       end do
    end do
    call get_free_dos(epsik,wt,file='DOS_free.kgrid')
  end subroutine get_cubic_k
  

end program ed_hm_2b_cubic




>>>>>>> 24c73b05581aca0b17af06ea30763e8ee244cf71
