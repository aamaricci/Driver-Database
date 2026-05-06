program ed_hm_square
  USE EDIPACK2
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                       :: iloop,Nb,Lk,Nx,Nso,ik,iorb,irepl
  logical                                       :: converged
  real(8)                                       :: wband,wmixing,onsite
  real(8),dimension(5)                          :: ts,Dband
  real(8),dimension(:),allocatable              :: dens
  !Bath:
  real(8),allocatable                           :: Bath(:),Bath_(:)
  real(8),dimension(:,:),allocatable            :: lambdasym_vectors
  complex(8),dimension(:,:,:,:,:),allocatable   :: Hsym_basis
  !The local hybridization function:
  complex(8),allocatable                        :: Hloc(:,:,:,:)
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Weiss,Weiss_
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gkmats
  complex(8),allocatable,dimension(:)           :: Gtest
  !Luttinger invariants:
  real(8),allocatable                           :: luttinger(:)
  !
  character(len=16)                             :: finput,foutput
  complex(8),allocatable                        :: Hk(:,:,:)
  real(8),allocatable                           :: Wt(:)
  !
  integer                                       :: comm,rank
  logical                                       :: master
  logical                                       :: mixG0,symOrbs
  !
  real(8)                                       :: Vel

  call init_MPI(comm,.true.)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputHM.in')
  call parse_input_variable(wmixing,"wmixing",finput,default=0.5d0,comment="Mixing bath parameter")
  call parse_input_variable(ts,"TS",finput,default=[0.25d0,0d0,0d0,0d0,0d0],comment="hopping parameter")
  call parse_input_variable(Dband,"Dband",finput,default=[0d0,0d0,0d0,0d0,0d0],comment="cystal field splittig (bands shift)")
  call parse_input_variable(Nx,"Nx",finput,default=100,comment="Number of kx point for 2d BZ integration")
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  call parse_input_variable(symOrbs,"symOrbs",finput,default=.false.)
  call parse_input_variable(Vel,"Vel",finput,default=0d0,comment="hopping parameter")
  !
  call ed_read_input(trim(finput))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")
  call add_ctrl_var(nbath,"nbath")
  call add_ctrl_var(ed_hw_bath,"ed_hw_bath")

  if(Nspin/=1.OR.Norb>5)stop "Wrong setup from input file: Nspin/=1 OR Norb>5"
  Nso=Nspin*Norb

  !Allocate Fields:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats),Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats),Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  allocate(dens(Norb))
  allocate(Gtest(Lmats))

  !Build Hk
  call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
  Lk = Nx*Nx
  allocate(Hk(Nso,Nso,Lk))
  call TB_build_model(Hk(:,:,:),hk_model,Nso,[Nx,Nx])
  Hloc   = zero
  Hloc(1,1,:,:) = sum(Hk,dim=3)/Lk
  where(abs(dreal(Hloc))<1.d-6) Hloc=0d0



  !Setup Bath
  select case(bath_type)
  case default
     stop "Wrong setup from input file: bath_type has to be 'normal' or 'replica'"
  case("normal")
     !
     Nb=ed_get_bath_dimension()
     !
  case("replica")
     !
     allocate(lambdasym_vectors(Nbath,1)) !Nsym=1
     allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,1))
     !
     Hsym_basis(:,:,:,:,1) = so2nn(zeye(Nso)) !Replica onsite energy
     !
     write(*,*) "Replica initialization: ed_hw_bath="//str(ed_hw_bath)
     !
     do irepl=1,Nbath
        onsite = irepl -1 - (Nbath-1)/2d0        ![-(Nbath-1)/2:(Nbath-1)/2]
        onsite = onsite * 2*ed_hw_bath/(Nbath-1) !P-H symmetric band, -ed_hw_bath:ed_hw_bath
        lambdasym_vectors(irepl,1) = onsite      !Multiplies the suitable identity 
     enddo
     !
     if(mod(Nbath,2)==0)then
        lambdasym_vectors(Nbath/2,1) = -1d-1    !Much needed small energies around
        lambdasym_vectors(Nbath/2+1,1) = 1d-1   !the fermi level. (for even Nbath)
     endif
     !
     call ed_set_Hreplica(Hsym_basis,lambdasym_vectors)
     Nb=ed_get_bath_dimension(1)
     !
  end select
  !


  call ed_set_hloc(Hloc)


  allocate(bath(Nb))
  allocate(bath_(Nb))
  bath_ = zero
  !
  !Setup solver
  call ed_init_solver(bath)




  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath) 
     call ed_get_sigma(Smats,axis='mats')
     call ed_get_sigma(Sreal,axis='real')

     !Compute the local gfs on the imaginary axis:
     call dmft_get_gloc(Hk,Gmats,Smats,axis='mats')
     call dmft_write_gf(Gmats,"Gloc",axis='mats',iprint=1)


     !Get the Weiss field/Delta function to be fitted
     select case(cg_scheme)
     case default;stop "cg_scheme unsupported"
     case('weiss')
        call dmft_self_consistency(Gmats,Smats,Weiss)
        call dmft_write_gf(Weiss,"Weiss",axis='mats',iprint=1)
     case('delta')
        call dmft_self_consistency(Gmats,Smats,Weiss,Hloc)
        call dmft_write_gf(Weiss,"Delta",axis='mats',iprint=1)
     end select


     if(mixG0)then
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
        Weiss_=Weiss
     endif

     !Perform the SELF-CONSISTENCY by fitting the new bath
     if(symOrbs)then
        call ed_chi2_fitgf(Weiss,bath,ispin=1,iorb=1)
        call ed_orb_equality_bath(bath,save=.true.)
     else
        call ed_chi2_fitgf(Weiss,bath,ispin=1)
     endif

     !MIXING:
     if(.not.mixG0)then
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
        Bath_=Bath
     endif

     !Check convergence (if required change chemical potential)     
     Gtest=zero
     do iorb=1,Norb
        Gtest=Gtest+Weiss(1,1,iorb,iorb,:)/Norb
     enddo
     converged = check_convergence(Gtest,dmft_error,nsuccess,nloop,reset=.false.)
     call ed_get_dens(dens)
     if(nread/=0d0)call ed_search_variable(xmu,sum(dens),converged)

     call end_loop
  enddo


  !Compute the local gfs on the real axis:
  call dmft_get_gloc(Hk,Greal,Sreal,axis='real')
  call dmft_write_gf(Greal,"Gloc",axis='real',iprint=1)


  !Get kinetic energy:
  call dmft_kinetic_energy(Hk,Smats)



  ! !Compute the Luttinger invariants:
  ! if(master)then
  !    allocate(luttinger(Norb))
  !    call luttinger_integral(luttinger,Greal,Sreal)
  !    call print_luttinger(luttinger)
  !    deallocate(luttinger)
  ! endif

  call finalize_MPI()

contains

  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Hk model for the 2d square lattice
  !-------------------------------------------------------------------------------------------
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:) :: kpoint
    integer              :: N,ih
    real(8)              :: kx,ky
    complex(8)           :: hk(N,N)
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = zero
    do ih=1,N
       Hk(ih,ih) = -one*2d0*ts(ih)*(cos(kx)+cos(ky)) + Dband(ih)
    enddo
  end function hk_model

  !+---------------------------------------------------------------------------+
  !PURPOSE : Compute the Luttinger integral IL = 1/π * |Im(∫dw(Gloc*dSloc/dw)|
  !          > Cfr. J. Phys. Cond. Mat. 28, 02560 and PRB B 102, 081110(R)
  ! NB) At ph-symmetry IL should be zero, but only if the T->0 limit is taken 
  !     before the mu->0 limit, and here we apparently invert the order of the 
  !     limits because of the numerical discretization on the imaginary axis. 
  !     > Look at [53] (report of a private communication) in the Phys. Rev. B
  !+---------------------------------------------------------------------------+
  subroutine luttinger_integral(IL,Gloc,Sloc)
    real(8),allocatable,dimension(:),intent(out)            :: IL
    complex(8),allocatable,dimension(:,:,:,:,:),intent(in)  :: Gloc,Sloc
    real(8),dimension(Lreal)                                :: dSreal,dSimag,integrand
    integer                                                 :: iorb,ispin
    !
    if(.not.allocated(IL)) allocate(IL(Norb))
    !
    do iorb=1,Norb
       do ispin=1,Nspin !Nspin is constrained to be equal to 1 above.
          dSreal(:) = deriv(dreal(Sloc(ispin,ispin,iorb,iorb,:)),1d0) !No need to include 1/dw if
          dSimag(:) = deriv(dimag(Sloc(ispin,ispin,iorb,iorb,:)),1d0) !we are integrating in dw...
          integrand = dimag(Gloc(ispin,ispin,iorb,iorb,:)*cmplx(dSreal,dSimag))
          IL(iorb) = sum(integrand)       !Naive Riemann-sum (no need of trapezoidal-sums with typical Lreal)
          IL(iorb) = 1/pi * abs(IL(iorb)) !The theoretical sign is determined by sign(mu)...
       enddo
    enddo
    !
  end subroutine luttinger_integral

  !+---------------------------------------------------------------------------+
  !PURPOSE : print to file (and stdout) the Luttinger invariants IL(Norb)
  !+---------------------------------------------------------------------------+
  subroutine print_luttinger(IL)
    real(8),allocatable,dimension(:),intent(in)  :: IL
    integer                                      :: iorb
    integer                                      :: unit
    !
    if(ed_verbose>0)then
       write(LOGfile,*) " "
       write(LOGfile,*) "Luttinger invariants:"
    endif
    !
    do iorb=1,Norb
       if(ed_verbose>0) write(LOGfile,*) iorb, IL(iorb)
       unit = free_unit()
       foutput = "luttinger_l"//str(iorb)//".dat"
       open(unit,file=foutput,action="write",position="rewind",status='unknown')
       write(unit,*) iorb, IL(iorb)
       close(unit)
    enddo
    !
    if(ed_verbose>0)then
       write(LOGfile,*) "          iorb        IL"
       write(LOGfile,*) " "
    endif
    !
  end subroutine print_luttinger



  !+---------------------------------------------------------------------------+
  !PURPOSE : reshape from [Nso]x[Nso] matrix to [Nspin,Nspin,Norb,Norb] array
  !+---------------------------------------------------------------------------+
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

  !+---------------------------------------------------------------------------+
  !PURPOSE : reshape from [Nspin,Nspin,Norb,Norb] array to [Nso]x[Nso] matrix
  !+---------------------------------------------------------------------------+
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

end program ed_hm_square





