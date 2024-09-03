!> Solve the 2D BHZ + Kanamori n.n + Disorder model using Mean-Field
program MF_bhz_2d_Anderson
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer,parameter                             :: Norb=2,Nspin=2,Nso=Nspin*Norb
  integer                                       :: Nk,Nktot,Nkx,Nky,Nkz,L
  integer                                       :: Nlat,Nx,Ny,Nz
  integer                                       :: Nkpath,Npts,z2(4)
  integer                                       :: i,j,k,ik
  integer                                       :: io,jo,ispin
  integer                                       :: iorb,ilat,jorb,jlat
  integer                                       :: ix,iy,iz
  real(8)                                       :: kx,ky,kz
  real(8),dimension(:),allocatable              :: Wtk,Evals,rhoDiag,Erandom
  real(8),dimension(:,:),allocatable            :: kpath,ktrims,Kgrid,Rgrid,Ldens
  integer,dimension(:,:),allocatable          :: Links
  complex(8),dimension(:,:,:),allocatable       :: Hij,HijMF
  complex(8),dimension(:,:,:,:),allocatable     :: Hlat
  integer                                       :: Iter,MaxIter,Nsuccess=2,Idum,Nblock
  real(8)                                       :: chern,Uloc,Jh,JU,Sz,Tz,Rz,Ntot
  real(8)                                       :: mh,lambda,delta,lz,Wdis
  real(8)                                       :: xmu,beta,eps,Ekin,Eloc
  real(8)                                       :: n(Nso),arg,dens(Nso),wmix,it_error,sb_field
  complex(8)                                    :: Hloc(Nso,Nso),Htmp(Nso,Nso)
  complex(8),dimension(:,:,:,:,:,:),allocatable :: GLmats,GLreal
  character(len=20)                             :: file
  logical                                       :: iexist,converged,with_mats_gf,with_real_gf,bool
  complex(8),dimension(Nso,Nso)                 :: Gamma0,Gamma5,GammaX,GammaY,GammaZ,GammaS
  real(8),dimension(3)                          :: vecK,vecRi,vecRj
  complex(8),dimension(:,:),allocatable         :: rhoH
  real(8),dimension(:,:),allocatable            :: params,params_prev,global_params
  integer                                       :: disorder_type !0=scalar disorder, 1=magnetic_disorder


  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(L,"L","inputBHZ.conf",default=2048)
  call parse_input_variable(Wdis,"WDIS","inputBHZ.conf",default=0d0)
  call parse_input_variable(idum,"IDUM","inputBHZ.conf",default=1234567)
  call parse_input_variable(nblock,"NBLOCK","inputBHZ.conf",default=2)
  call parse_input_variable(disorder_type,"DISORDER_TYPE","inputBHZ.conf",default=0)
  call parse_input_variable(Uloc,"ULOC","inputBHZ.conf",default=1d0)
  call parse_input_variable(Jh,"Jh","inputBHZ.conf",default=0.125d0)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=3.d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=1000.d0)
  call parse_input_variable(wmix,"WMIX","inputBHZ.conf",default=0.5d0)
  call parse_input_variable(sb_field,"SB_FIELD","inputBHZ.conf",default=0.01d0)
  call parse_input_variable(it_error,"IT_ERROR","inputBHZ.conf",default=1d-5)
  call parse_input_variable(maxiter,"MAXITER","inputBHZ.conf",default=100)
  call parse_input_variable(with_mats_gf,"WITH_MATS_GF","inputBHZ.conf",default=.false.)
  call parse_input_variable(with_real_gf,"WITH_REAL_GF","inputBHZ.conf",default=.false.)
  call save_input_file("inputBHZ.conf")
  call print_input()
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-10d0,"wini")
  call add_ctrl_var(10d0,"wfin")
  call add_ctrl_var(eps,"eps")

  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  !< Momentum space:
  Nky = Nkx
  Nktot=Nkx*Nky
  !
  !< Real space
  Nx = Nkx
  Ny = Nkx
  Nlat = Nx*Ny

  !SETUP THE GAMMA MATRICES:
  gamma0=kron_pauli( pauli_tau_0, pauli_sigma_0)
  gammaX=kron_pauli( pauli_tau_z, pauli_sigma_x)
  gammaY=kron_pauli( pauli_tau_0,-pauli_sigma_y)
  gammaZ=kron_pauli( pauli_tau_x, pauli_sigma_x)
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z)
  gammaS=kron_pauli( pauli_tau_z, pauli_sigma_0)

  !< Build disorder:
  allocate(erandom(Nlat))
  call mersenne_init(idum)
  call mt_random(erandom)
  erandom=(2d0*erandom-1d0)*Wdis/2d0
  inquire(file='erandom_'//str(idum)//'.restart',exist=bool)
  if(bool)then
     if(file_length('erandom_'//str(idum)//'.restart')/=Nlat)&
          stop "mf_bhz_2d_anderson error: found erandom.restart with length different from Nlat"
     call read_array('erandom_'//str(idum)//'.restart',erandom)
  endif
  call save_array('erandom_'//str(idum)//'.used',erandom)
  !



  !< start using TB_procedures:
  !< 1st set up the direct and momentum space lattice basis
  call TB_set_ei([1d0,0d0],[0d0,1d0])
  call TB_set_bk([pi2,0d0],[0d0,pi2])


  allocate(Hlat(Nso,Nso,Nlat,Nlat))
  allocate(Links(4,2))
  !
  Links(1,:) = [1,0]
  Links(2,:) = [0,1]
  Links(3,:) = [-1,0]
  Links(4,:) = [0,-1]
  print*,"Building up Hlat model"
  call TB_build_model(Hlat,ts_model,Nso,[Nkx,Nky],Links,pbc=.true.)
  !
  !Add DISORDER TERM TO H_0
  select case(disorder_type)
  case(0)
     do ilat=1,Nlat
        Hlat(:,:,ilat,ilat) = Hlat(:,:,ilat,ilat) + erandom(ilat)*Gamma0
     enddo
  case(1)
     do ilat=1,Nlat
        Hlat(:,:,ilat,ilat) = Hlat(:,:,ilat,ilat) + erandom(ilat)*GammaS
     enddo
  case(2)
     do ilat=1,Nlat
        Hlat(:,:,ilat,ilat) = Hlat(:,:,ilat,ilat) + erandom(ilat)*Gamma5
     enddo
  case default
     stop "disorder_type not [0:2]"
  end select
  !
  !Dump Hlat into H* with shape(Nlso,Nlso,Nk==1) to be used later
  allocate(Hij(Nlat*Nso,Nlat*Nso,1))
  Hij = zero
  do ilat=1,Nlat
     do jlat=1,Nlat
        do io=1,Nso
           do jo=1,Nso
              i = io + (ilat-1)*Nso
              j = jo + (jlat-1)*Nso
              Hij(i,j,1) = Hlat(io,jo,ilat,jlat)
           enddo
        enddo
     enddo
  enddo




  !Start MF HERE:
  allocate(params(Nlat,2))
  allocate(params_prev(Nlat,2))
  do ilat=1,Nlat
     params(ilat,:)        = [sb_field,sb_field]   ![Tz,Sz]
  enddo
  inquire(file="params.restart",exist=iexist)
  if(iexist)then
     call read_array("params.restart",params)     
     params(:,2)=params(:,2)+sb_field
  endif
  call save_array("params.init",params)
  !
  !DO MEAN-FIELD CYCLE
  open(100,file="tz_sz.dat")
  open(101,file="dens_1.dat")
  open(102,file="dens_2.dat")
  converged=.false. ; iter=0
  do while(.not.converged.AND.iter<maxiter)
     iter=iter+1
     call start_loop(iter,maxiter,"MF-loop")
     !
     call symmetrize_params(params)
     call solve_MF_bhz(iter,params)
     if(iter>1)params = wmix*params + (1d0-wmix)*params_prev
     params_prev = params
     !
     converged = check_convergence_local(params,it_error,nsuccess,maxiter) 
     !
     call end_loop
  end do
  call save_array("params.restart",params)
  close(100)
  close(101)
  close(102)


  !< BUILD THE LOCAL GF
  if(with_mats_gf)then
     allocate(GLmats(Nlat,Nspin,Nspin,Norb,Norb,L))
     allocate(HijMF(Nlat*Nso,Nlat*Nso,1))
     HijMF(:,:,1) = get_Hmf(params,Hij)
     call dmft_gloc_matsubara(HijMF,GLmats,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
     call dmft_print_gf_matsubara(GLmats,"Gloc",iprint=1)
  endif

  if(with_real_gf)then
     allocate(GLreal(Nlat,Nspin,Nspin,Norb,Norb,L))
     allocate(HijMF(Nlat*Nso,Nlat*Nso,1))
     HijMF(:,:,1) = get_Hmf(params,Hij)
     call dmft_gloc_realaxis(HijMF,GLreal,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
     call dmft_print_gf_realaxis(GLreal,"Gloc",iprint=1)
  endif


  

  
contains




  subroutine solve_MF_bhz(iter,a)
    real(8),dimension(Nlat,2),intent(inout) :: a
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Hmat,rhoDiag
    real(8),dimension(Nlat*Nso)             :: Evals!,rhoDiag
    real(8),dimension(Nlat*Nso,Nlat*Nso)    :: rhoH!,rhoDiag
    real(8)                                 :: dens(Nlat,Nspin,Norb)
    integer                                 :: iter,iorb,ispin
    !
    Hmat = Hij(:,:,1) + mf_Hij_correction(a)
    !
#ifdef _W_SCALAPACK
    call p_eigh(Hmat,Evals,Nblock)       !diag Hij -> Ea_i
#else
    call eigh(Hmat,Evals)       !diag Hij -> Ea_i
#endif

    rhoDiag = diag(fermi(Evals,beta))
#ifdef _W_SCALAPACK
    rhoH    = ( Hmat.px.rhoDiag ).px.conjg(transpose(Hmat))
#else
    rhoH    = matmul(Hmat , matmul(rhoDiag, conjg(transpose(Hmat))) )
#endif
    !
    dens = 0d0        
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb+(ispin-1)*Norb+(ilat-1)*Nspin*Norb
             dens(ilat,ispin,iorb) = rhoH(io,io)
          enddo
       enddo
       a(ilat,1) = 0.5d0*sum(dens(ilat,:,1)) - 0.5d0*sum(dens(ilat,:,2)) !Tz = N_1  - N_2
       a(ilat,2) = 0.5d0*sum(dens(ilat,1,:)) - 0.5d0*sum(dens(ilat,2,:)) !Sz = N_up - N_dw
    enddo
    !
    rewind(100);rewind(101);rewind(102)
    do ilat=1,Nlat
       write(100,*)a(ilat,1),a(ilat,2)
       write(101,*)(dens(ilat,ispin,1),ispin=1,Nspin)
       write(102,*)(dens(ilat,ispin,2),ispin=1,Nspin)
    enddo
  end subroutine solve_MF_bhz


  function mf_Hij_correction(a) result(HijMF)
    real(8),dimension(Nlat,2)               :: a
    complex(8),dimension(Nlat,Nso,Nso)      :: Htmp
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: HijMF
    integer                                 :: ilat,io,jo,i,j
    HijMF=zero
    do ilat=1,Nlat
       Htmp(ilat,:,:) = -a(ilat,1)*(Uloc-5d0*Jh)/2d0*Gamma5 -a(ilat,2)*(Uloc+Jh)/2d0*GammaS
       do io=1,Nso
          do jo=1,Nso
             i = io + (ilat-1)*Nso
             j = jo + (ilat-1)*Nso
             HijMF(i,j) = Htmp(ilat,io,jo)
          enddo
       enddo
    enddo
  end function mf_Hij_correction



  function ts_model(link,Nso) result(Hts)
    integer                       :: link
    integer                       :: Nso
    complex(8),dimension(Nso,Nso) :: Hts
    select case(link)
    case (0) !LOCAL PART
       Hts =  Mh*Gamma5 
    case (1) !RIGHT HOPPING
       Hts = -0.5d0*Gamma5 + xi*0.5d0*lambda*GammaX
    case (2) !UP HOPPING
       Hts = -0.5d0*Gamma5 + xi*0.5d0*lambda*GammaY
    case (3) !LEFT HOPPING
       Hts = -0.5d0*Gamma5 - xi*0.5d0*lambda*GammaX
    case (4) !DOWN HOPPING
       Hts = -0.5d0*Gamma5 - xi*0.5d0*lambda*GammaY
    case default 
       stop "ts_model ERROR: link index in {0..6}"
    end select
  end function ts_model


  function get_Hmf(a,Hij) result(Hmat)
    real(8),dimension(Nlat,2)                :: a
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Hij
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Hmat
    Hmat = Hij + mf_Hij_correction(a)
  end function get_Hmf



  subroutine symmetrize_params(a)
    real(8),dimension(Nlat,Nso) :: a
    integer :: ilat
    do ilat=2,Nlat
       a(ilat,:) = a(1,:)
    enddo
  end subroutine symmetrize_params



end program







! !Build up the real-space Hamiltonian thru FT:"
! allocate(Hij(Nlat*Nso,Nlat*Nso,1)) ; Hij = zero
! call start_timer
! do ilat=1,Nlat
!    vecRi = Rgrid(ilat,:)
!    do jlat=1,Nlat
!       vecRj = Rgrid(jlat,:)
!       !
!       Htmp = zero
!       do ik=1,Nktot
!          vecK = Kgrid(ik,:)
!          arg=dot_product(vecK,vecRj-vecRi)
!          Htmp(:,:)= Htmp(:,:) + exp(xi*arg)*hk_model(vecK,Nso)/Nktot
!       enddo
!       !
!       do io=1,Nso
!          i = io + (ilat-1)*Nso
!          do jo=1,Nso
!             j = jo + (jlat-1)*Nso
!             !
!             Hij(i,j,1) = Htmp(io,jo)
!             !
!          enddo
!       enddo
!       !
!    enddo
!    call eta(ilat,Nlat)
! enddo
! where(abs(Hij)<1.d-6)Hij=zero
! call stop_timer


! print*,Nlat*Nso
! allocate(Evals(Nlat*Nso))
! allocate(rhoH(Nlat*Nso,Nlat*Nso))
! call eigh(Hij(:,:,1),Evals)
! rhoDiag = fermi(Evals,beta)
! rhoH    = matmul(Hij(:,:,1) , matmul(diag(rhoDiag), conjg(transpose(Hij(:,:,1)))) )

! ilat=1
! do io=1,Nso
!    dens(io) = dreal(rhoH(io+(ilat-1)*Nso,io+(ilat-1)*Nso))
! enddo
! write(*,"(A,10F14.9)")"Occupations =",(dens(io),io=1,Nso),sum(dens)
! open(10,file="Robservables.nint")
! write(10,"(10F20.12)")(dens(io),io=1,Nso),sum(dens)
! close(10)


