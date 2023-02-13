!3bands + SOC Model with e0=2t=1, e1=2t_1=2.a.t (a<1)
! bands order: XY(1), XZ(2), YZ(3)
program ss_hmsoc
  USE SLAVE_SPINS
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  implicit none
  integer                                 :: Nlso
  complex(8),allocatable,dimension(:,:,:) :: Hk
  complex(8),allocatable,dimension(:,:)   :: Hloc
  real(8)                                 :: ts(3),Mh(3),lambda,alfa
  integer                                 :: Nkx,Nktot,Npts,Nkpath,ik
  real(8),dimension(:,:),allocatable      :: kpath

  call parse_input_variable(ts,"TS","inputSOC.conf",default=[0.5d0,0.5d0,0.5d0])
  call parse_input_variable(lambda,"LAMBDA","inputSOC.conf",default=0.d0)
  call parse_input_variable(mh,"MH","inputSOC.conf",default=[0d0,0.1d0,0.1d0])
  call parse_input_variable(Nkx,"Nkx","inputSOC.conf",default=20)
  call parse_input_variable(nkpath,"NKPATH","inputSOC.conf",default=500)
  call ss_read_input('inputSOC.conf')

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Nlat,"nlat")
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(alfa>1d0.OR.alfa<0d0)stop "ERROR: alfa must be [0:1] or die"
  if(Nspin/=2.OR.Norb/=3)stop "ERROR: Norb==3, Nspin==2 or die"
  Nlso  = Nlat*Nspin*Norb
  Nktot = Nkx**2


  allocate(Hk(Nlso,Nlso,Nktot),Hloc(Nlso,Nlso))
  call TB_set_bk([pi2,0d0],[0d0,pi2])
  call TB_build_model(Hk,hk_model,Nlso,[Nkx,Nkx])
  Hloc = sum(Hk,dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=0d0
  call TB_write_hloc(Hloc)

  !solve along a path in the 2d BZ.
  Npts = 7
  allocate(kpath(Npts,2))
  kpath(1,:)=[0d0,0d0]
  kpath(2,:)=[0.5d0,0d0]
  kpath(3,:)=[0.5d0,0.5d0]
  kpath(4,:)=[0d0,0d0]
  kpath(5,:)=[0d0,0.5d0]
  kpath(6,:)=[0.5d0,0.5d0]
  kpath(7,:)=[0d0,0d0]
  kpath = kpath*pi2  
  !Solve for the renormalized bands:
  call TB_Solve_model(hk_model,Nlso,kpath,Nkpath,&
       colors_name=[red,blue,green,darkred,darkblue,darkgreen],&
       points_name=[character(len=40) :: '{/Symbol} G', 'X', 'M', '{/Symbol} G','Y', 'M', '{/Symbol} G'],&
       file="Bands_SOC",iproject=.false.)



  call start_timer
  call ss_solve(Hk)
  call stop_timer("SS SOLUTION")



  !Solve for the renormalized bands:
  call TB_Solve_model(ss_hk_model,Nlso,kpath,Nkpath,&
       colors_name=[red,blue,green,darkred,darkblue,darkgreen],&
       points_name=[character(len=40) :: '{/Symbol} G', 'X', 'M', '{/Symbol} G','Y', 'M', '{/Symbol} G'],&
       file="zBands_SOC",iproject=.false.)


contains

  function hk_model(kvec,N) result(hk)
    real(8),dimension(:)      :: kvec
    integer                   :: N
    complex(8),dimension(N,N) :: hk
    real(8)                   :: kx,ky
    real(8)                   :: cx,cy,e,ek(3)
    if(N/=Nlso)stop "hk_modle error: N != Nspin*Norb == 6"
    kx = kvec(1)
    ky = kvec(2)
    cx = cos(kx)
    cy = cos(ky)
    ek = -2d0*ts*(cx+cy)
    Hk = zero
    Hk = diag([ek,ek]) + diag([Mh,Mh]) + lambda*H_ls()
  end function hk_model


  function H_LS() result(hls)
    complex(8),dimension(Nlso,Nlso) :: hls
    integer                         :: i,j
    hls=zero
    hls(1,2) = xi
    hls(1,6) = -1d0
    hls(2,6) = xi
    hls(3,4) = 1d0
    hls(3,5) = -xi
    hls(4,5) = -xi
    do i=1,Nlso
       do j=1,i
          hls(i,j)=conjg(hls(j,i))
       enddo
    enddo
  end function H_LS



  function ss_Hk_model(kvec,N) result(Hk)
    real(8),dimension(:)      :: kvec
    integer                   :: N
    complex(8),dimension(N,N) :: Hk
    Hk = hk_model(kvec,N)
    !< Build effective fermionic H*(k) 
    call ss_get_ssHk(Hk,Hloc)
  end function ss_Hk_model


end program ss_hmsoc
