!BHZ Model with e0=2t=1 
program ss_bhz2d
  USE SLAVE_SPINS
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  implicit none
  complex(8),allocatable,dimension(:,:,:) :: Hk
  complex(8),allocatable,dimension(:,:)   :: Hloc
  real(8)                                 :: mh,rh,dh,lambda
  integer                                 :: Nkx,Nktot,Nlso,ilat,Npts,Nkpath,ik
  real(8),dimension(:,:),allocatable      :: kpath
  complex(8),dimension(4,4)               :: Gamma1,Gamma2,Gamma5


  call parse_input_variable(lambda,"LAMBDA","inputSS.conf",default=0.3d0)
  call parse_input_variable(mh,"MH","inputSS.conf",default=0.d0)
  call parse_input_variable(rh,"RH","inputSS.conf",default=0.d0)
  call parse_input_variable(dh,"DH","inputSS.conf",default=0.d0)
  call parse_input_variable(Nkx,"Nkx","inputSS.conf",default=20)
  call parse_input_variable(nkpath,"NKPATH","inputSS.conf",default=500)
  call ss_read_input('inputSS.conf')

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Nlat,"nlat")
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  write(*,*)"Unit energy 2t=1"
  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nlso=Nspin*Norb

  Nktot = Nkx**2

  gamma1=kron_pauli( pauli_sigma_z, pauli_tau_x)
  gamma2=kron_pauli( pauli_sigma_0,-pauli_tau_y)
  gamma5=kron_pauli( pauli_sigma_0, pauli_tau_z) !GammaMH
  !

  allocate(Hk(Nlso,Nlso,Nktot),Hloc(Nlso,Nlso))
  call TB_set_bk([pi2,0d0],[0d0,pi2])
  call TB_build_model(Hk,hk_model,Nlso,[Nkx,Nkx])
  Hloc = sum(Hk,dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=0d0
  call TB_write_hloc(Hloc)

  !solve along a path in the 2d BZ.
  Npts = 4
  allocate(kpath(Npts,2))
  kpath(1,:)=[0d0,0d0]
  kpath(2,:)=[1d0,1d0]
  kpath(3,:)=[1d0,0d0]
  kpath(4,:)=[0d0,0d0]
  kpath = kpath*pi

  !Solve for the renormalized bands:
  call TB_Solve_model(hk_model,Nlso,kpath,Nkpath,&
       colors_name=[red,pink,yellow,blue],&
       points_name=[character(len=40) ::'G', 'M', 'X', 'G'],&
       file="Bands_BHZ",iproject=.false.)



  call start_timer
  call ss_solve(Hk)
  call stop_timer("SS SOLUTION")



  !Solve for the renormalized bands:
  call TB_Solve_model(ss_hk_model,Nlso,kpath,Nkpath,&
       colors_name=[red,pink,yellow,blue],&
       points_name=[character(len=40) ::'G', 'M', 'X', 'G'],&
       file="zBands_BHZ",iproject=.false.)


contains

  function hk_model(kvec,N) result(hk)
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: kx,ky,ek
    integer                   :: N,ii
    if(N/=Nlso)stop "hk_modle error: N != Nspin*Norb == 4"
    kx = kvec(1)
    ky = kvec(2)
    ek = -1d0*(cos(kx)+cos(ky)) !
    Hk = zero
    Hk = (Mh+ek)*Gamma5 + lambda*sin(kx)*Gamma1 + lambda*sin(ky)*Gamma2
    Hk(1,4) = -dh ; Hk(4,1)=-dh
    Hk(2,3) =  dh ; Hk(3,2)= dh
    Hk(1,3) = xi*rh*(sin(kx)-xi*sin(ky))
    Hk(3,1) =-xi*rh*(sin(kx)+xi*sin(ky))
  end function hk_model


  function ss_Hk_model(kvec,N) result(Hk)
    real(8),dimension(:)      :: kvec
    integer                   :: N
    complex(8),dimension(N,N) :: Hk
    !< Build H(k) from W90:
    Hk = hk_model(kvec,N)
    !< Build effective fermionic H*(k) 
    call ss_get_ssHk(Hk,Hloc)
  end function ss_Hk_model


end program ss_bhz2d
