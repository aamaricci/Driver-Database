program mf_bhz_1d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none


  integer                                     :: Nparams=2
  integer,parameter                           :: Norb=2
  integer,parameter                           :: Nspin=2
  integer,parameter                           :: Nso=Norb*Nspin,L=2048
  integer,parameter                           :: Nkpath=500
  integer                                     :: Nk
  integer                                     :: Nlat,Nx
  integer                                     :: i,j,k,ik,iorb,jorb,ispin,io,jo
  integer                                     :: ilat,jlat
  integer                                     :: ix
  real(8)                                     :: kx
  real(8),dimension(:,:),allocatable          :: kgrid,kpath,ktrims,Rgrid
  integer,dimension(:,:),allocatable          :: Links
  complex(8),dimension(:,:,:),allocatable     :: Hk
  complex(8),dimension(:,:),allocatable       :: Hij,mfHij_glob
  complex(8),dimension(:,:,:,:),allocatable   :: Hlat
  integer 				      :: Iter,MaxIter,Nsuccess=2
  real(8)                                     :: chern,z2,Uloc,Jratio,Jh,Sz,Tz,Rz,Ntot
  real(8)                                     :: mh,lambda
  real(8)                                     :: xmu,beta,eps
  real(8)                                     :: wmix,it_error,sb_field
  complex(8)                                  :: Hloc(Nso,Nso),mfHk_glob(Nso,Nso)
  complex(8),dimension(:,:,:,:,:),allocatable :: Gmats,Greal
  character(len=20)                           :: Finput
  logical                                     :: iexist,converged,withgf
  complex(8),dimension(Nso,Nso)               :: Gamma1,Gamma2,Gamma5,GammaS
  complex(8),dimension(Nso,Nso)               :: GammaN,GammaTz,GammaSz,GammaRz
  complex(8),dimension(Nso,Nso)               :: GammaE0,GammaEx,GammaEy,GammaEz
  real(8),dimension(:),allocatable            :: params,params_prev

  call parse_cmd_variable(Finput,"FINPUT",default="input.conf")
  call parse_input_variable(Nx,"Nx",Finput,default=100)
  call parse_input_variable(mh,"MH","input.conf",default=0d0)
  call parse_input_variable(lambda,"LAMBDA","input.conf",default=0.3d0)
  call parse_input_variable(Uloc,"ULOC",Finput,default=1d0)
  call parse_input_variable(Jratio,"JRATIO",Finput,default=0.25d0)
  call parse_input_variable(xmu,"XMU",Finput,default=0.d0)
  call parse_input_variable(eps,"EPS",Finput,default=4.d-2)
  call parse_input_variable(beta,"BETA",Finput,default=1000.d0)
  call parse_input_variable(wmix,"WMIX",Finput,default=0.5d0)
  call parse_input_variable(sb_field,"SB_FIELD",Finput,default=0.1d0)
  call parse_input_variable(it_error,"IT_ERROR",Finput,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",Finput,default=100)
  call parse_input_variable(withgf,"WITHGF",Finput,default=.false.)
  !
  call print_input(trim(Finput))
  call save_input_file(trim(Finput))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-6d0,"wini")
  call add_ctrl_var(6d0,"wfin")
  call add_ctrl_var(eps,"eps")

  Jh = Jratio*Uloc

  Nlat = Nx

  write(*,*)"Solving BHZ with Hartree terms: [Tz,Sz]"

  !
  allocate( params(Nparams),params_prev(Nparams) )
  !
  !SETUP THE GAMMA MATRICES:
  gamma1=kron( pauli_sigma_z, pauli_tau_x)
  gamma2=kron( pauli_sigma_0,-pauli_tau_y)
  gamma5=kron( pauli_sigma_0, pauli_tau_z)
  gammaS=kron( pauli_sigma_z, pauli_tau_0)
  !
  gammaN=kron( pauli_sigma_0, pauli_tau_0 )
  gammaTz=kron( pauli_sigma_0, pauli_tau_z )
  gammaSz=kron( pauli_sigma_z, pauli_tau_0 )
  gammaRz=kron( pauli_sigma_z, pauli_tau_z )
  !
  gammaE0=kron( pauli_sigma_0, pauli_tau_x )
  gammaEx=kron( pauli_sigma_x, pauli_tau_x )
  gammaEy=kron( pauli_sigma_y, pauli_tau_x )
  gammaEz=kron( pauli_sigma_z, pauli_tau_x )


  !Setup the lattice basis:
  call TB_set_bk([pi2,0d0])
  call TB_set_ei([1d0,0d0])
  !
  !BUILD H(k) 
  allocate(Hk(Nso,Nso,Nlat))
  Hk = zero
  mfHk_glob = zero
  call TB_build_model(Hk,hk_model,Nso,[Nlat])
  !
  !
  params   = sb_field      ![Tz,Sz]
  inquire(file="params.restart",exist=iexist)
  if(iexist)then
     call read_array("params.restart",params)
     params(2)=params(2)+sb_field
  endif
  call save_array("params.init",params)
  !
  !
  open(100,file="tz_szVSu_hk.dat")
  converged=.false. ; iter=0
  do while(.not.converged.AND.iter<maxiter)
     iter=iter+1
     call start_loop(iter,maxiter,"MF-loop")
     !
     call solve_MF_Hk(iter,params)
     !
     if(iter>1)params = wmix*params + (1d0-wmix)*params_prev
     params_prev = params
     !
     converged = check_convergence_local(params,it_error,nsuccess,maxiter)
     !
     call end_loop
  end do
  call save_array("params.restart",params)
  close(100)
  !
  !Update Global Mean-Field Hamiltonian correction:
  mfHk_glob =  mf_hk_correction(params)
  call TB_write_hloc(mfHk_glob)

contains

  subroutine solve_MF_Hk(iter,a)
    real(8),dimension(:),intent(inout) :: a
    complex(8),dimension(Nso,Nso)      :: H,Hmf
    real(8),dimension(Nso)             :: Ek,rhoDiag
    complex(8),dimension(Nso,Nso)      :: rhoHk,rhoH
    real(8)                            :: N,Tz,Sz,Rz
    real(8)                            :: E0,Ez,Ex,Ey
    integer                            :: ik,iter
    !
    rewind(100)
    !
    Tz = 0d0
    Sz = 0d0
    !
    Hmf=mf_hk_correction(a)
    !
    rhoH = zero
    do ik=1,Nx
       H   = Hk(:,:,ik) + Hmf
       !
       call eigh(H,Ek)       !diag Hk --> Ek
       !
       rhoDiag = fermi(Ek,beta)
       rhoHk   = matmul( H, matmul(diag(rhoDiag),conjg(transpose(H))) )
       rhoH    = rhoH + rhoHk/Nx
       !
    enddo
    !
    Tz = Tz + sum( GammaTz*rhoH )
    Sz = Sz + sum( GammaSz*rhoH )
    !
    write(*,*)iter,Tz,Sz
    write(100,"(3F21.12)")uloc,Tz,Sz
    a = [Tz,Sz]
    return
  end subroutine solve_MF_PBC



  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: ek
    real(8)                   :: kx
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    !
    kx=kpoint(1)
    ek = -1d0*cos(kx)
    Hk = (Mh+ek)*Gamma5 + lambda*sin(kx)*Gamma1
    !
    Hk  = Hk + mfHk_glob
  end function hk_model



  function mf_Hk_correction(a) result(HkMF)
    real(8),dimension(:)          :: a
    complex(8),dimension(Nso,Nso) :: HkMF
    HkMF = -a(1)*(Uloc-5d0*Jh)/4d0*Gamma5 &
         -a(2)*(Uloc+Jh)/4d0*GammaS    
  end function mf_Hk_correction



end program mf_bhz_1d


