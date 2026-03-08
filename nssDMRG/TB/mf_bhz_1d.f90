program mf_bhz_1d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none


  integer                                   :: Nparams=2
  integer,parameter                         :: Norb=2
  integer,parameter                         :: Nspin=2
  integer,parameter                         :: Nso=Norb*Nspin
  integer                                   :: Nlat,Nx
  integer                                   :: i,j,io,jo
  integer                                   :: ilat,jlat
  integer,dimension(:,:),allocatable        :: Links
  complex(8),dimension(:,:),allocatable     :: Hij,mfHij_glob
  complex(8),dimension(:,:,:,:),allocatable :: Hlat
  integer                                   :: Iter,MaxIter,Nsuccess=2
  real(8)                                   :: Uloc,Jratio,Jh
  real(8)                                   :: mh,lambda
  real(8)                                   :: xmu,beta,eps
  real(8)                                   :: wmix,it_error,sb_field
  character(len=20)                         :: Finput
  logical                                   :: iexist,converged,withgf
  complex(8),dimension(Nso,Nso)             :: Gamma1,Gamma2,Gamma5,GammaS
  complex(8),dimension(Nso,Nso)             :: GammaN,GammaTz,GammaSz,GammaRz
  complex(8),dimension(Nso,Nso)             :: GammaE0,GammaEx,GammaEy,GammaEz
  logical                                   :: pbc
  real(8),dimension(:),allocatable          :: params,params_prev

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
  call parse_input_variable(pbc,"pbc","input.conf",default=.true.)
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
  Nparams=Nlat*Nparams
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
  allocate(Links(2,1))          !Links: right,left
  Links(1,:) = [1]
  Links(2,:) =-[1]
  !
  !BUILD H(k) or H(i,j)
  allocate(Hlat(Nso,Nso,Nlat,Nlat))
  allocate(Hij(Nlat*Nso,Nlat*Nso))
  allocate(mfHij_glob(Nlat*Nso,Nlat*Nso))
  mfHij_glob= zero
  Hlat      = zero
  Hij       = zero
  call TB_build_model(Hlat,ts_model,Nso,[Nlat],Links,pbc=pbc)
  do concurrent(ilat=1:Nlat,jlat=1:Nlat,io=1:Nso,jo=1:Nso)
     i = io + (ilat-1)*Nso        
     j = jo + (jlat-1)*Nso
     Hij(i,j) = Hlat(io,jo,ilat,jlat)
  enddo
  !
  !
  !
  params = 0d0
  ! params(1:Nlat)   = sb_field
  ! do ilat=1,Nlat
  !    params(Nlat+ilat)= (-1d0)**(ilat+1)*sb_field
  ! enddo
  inquire(file="params.restart",exist=iexist)
  if(iexist)then
     call read_array("params.restart",params)
  endif
  call save_array("params.init",params)
  !
  !
  if(pbc)then
     open(100,file="tz_szVSu_PBC.dat")
  else
     open(100,file="tz_szVSu_OBC.dat")
  endif
  converged=.false. ; iter=0
  do while(.not.converged.AND.iter<maxiter)
     iter=iter+1
     call start_loop(iter,maxiter,"MF-loop")
     !
     call solve_MF_Hij(iter,params)
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
  mfHij_glob= mf_hij_correction(params)
  if(pbc)then
     open(100,file="tz_szVSiVSu_PBC.dat")
  else
     open(100,file="tz_szVSiVSu_OBC.dat")
  endif
  do ilat=1,Nlat
     write(100,*)ilat,params(ilat),params(Nlat+ilat)
  enddo
  close(100)


contains



  subroutine solve_MF_Hij(iter,a)
    real(8),dimension(:),intent(inout)      :: a
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: H
    real(8),dimension(Nlat*Nso)             :: E,rhoDiag
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: rhoH
    real(8),dimension(Nlat,Nspin,Norb)      :: dens
    real(8),dimension(Nlat)                 :: Tz,Sz,n1,n2,nup,ndw
    integer                                 :: iter,ilat,ispin,iorb
    !
    rewind(100)
    !
    Tz = 0d0
    Sz = 0d0
    !
    H    = Hij + mf_hij_correction(a)
    !
    call eigh(H,E)
    !
    rhoDiag = fermi(E,beta)
    rhoH    = matmul(H , matmul(diag(one*rhoDiag), conjg(transpose(H))) ) 
    !
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,iorb=1:Norb)
       i = iorb + (ispin-1)*Norb + + (ilat-1)*Nso
       dens(ilat,ispin,iorb) = dreal(rhoH(i,i))
    enddo
    !
    N1 = dens(:,1,1)+ dens(:,2,1)
    N2 = dens(:,1,2)+ dens(:,2,2)
    Nup= dens(:,1,1)+ dens(:,1,2)
    Ndw= dens(:,2,1)+ dens(:,2,2)
    Tz = N1 - N2
    Sz = Nup- Ndw
    !
    write(*,*)iter,sum(Tz)/Nlat,sum(abs(Sz))/Nlat
    write(100,"(3F21.12)")uloc,sum(Tz)/Nlat,sum(abs(Sz))/Nlat
    a = [Tz,Sz]
    return
  end subroutine solve_MF_Hij




  function ts_model(link,N) result(Hts)
    integer                   :: link
    integer                   :: N
    complex(8),dimension(N,N) :: Hts
    select case(link)
    case (0) !LOCAL PART
       Hts =  Mh*Gamma5 !+ mfHk_glob 
    case (1) !RIGHT HOPPING
       Hts = -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma1
    case (2) !LEFT HOPPING
       Hts = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma1
    case default 
       stop "ts_model ERROR: link != [0:2]"
    end select
  end function ts_model



  function mf_Hk_correction(a) result(HkMF)
    real(8),dimension(:)          :: a
    complex(8),dimension(Nso,Nso) :: HkMF
    HkMF = -a(1)*(Uloc-5d0*Jh)/4d0*Gamma5 &
         -a(2)*(Uloc+Jh)/4d0*GammaS    
  end function mf_Hk_correction



  function mf_Hij_correction(a) result(HijMF)
    real(8),dimension(:)                    :: a
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: HijMF
    real(8)                                 :: a_(2)
    complex(8),dimension(Nso,Nso)           :: H_
    integer                                 :: ilat,io,jo
    !
    HijMF=zero
    do ilat=1,Nlat
       a_ = [a(ilat),a(Nlat+ilat)]
       H_ = mf_Hk_correction(a_)
       do concurrent(io=1:Nso,jo=1:Nso)
          i = io + (ilat-1)*Nso        
          j = jo + (ilat-1)*Nso
          HijMF(i,j) = H_(io,jo)
       enddo
    enddo
    !
  end function mf_Hij_correction



end program mf_bhz_1d


