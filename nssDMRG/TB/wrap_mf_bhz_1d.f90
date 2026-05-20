subroutine mf_bhz_1d(uloc,mh,res)
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  real(8),intent(in)                      :: uloc,mh
  integer,intent(out)                     :: res
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
  call parse_input_variable(Jratio,"JRATIO",Finput,default=0.25d0)
  call parse_input_variable(xmu,"XMU",Finput,default=0.d0)
  call parse_input_variable(eps,"EPS",Finput,default=4.d-2)
  call parse_input_variable(beta,"BETA",Finput,default=1000.d0)
  call parse_input_variable(wmix,"WMIX",Finput,default=0.5d0)
  call parse_input_variable(sb_field,"SB_FIELD",Finput,default=0.1d0)
  call parse_input_variable(it_error,"IT_ERROR",Finput,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",Finput,default=100)
  call parse_input_variable(pbc,"pbc","input.conf",default=.true.)
  !
  !
  Jh = Jratio*Uloc

  Nlat = Nx
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


  params_prev= 0d0
  params     = sb_field      ![Tz,Sz,S0,Ex,Ey,Ez]
  !
  write(*,"(2F12.7)",advance="no")Uloc,Mh
  converged=.false. ; iter=0
  do while(.not.converged.AND.iter<maxiter)
     iter=iter+1
     call solve_MF_Hij(iter,params)
     if(iter>1)params = wmix*params + (1d0-wmix)*params_prev
     params_prev = params
     converged = check_error(params,it_error,1,maxiter,error)
  end do
  !
  !
  threshold=1d-2
  res=1
  Tz = params(1)
  Sz = params(2)
  if(abs(Sz)>threshold)res=2
  write(*,"(I4,2F12.7)",advance="no")res,Tz,Sz
  inquire(file="params.run",exist=iexist)
  if(iexist)then
     open(100,file="params.run",status="old", position="append", action="write")
  else
     open(100,file="params.run",status="new", action="write")
  endif
  write(100,*)Uloc,Mh,Tz,Sz,res,jhratio
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

function check_error(Xnew,eps,N1,N2,error) result(convergence)
    real(8),intent(in)            :: Xnew(:)
    real(8)                       :: eps
    integer                       :: N1,N2
    integer                       :: Msize1
    logical                       :: convergence  
    real(8)                       :: err,error
    real(8),dimension(size(Xnew)) :: Verror
    real(8),save,allocatable      :: Xold(:,:)
    integer,save                  :: success=0,check=1
    Msize1=size(Xnew)
    if(.not.allocated(Xold))then
       allocate(Xold(1,Msize1))
       Xold=0.d0
    endif
    Verror=abs(Xnew-Xold(1,:))
    if(check==1)Verror=1d0
    err=sum(Verror)/dble(size(Verror))
    Xold(1,:)=Xnew
    if(err < eps)then
       success=success+1
    else
       success=0
    endif
    convergence=.false.
    if(success > N1)convergence=.true.
    ! if(convergence)then
    !    write(*,"(A,ES15.7,I8)")bold_green("    error="),err
    ! else
    !    if(err < eps)then
    !       write(*,"(A,ES15.7,I8)")bold_yellow("    error="),err
    !    else
    !       write(*,"(A,ES15.7,I8)")bold_red("    error="),err
    !    endif
    ! endif
    error=err
    ! if(check>=N2)convergence=.true.
    check=check+1

  end function check_error

end subroutine mf_bhz_1d


