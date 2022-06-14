program bhz_2d
  !LIBRARIES:
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  !
  !SET THE DIMENSION OF THE PROBLEM
  integer,parameter                           :: Norb=2
  integer,parameter                           :: Nspin=2
  integer,parameter                           :: Nso=Nspin*Norb
  !
  !INPUT VARIABLES
  integer                                     :: Nkx,Nky,Nktot
  integer                                     :: L,Lf,Lb
  integer                                     :: Iter
  integer                                     :: MaxIter
  integer                                     :: Nsuccess=0
  integer                                     :: OrderInt
  real(8)                                     :: gt
  real(8)                                     :: gn
  real(8)                                     :: mh
  real(8)                                     :: lambda
  real(8)                                     :: xmu
  real(8)                                     :: beta
  real(8)                                     :: eps
  real(8)                                     :: wmix
  real(8)                                     :: it_error
  real(8)                                     :: tz0
  real(8)                                     :: dtz0
  real(8)                                     :: dplus0,dminus0
  logical                                     :: withgf,check_nel
  !
  !Gamma matrices:
  complex(8),dimension(Nso,Nso)               :: Gamma1,Gamma2,Gamma5,GammaN
  !
  !GLOBALLY SHARED VARIABLES:
  real(8)                                     :: Tz !the orbital polarization
  real(8)                                     :: Ne !the total density
  real(8),dimension(:,:),allocatable          :: kgrid !K vector grid
  !
  !WORK ARRAYS:
  real(8),dimension(:),allocatable            :: p_work
  real(8)                                     :: qvec_work(2)
  integer                                     :: m_work
  !
  integer                                     :: Nparams
  integer                                     :: i,j,k,ik,m
  integer                                     :: info,unit
  real(8)                                     :: x(1),dx(1)
  logical                                     :: iexist
  logical                                     :: converged
  complex(8),dimension(:,:,:),allocatable     :: Hk
  complex(8),dimension(:,:,:,:,:),allocatable :: Smats
  complex(8),dimension(:,:,:,:,:),allocatable :: Gmats
  character(len=20)                           :: Finput
  real(8),dimension(:),allocatable            :: params
  real(8),dimension(:),allocatable            :: params_prev
  real(8),dimension(2,2)                      :: Chi_qv

  call parse_cmd_variable(Finput,"FINPUT",default="input.conf")
  call parse_input_variable(gt,"GT",Finput,default=1d0)
  call parse_input_variable(gn,"GN",Finput,default=1d0)
  call parse_input_variable(nkx,"NKX",Finput,default=10)
  call parse_input_variable(Lf,"LF",Finput,default=256,comment="# of fermionic Mats frequencies, L=Lf+Lb")
  call parse_input_variable(Lb,"LB",Finput,default=64,comment="# of bosonix Mats frequencies, L=Lf+Lb")
  call parse_input_variable(Mh,"MH",Finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",Finput,default=0.3d0)
  call parse_input_variable(xmu,"XMU",Finput,default=0.d0)  
  call parse_input_variable(beta,"BETA",Finput,default=1000.d0)
  call parse_input_variable(tz0,"tz0",Finput,default=-0.1d0,comment="Guess for MF search of Tz (tz0<0)")
  call parse_input_variable(dtz0,"dtz0",Finput,default=0.1d0,comment="Guess for dTz fluctuations (dtz0>0)")
  call parse_input_variable(dplus0,"dplus0",Finput,default=0.1d0,comment="Guess for d+0 fluctuations (d+0>0)")
  call parse_input_variable(dminus0,"dminus0",Finput,default=0.1d0,comment="Guess for d-0 fluctuations (d-0>0)")
  call parse_input_variable(it_error,"IT_ERROR",Finput,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",Finput,default=100)    
  call parse_input_variable(eps,"EPS",Finput,default=4.d-2)
  call parse_input_variable(wmix,"WMIX",Finput,default=1d0)
  call parse_input_variable(withgf,"WITHGF",Finput,default=.false.)
  call parse_input_variable(check_nel,"check_nel",Finput,default=.false.)
  call parse_input_variable(OrderInt,"ORDERINT",Finput,default=0)
  call print_input(trim(Finput))
  call save_input_file(trim(Finput))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-10d0,"wini")
  call add_ctrl_var(10d0,"wfin")
  call add_ctrl_var(eps,"eps")

  if(gn==0d0.OR.gt==0d0)stop "ERROR: gt AND gn should be non zero (check in progress)"

  gamma1 = kron_pauli( pauli_sigma_z, pauli_tau_x)
  gamma2 = kron_pauli( pauli_sigma_0,-pauli_tau_y)
  gamma5 = kron_pauli( pauli_sigma_0, pauli_tau_z)
  gammaN = kron_pauli( pauli_sigma_0, pauli_tau_0)
  Nky    = Nkx
  Nktot  = Nkx*Nky
  print*,Nktot
  !
  L      = Lf+Lb
  write(*,*)"Using L freq.=",L


  !>SOLVE MF PROBLEM 1st: >>ACTHUNG<< This solution does not use BZ basis defined later!!
  call start_timer()
  x(1)=-abs(tz0)
  dx(1)=0.1d0
  call fmin(bhz_f,x,lambda=dx)
  Tz=x(1)
  open(free_unit(unit),file="mf_TzVSg.dat")
  write(unit,*)gt,gn,Tz
  close(unit)
  write(*,*) "Tz=",Tz
  call stop_timer(" Mean-Field")
  call sleep(1)



  !> SOLVE FLUCTUATIONS:
  !Setup the k-space lattice basis:
  call TB_set_bk([pi2,0d0],[0d0,pi2])

  allocate(kgrid(Nktot,2))      !Nktot=# tot kpoints, 2= 2D
  call TB_build_kgrid([Nkx,Nky],kgrid)

  !+ ReSigma(iw_n)[L] + ImSigma(iw_n)[L] + Tz[1] + <|d+|**2>[1] + <|d-|**2> [1]
  Nparams = 3 + 2*L
  allocate( params(Nparams), params_prev(Nparams), p_work(Nparams))

  !Start from MF solution
  params = [dble(zeros(L)),dble(zeros(L)),Tz,abs(dplus0),abs(dminus0)]
  inquire(file="params.restart",exist=iexist)
  if(iexist)call read_array("params.restart",params)
  call save_array("params.init",params)

  open(free_unit(unit),file="chi0.dat ")
  do m=1,Lb
     p_work   = params
     qvec_work=[0d0,0d0]
     call Intk_SumMats_Chi_qv(m,Chi_qv)
     write(unit,*) 2d0*pi*m/beta, Chi_qv(1,1),Chi_qv(2,2),Chi_qv(1,2)
  end do
  close(unit)     


  converged=.false. ; iter=0
  do while(.not.converged.AND.iter<maxiter)
     iter=iter+1
     call start_loop(iter,maxiter,"SC-loop")
     !
     !>SOLVE 4 EQUATIONS For Delta,<|dDelta|**2>,Sigma
     call solve_eqs(params)
     if(iter>1)params = wmix*params + (1d0-wmix)*params_prev
     params_prev = params
     !
     converged = check_convergence_local(params,it_error,nsuccess,maxiter) 
     !
     call end_loop
  end do
  call save_array("params.restart",params)	!ok forse va salvato anche dSigma, ma only last step(?)
  !
  open(free_unit(unit),file="tz_dtzVSg.dat")
  write(unit,*)gt,gn,params(2*L+1),params(2*L+2),params(2*L+3)
  close(unit)
  write(*,*) "Tz,d+,d-=",params(2*L+1),params(2*L+2),params(2*L+3)


  allocate(Smats(Nspin,Nspin,Norb,Norb,L))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,L))
  allocate(Hk(Nso,Nso,Nktot));Hk=zero
  !
  call TB_build_model(Hk,hk_bhz,Nso,[Nkx,Nky])
  call build_self_energy(params,Smats)
  call dmft_gloc_matsubara(Hk,Gmats,Smats)
  call dmft_print_gf_matsubara(Smats,"Smats",iprint=1)
  call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)


  open(free_unit(unit),file="chi.dat ")
  do m=1,Lb
     p_work   = params
     qvec_work=[0d0,0d0]
     call Intk_SumMats_Chi_qv(m,Chi_qv)
     write(unit,*) 2d0*pi*m/beta, Chi_qv(1,1),Chi_qv(2,2),Chi_qv(1,2)
  end do
  close(unit)     


contains



  !For MF calculations:
  !#################################
  function bhz_f(a) result(f)
    real(8),dimension(:) :: a
    real(8)              :: f
    real(8)              :: integral
    Tz = a(1)
    call gauss_quad(Fmf_bhz,[0d0,0d0],[pi2,pi2],integral)
    f = gt*(Tz**2) - 2d0*integral/pi2/pi2
  end function bhz_f

  function Fmf_bhz(kvec) result(Fk)
    real(8),dimension(:) :: kvec
    real(8)              :: Fk
    real(8)              :: ek,x2,y2,kx,ky
    kx  = kvec(1)
    ky  = kvec(2)
    ek  = -1d0*(cos(kx)+cos(ky))
    x2  =  lambda*sin(kx) ;x2=x2**2
    y2  =  lambda*sin(ky) ;y2=y2**2
    Fk  = sqrt( (mh - gt*Tz + ek)**2 + (x2+y2) )
  end function fmf_bhz
  !#################################








  subroutine solve_eqs(p)
    real(8),dimension(:),intent(inout) :: p ![3+2L]
    real(8),dimension(Nparams)         :: integral,ints
    real(8)                            :: Tz,dTz
    real(8)                            :: TzTmp,dPlus,dMinus
    real(8),dimension(L)               :: ReSigma,ImSigma
    real(8)                            :: N_el 
    complex(8)                         :: Smats(Nspin,Nspin,Norb,Norb,L)
    complex(8)                         :: Gkmats(Nso,Nso,L)
    real(8)                            :: n_k(Nktot,Nso),kvec(2)
    !
    p_work = p
    !
    print*,"Integrating F_System:",Nktot
    select case(OrderInt)
    case default
       integral  = 0d0
       do ik=1,Nktot
          kvec   = Kgrid(ik,:)
          integral = integral + fk_system(kvec,Nparams)/Nktot
       enddo
    case(1)
       integral  = trapz2d_system(Nparams,fk_system,[0d0,pi2],[0d0,pi2],Nkx,Nkx)/pi2/pi2
    end select
    !
    resigma   = integral(1:L)
    ImSigma   = integral(L+1:2*L)
    Tz        = integral(2*L+1)
    dPlus     = integral(2*L+2)
    dMinus    = integral(2*L+3)    
    !
    !Update params:
    p(1:L)     = ReSigma
    p(L+1:2*L) = ImSigma
    p(2*L+1)   = Tz
    p(2*L+2)   = dPlus
    p(2*L+3)   = dMinus
    !
    if(check_nel)then
       call build_self_energy(p,Smats)
       do ik=1,Nktot
          kvec = Kgrid(ik,:)
          do i=1,L              
             Gkmats(:,:,i) = get_simplified_gf(kvec,i,Nso,so2j(Smats(:,:,:,:,i)))          
          end do
          do i=1,Nso
             n_k(ik,i) = fft_get_density(Gkmats(i,i,:),beta)
          end do
       end do
       N_el = sum(n_k)/Nktot
       write(*,*)N_el
       call splot("Sigma_iw_iter"//str(iter,3)//".dat",pi/beta*(2*arange(1,L)-1),dcmplx(ReSigma(:),ImSigma(:)))
    endif
    write(*,*)iter,Tz,dPlus,dMinus,ReSigma(1),ImSigma(1)
    return
  end subroutine solve_eqs




  function fk_system(kvec,dim) result(integral)
    real(8),dimension(:)   :: kvec
    integer                :: dim
    real(8),dimension(dim) :: integral
    !
    integer                :: n,m,ik,iq
    real(8)                :: kx,ky
    real(8),dimension(L)   :: ReSigma,ImSigma
    real(8)                :: ReS(L),ImS(L)
    real(8)                :: Tz,dTz,dPlus,dMinus
    real(8)                :: TzTmp,dTzTmp,dPlusTmp,dMinusTmp
    real(8)                :: z,simZ
    real(8)                :: Ek,simEk
    real(8)                :: Meff
    real(8)                :: xk,yk
    real(8)                :: Den,gplus,gminus
    real(8)                :: ChiTmp(2,2)
    !
    if(dim/=Nparams)stop "f_system ERROR: Dim != Nparams"
    !
    !split params as required:
    ReSigma   = p_work(1:L)
    ImSigma   = p_work(L+1:2*L)
    Tz        = p_work(2*L+1)
    dPlus     = p_work(2*L+2)
    dMinus    = p_work(2*L+3)
    qvec_work = kvec            !q==k
    !
    kx   = kvec(1)
    ky   = kvec(2)
    ek   = -1d0*(cos(kx)+cos(ky))
    xk   = lambda*sin(kx)
    yk   = lambda*sin(ky)
    !
    gplus    =  1d0/gt + 1d0/gn 
    gminus   = -1d0/gt + 1d0/gn
    !
    !Sum over Matsubara frequencies first:
    ReS       = 0d0
    ImS       = 0d0
    TzTmp     = 0d0
    dPlusTmp  = 0d0
    dMinusTmp = 0d0
    do n=1,L
       z      = pi/beta*(2*n-1)
       simZ   = z - ImSigma(n)  !modulo a factor xi
       Meff   = Mh - Tz*gt + ek + ReSigma(n)
       simEk  = Meff**2 + xk**2 + yk**2
       Den    = z**2d0 + simEk
       !
       TzTmp  = TzTmp  + Meff/Den !Sum over n
       ReS(n) = -simZ/Den*dMinus+ Meff/Den*dPlus
       ImS(n) =  simZ/Den*dPlus + Meff/Den*dMinus
    enddo
    !
    !Int(k) Chi({q,m};k) = Int(k)Sum(n) Chi({q,m};{k,n})
    do m=1,Lb
       call Intk_SumMats_Chi_qv(m,ChiTmp)
       Den       = (ChiTmp(1,1) + gminus)**2 - ChiTmp(1,2)**2 - (ChiTmp(2,2) + gplus)**2
       dPlusTmp  = dPlusTmp  - 2d0*(ChiTmp(1,1) + gminus)/Den - (gt-gn)/2d0  
       dMinusTmp = dMinusTmp + 2d0*ChiTmp(1,2)/Den
    enddo
    !
    integral(1:L)     = -ReS                !ReSigma(k)
    integral(L+1:2*L) = -ImS                !ImSigma(k)
    integral(2*L+1)   = -4d0*TzTmp/beta     !Tz(k)
    integral(2*L+2)   =  2d0*dPlusTmp/beta  !d+(q)
    integral(2*L+3)   =  2d0*dMinusTmp/beta !d-(q)
  end function fk_system



  subroutine Intk_SumMats_Chi_qv(m,chi)
    integer              :: m
    real(8)              :: chi(2,2)
    real(8),dimension(2) :: kvec
    real(8)              :: Tz,dPlus,dMinus
    real(8),dimension(L) :: ReSigma,ImSigma
    real(8)              :: wn,wnm
    real(8)              :: zn,zn_plus_m
    real(8)              :: Mk,Mk_plus_q
    real(8)              :: ek,ek_plus_q
    real(8)              :: xk,xk_plus_q
    real(8)              :: yk,yk_plus_q
    real(8)              :: SimEk,SimEk_plus_q
    real(8)              :: Dk,Dk_plus_q
    real(8)              :: num(2,2),den
    real(8)              :: num_t(2,2),den_t, Mk_t, Mkq_t, Ek2, Ekq2    
    real(8)              :: tail_num(2,2), tail_den, tail, nu       
    integer              :: ik,n
    real(8)              :: kx,ky,qx,qy,vkq
    !
    ReSigma   = p_work(1:L)
    ImSigma   = p_work(L+1:2*L)
    Tz        = p_work(2*L+1)
    dPlus     = p_work(2*L+2)
    dMinus    = p_work(2*L+3)
    !
    Chi       = 0d0
    qx        = qvec_work(1)
    qy        = qvec_work(2)
    nu        = 2*m*pi/beta
    do ik=1,Nktot
       kvec      = Kgrid(ik,:)
       kx        = kvec(1)
       ky        = kvec(2)
       !
       ek        = -1d0*(cos(kx)+cos(ky))
       ek_plus_q = -1d0*(cos(kx+qx)+cos(ky+qy))
       xk        = lambda*sin(kx)
       yk        = lambda*sin(ky)
       xk_plus_q = lambda*sin(kx+qx)
       yk_plus_q = lambda*sin(ky+qy)
       vkq       = xk*xk_plus_q + yk*yk_plus_q
       !
       Mk_t      = Mh - Tz*gt + ek
       Mkq_t     = Mh - Tz*gt + ek_plus_q          
       !
       Ek2       = Mk_t**2 + xk**2 + yk**2
       Ekq2      = Mkq_t**2 + xk_plus_q**2 + yk_plus_q**2          
       !
       do n=1,Lf
          wn           = pi/beta*(2*n-1)
          wnm          = pi/beta*(2*(n+m)-1)
          !
          zn           = wn-ImSigma(n) 			
          zn_plus_m    = wnm-ImSigma(n+m)
          !
          Mk           = Mh - Tz*gt + ReSigma(n) + ek
          Mk_plus_q    = Mh - Tz*gt + ReSigma(n+m) + ek_plus_q
          !
          simEk        = Mk**2 + xk**2 + yk**2
          simEk_plus_q = Mk_plus_q**2 + xk_plus_q**2 + yk_plus_q**2
          !
          Dk           = zn**2d0 + simEk
          Dk_plus_q    = zn_plus_m**2 + simEk_plus_q
          !
          den          = Dk*Dk_plus_q
          num(1,1)     = -zn*zn_plus_m + Mk*Mk_plus_q
          num(2,2)     = vkq
          num(1,2)     = zn*Mk_plus_q  + zn_plus_m*Mk
          num(2,1)     = num(1,2)
          !
          den_t        = (wnm**2d0 + Ekq2)*(wn**2d0 + Ek2)
          num_t(1,1)   = -wn*wnm + Mk_t*Mkq_t
          num_t(2,2)   =  vkq
          num_t(1,2)   = zero
          num_t(2,1)   = zero
          !
          Chi(1,1)     = Chi(1,1) - 4d0*(num(1,1)/den-num_t(1,1)/den_t)/beta/Nktot
          Chi(2,2)     = Chi(2,2) - 4d0*(num(2,2)/den-num_t(2,2)/den_t)/beta/Nktot
          Chi(1,2)     = Chi(1,2) - 8d0*num(1,2)/den/beta/Nktot
          Chi(2,1)     = Chi(1,2)
       enddo
       tail_den      = sqrt(Ek2)*((Ekq2-Ek2)**2d0 + 2d0*(Ekq2+Ek2)*nu**2d0 + nu**4d0)
       tail_num(1,1) = (Ek2-Ekq2)*(Ek2 + Mk_t*Mkq_t) + (Ek2 - Mk_t*Mkq_t)*nu**2d0
       tail_num(2,2) = vkq*(Ek2 - Ekq2 - nu**2d0)
       tail_num(1,2) = zero
       tail_num(2,1) = zero
       !
       !Chi = Chi + Tail
       Chi      = Chi + 2d0*tanh(0.5d0*beta*sqrt(Ek2))*tail_num/tail_den/Nktot
    enddo
  end subroutine Intk_SumMats_Chi_qv







  subroutine build_self_energy(p,Sigma)
    real(8),dimension(:)                          :: p
    complex(8),dimension(Nspin,Nspin,Norb,Norb,L) :: Sigma
    complex(8),dimension(Nspin*Norb,Nspin*Norb,L) :: Self
    do i=1,L
       Self(:,:,i)      = -gt*p(2*L+1)*Gamma5 + p(i)*Gamma5 + xi*p(L+i)*GammaN
       Sigma(:,:,:,:,i) = j2so(Self(:,:,i))
    enddo
  end subroutine build_self_energy







  !#################################
  !#################################
  !#################################
  !#################################



  function hk_bhz(kvec,N) result(hk)
    integer                   :: N
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: ek,kx,ky
    integer                   :: ii
    if(N/=Nso)stop "hk_bhz error: N != Nspin*Norb == 4"
    kx=kvec(1)
    ky=kvec(2)
    ek = -1d0*(cos(kx)+cos(ky))
    Hk = (Mh+ek)*Gamma5 + lambda*sin(kx)*Gamma1 + lambda*sin(ky)*Gamma2
  end function hk_bhz



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



  function get_simplified_gf(kvec,i,N,Sigma) result(gk)
    real(8),dimension(:)      :: kvec
    complex(8)                :: z
    integer                   :: i,N
    complex(8),dimension(N,N) :: gk,sigma
    real(8)                   :: kx,ky
    real(8)                   :: w_,M_,x_,y,Ek

    kx = kvec(1)
    ky = kvec(2)
    !
    w_ = pi/beta*(2*i-1) - dimag(sigma(1,1))
    M_ = Mh-1d0*(cos(kx)+cos(ky)) + dreal(sigma(1,1))
    x_ = lambda*sin(kx) 
    y  = lambda*sin(ky)
    !
    Ek = M_**2 + x_**2 + y**2 
    Gk = xi*w_*eye(Nso) + M_*Gamma5 + x_*Gamma1  + y*Gamma2 
    Gk = Gk/(-w_**2-Ek)
  end function get_simplified_gf





  function trapz2d_system(m,func,xrange,yrange,Nx,Ny) result(int)
    integer :: m
    interface
       function func(x,m)
         real(8),dimension(:) :: x
         integer              :: m
         real(8)              :: func(m)
       end function func
    end interface
    real(8),dimension(2) :: xrange,yrange
    integer              :: Nx,Ny,i,j
    real(8)              :: xx(Nx),yy(Ny)
    real(8)              :: hx,hy
    real(8)              :: int(m)
    hx=xrange(2)-xrange(1)
    hx=hx/Nx
    hy=yrange(2)-yrange(1)
    hy=hy/Ny
    int=&
         func([xrange(1),yrange(1)],m)+&
         func([xrange(1),yrange(2)],m)+&
         func([xrange(2),yrange(1)],m)+&
         func([xrange(2),yrange(2)],m)
    xx=linspace(xrange(1),xrange(2),Nx,iend=.false.)
    yy=linspace(yrange(1),yrange(2),Ny,iend=.false.)
    do i=2,Nx
       do j=2,Ny
          int=int+4d0*func([xx(i),yy(j)],m)
       enddo
    enddo
    do j=2,Ny
       int=int+2d0*( func([xrange(1),yy(j)],m) + func([xrange(2),yy(j)],m) )
    enddo
    do i=2,Nx
       int=int+2d0*( func([xx(i),yrange(1)],m) + func([xx(i),yrange(2)],m) )
    enddo
    int=int*hx*hy/4d0
  end function trapz2d_system








  function trapz2d_func(func,xrange,yrange,Nx,Ny) result(int)
    interface
       function func(x)
         real(8),dimension(:) :: x
         real(8)              :: func
       end function func
    end interface
    real(8),dimension(2) :: xrange,yrange
    integer              :: Nx,Ny,i,j
    real(8)              :: xx(Nx),yy(Ny)
    real(8)              :: hx,hy
    real(8)              :: int
    hx=xrange(2)-xrange(1)
    hx=hx/Nx
    hy=yrange(2)-yrange(1)
    hy=hy/Ny
    int=&
         func([xrange(1),yrange(1)])+&
         func([xrange(1),yrange(2)])+&
         func([xrange(2),yrange(1)])+&
         func([xrange(2),yrange(2)])
    xx=linspace(xrange(1),xrange(2),Nx,iend=.false.)
    yy=linspace(yrange(1),yrange(2),Ny,iend=.false.)
    do i=2,Nx
       do j=2,Ny
          int=int+4d0*func([xx(i),yy(j)])
       enddo
    enddo
    do j=2,Ny
       int=int+2d0*( func([xrange(1),yy(j)]) + func([xrange(2),yy(j)]) )
    enddo
    do i=2,Nx
       int=int+2d0*( func([xx(i),yrange(1)]) + func([xx(i),yrange(2)]) )
    enddo
    int=int*hx*hy/4d0
  end function trapz2d_func


end program bhz_2d








! subroutine solve_eqs(p) 
!   real(8),dimension(:),intent(inout) :: p ![2+2L]
!   real(8)                            :: Tz,dTz,TzTmp,dTzTmp
!   real(8),dimension(L)               :: ReSigma,ImSigma
!   real(8)                            :: kvec(2),qvec(2)
!   real(8)                            :: wn
!   real(8)                            :: ReS(L),ImS(L)
!   real(8)                            :: Meff,ek,xk,yk,SimEk
!   real(8)                            :: Den,ChiTmp
!   complex(8)                         :: N_el 
!   complex(8)                         :: Smats(Nspin,Nspin,Norb,Norb,L)
!   complex(8)                         :: Gkmats(Nso,Nso,L)
!   real(8)                            :: n_k(Nktot,Nso)
!   integer                            :: ik,n,m


!   !split params as required:
!   ReSigma = p(1:L)
!   ImSigma = p(L+1:2*L)
!   Tz      = p(2*L+1)
!   dTz     = p(2*L+2)
!   !
!   ReS     = 0d0
!   ImS     = 0d0
!   TzTmp   = 0d0
!   dTzTmp  = 0d0
!   !
!   do ik=1,Nktot
!      kvec = Kgrid(ik,:)
!      qvec = Kgrid(ik,:)
!      ek   = -1d0*(cos(kvec(1))+cos(kvec(2)))
!      xk   = lambda*sin(kvec(1))
!      yk   = lambda*sin(kvec(2))
!      !
!      do n=1,L
!         wn     = pi/beta*(2*n-1)-ImSigma(n) 		
!         Meff   = Mh - Tz*gt + ReSigma(n)
!         simEk  = (Meff + ek)**2 + xk**2 + yk**2
!         Den    = wn**2d0 + simEk
!         ReS(n) = ReS(n) + (Meff + ek)/Den
!         ImS(n) = ImS(n) + wn/Den
!         TzTmp  = TzTmp  + (Meff + ek)/Den  	
!      enddo
!      !
!      do m=1,Lb
!         ChiTmp = Chi_qv(p,qvec,m)
!         dTzTmp = dTzTmp + ChiTmp/(1d0-gt*ChiTmp)                   
!      enddo
!   enddo
!   !
!   ReSigma = -ReS*dTz/Nktot
!   ImSigma = -ImS*dTz/Nktot
!   Tz      = -4d0*TzTmp/beta/Nktot 
!   dTz     =  2d0*dTzTmp*gt**2d0/beta/Nktot
!   !
!   !Update params:
!   p(1:L)     = ReSigma
!   p(L+1:2*L) = ImSigma
!   p(2*L+1)   = Tz
!   p(2*L+2)   = dTz


!   call build_self_energy(p,Smats)
!   do ik=1,Nktot
!      kvec = Kgrid(ik,:)
!      do i=1,L              
!         Gkmats(:,:,i) = get_simplified_gf(kvec,i,Nso,so2j(Smats(:,:,:,:,i)))          
!      end do
!      do i=1,Nso
!         n_k(ik,i) = fft_get_density(Gkmats(i,i,:),beta)
!      end do
!   end do

!   N_el = sum(n_k)/Nktot

!   write(*,*)iter,Tz,dTz,ReSigma(1),ImSigma(1), N_el
!   call splot("Sigma_iw_iter"//str(iter,3)//".dat",pi/beta*(2*arange(1,L)-1),dcmplx(ReSigma(:),ImSigma(:)))
!   return
! end subroutine solve_eqs



! function Chi_qv(p,qvec,m) result(chi)
!   real(8),dimension(:),intent(in) :: p
!   real(8),dimension(:),intent(in) :: qvec
!   integer,intent(in)              :: m
!   real(8)		   	    :: chi
!   real(8)                         :: Tz,dTz
!   real(8),dimension(L)            :: ReSigma,ImSigma
!   real(8)                         :: kvec(2)
!   real(8)                         :: wn,wn_plus_m
!   real(8)                         :: Mk,Mk_plus_q
!   real(8)                         :: ek,ek_plus_q
!   real(8)                         :: xk,xk_plus_q
!   real(8)                         :: yk,yk_plus_q
!   real(8)                         :: SimEk,SimEk_plus_q
!   real(8)                         :: Dk,Dk_plus_q
!   real(8)                         :: num,den
!   real(8)                         :: num_t,den_t, Mk_t, Mkq_t, Ek2, Ekq2    
!   real(8)                         :: tail_num, tail_den, tail, nu       
!   integer                         :: ik,n
!   real(8)                         :: kx,ky,qx,qy,vkq

!   ReSigma = p(1:L)
!   ImSigma = p(L+1:2*L)
!   Tz      = p(2*L+1)
!   dTz     = p(2*L+2)
!   Chi     = 0d0
!   tail    = 0d0

!   qx = qvec(1)
!   qy = qvec(2)
!   nu = 2*m*pi/beta

!   do ik=1,Nktot
!      kx        = Kgrid(ik,1)
!      ky        = Kgrid(ik,2)     
!      ek        = -1d0*(cos(kx)+cos(ky))
!      ek_plus_q = -1d0*(cos(kx+qx)+cos(ky+qy))
!      xk        = lambda*sin(kx)
!      yk        = lambda*sin(ky)
!      xk_plus_q = lambda*sin(kx+qx)
!      yk_plus_q = lambda*sin(ky+qy)
!      vkq       = xk*xk_plus_q + yk*yk_plus_q
!      !
!      Mk_t      = Mh - Tz*gt + ek
!      Mkq_t     = Mh - Tz*gt + ek_plus_q          
!      !
!      Ek2       = Mk_t**2 + xk**2 + yk**2
!      Ekq2      = Mkq_t**2 + xk_plus_q**2 + yk_plus_q**2          
!      !       
!      do n=1,Lf
!         wn            = pi/beta*(2*n-1)-ImSigma(n) 			
!         wn_plus_m     = pi/beta*(2*(n+m)-1)-ImSigma(n+m) 	 
!         !
!         Mk            = Mh - Tz*gt + ReSigma(n) + ek
!         Mk_plus_q     = Mh - Tz*gt + ReSigma(n+m) + ek_plus_q
!         !
!         simEk         = Mk**2 + xk**2 + yk**2
!         simEk_plus_q  = Mk_plus_q**2 + xk_plus_q**2 + yk_plus_q**2
!         !
!         Dk            = wn**2d0 + simEk
!         Dk_plus_q     = wn_plus_m**2 + simEk_plus_q
!         !
!         num = wn*wn_plus_m - Mk*Mk_plus_q + vkq
!         den = Dk*Dk_plus_q
!         !
!         num_t  = (pi/beta*(2*n-1))*(pi/beta*(2*(n+m)-1)) - Mk_t*Mkq_t + vkq
!         den_t  = ((pi/beta*(2*(n+m)-1))**2d0 + Ekq2)*((pi/beta*(2*n-1))**2d0 + Ek2)
!         !
!         Chi = Chi + 2d0*(num/den-num_t/den_t)/beta/Nktot
!      enddo

!      tail_num = (Ekq2-Ek2)*(vkq-Mk_t*Mkq_t-Ek2) + (vkq-Mk_t*Mkq_t+Ek2)*nu**2d0
!      tail_den = sqrt(Ek2)*((Ekq2-Ek2)**2d0 + 2d0*(Ekq2+Ek2)*nu**2d0 + nu**4d0)

!      tail = tail + tanh(0.5d0*beta*sqrt(Ek2))*tail_num/tail_den/Nktot                  
!   enddo

!   Chi = Chi + tail

! end function Chi_qv









! function SumMats_Chi_kv_q(kvec) result(chi)
!   real(8),dimension(:) :: kvec
!   real(8)              :: chi
!   real(8)              :: Tz,dTz
!   real(8),dimension(L) :: ReSigma,ImSigma
!   real(8)              :: iwn,iwnm
!   real(8)              :: wn,wn_plus_m
!   real(8)              :: Mk,Mk_plus_q
!   real(8)              :: ek,ek_plus_q
!   real(8)              :: xk,xk_plus_q
!   real(8)              :: yk,yk_plus_q
!   real(8)              :: SimEk,SimEk_plus_q
!   real(8)              :: Dk,Dk_plus_q
!   real(8)              :: num,den
!   real(8)              :: num_t,den_t, Mk_t, Mkq_t, Ek2, Ekq2    
!   real(8)              :: tail_num, tail_den, tail, nu       
!   integer              :: ik,n
!   real(8)              :: kx,ky,qx,qy,vkq
!   !
!   ReSigma   = p_work(1:L)
!   ImSigma   = p_work(L+1:2*L)
!   Tz        = p_work(2*L+1)
!   dTz       = p_work(2*L+2)
!   m         = m_work
!   !
!   kx        = kvec(1)
!   ky        = kvec(2)
!   qx        = qvec_work(1)
!   qy        = qvec_work(2)
!   nu        = 2*m*pi/beta
!   !
!   ek        = -1d0*(cos(kx)+cos(ky))
!   ek_plus_q = -1d0*(cos(kx+qx)+cos(ky+qy))
!   xk        = lambda*sin(kx)
!   yk        = lambda*sin(ky)
!   xk_plus_q = lambda*sin(kx+qx)
!   yk_plus_q = lambda*sin(ky+qy)
!   vkq       = xk*xk_plus_q + yk*yk_plus_q
!   !
!   Mk_t      = Mh - Tz*gt + ek
!   Mkq_t     = Mh - Tz*gt + ek_plus_q          
!   !
!   Ek2       = Mk_t**2 + xk**2 + yk**2
!   Ekq2      = Mkq_t**2 + xk_plus_q**2 + yk_plus_q**2          
!   !
!   Chi = 0d0
!   do n=1,Lf
!      iwn          = pi/beta*(2*n-1)
!      iwnm         = pi/beta*(2*(n+m)-1)
!      !
!      wn           = iwn-ImSigma(n) 			
!      wn_plus_m    = iwnm-ImSigma(n+m)
!      !
!      Mk           = Mh - Tz*gt + ReSigma(n) + ek
!      Mk_plus_q    = Mh - Tz*gt + ReSigma(n+m) + ek_plus_q
!      !
!      simEk        = Mk**2 + xk**2 + yk**2
!      simEk_plus_q = Mk_plus_q**2 + xk_plus_q**2 + yk_plus_q**2
!      !
!      Dk           = wn**2d0 + simEk
!      Dk_plus_q    = wn_plus_m**2 + simEk_plus_q
!      !
!      num          = wn*wn_plus_m - Mk*Mk_plus_q + vkq
!      den          = Dk*Dk_plus_q
!      !
!      num_t        = iwn*iwnm - Mk_t*Mkq_t + vkq
!      den_t        = (iwnm**2d0 + Ekq2)*(iwn**2d0 + Ek2)
!      !
!      Chi          = Chi + 2d0*(num/den-num_t/den_t)/beta
!   enddo
!   tail_num = (Ekq2-Ek2)*(vkq-Mk_t*Mkq_t-Ek2) + (vkq-Mk_t*Mkq_t+Ek2)*nu**2d0
!   tail_den = sqrt(Ek2)*((Ekq2-Ek2)**2d0 + 2d0*(Ekq2+Ek2)*nu**2d0 + nu**4d0)
!   !Chi = Chi + Tail
!   Chi      = Chi + tanh(0.5d0*beta*sqrt(Ek2))*tail_num/tail_den
! end function SumMats_Chi_kv_q





