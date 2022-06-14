MODULE GLOBAL_VARIABLES
  implicit none
  !SET THE DIMENSION OF THE PROBLEM
  integer,parameter                  :: Norb=2
  integer,parameter                  :: Nspin=2
  integer,parameter                  :: Nso=Nspin*Norb
  !
  !INPUT VARIABLES
  integer                            :: Nkx,Nky,Nktot
  integer                            :: L,Lf,Lb
  integer 			     :: Iter
  integer 			     :: MaxIter
  integer 			     :: Nsuccess=2  
  real(8)                            :: gt
  real(8)                            :: gn
  real(8)                            :: mh
  real(8)                            :: lambda
  real(8)                            :: xmu
  real(8)                            :: beta
  real(8)                            :: eps
  real(8)                            :: wmix
  real(8)                            :: it_error
  real(8)                            :: tz0,dtz0
  logical                            :: withgf
  !
  !Gamma matrices:
  complex(8),dimension(Nso,Nso)      :: Gamma1,Gamma2,Gamma5,GammaN
  !
  !GLOBALLY SHARED VARIABLES:
  real(8)                            :: Tz !the orbital polarization
  real(8)                            :: Ne !the total density
  real(8),dimension(:,:),allocatable :: kgrid !K vector grid
  !
  !WORK ARRAYS:
  real(8),dimension(:),allocatable   :: pwork
  real(8)                            :: qvec_work(2)
  integer                            :: m_work
END MODULE GLOBAL_VARIABLES





program bhz_2d
  !LOCAL:
  USE GLOBAL_VARIABLES
  !LIBRARIES:
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
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
  real(8),dimension(:),allocatable            :: wmats


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
  call parse_input_variable(it_error,"IT_ERROR",Finput,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",Finput,default=100)    
  call parse_input_variable(eps,"EPS",Finput,default=4.d-2)
  call parse_input_variable(wmix,"WMIX",Finput,default=1d0)
  call parse_input_variable(withgf,"WITHGF",Finput,default=.false.)
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


  gamma1 = kron_pauli( pauli_sigma_z, pauli_tau_x)
  gamma2 = kron_pauli( pauli_sigma_0,-pauli_tau_y)
  gamma5 = kron_pauli( pauli_sigma_0, pauli_tau_z)
  gammaN = kron_pauli( pauli_sigma_0, pauli_tau_0)
  Nky    = Nkx
  Nktot  = Nkx*Nky
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




  !> SOLVE FLUCTUATIONS:
  !Setup the k-space lattice basis:
  call TB_set_bk([pi2,0d0],[0d0,pi2])

  allocate(kgrid(Nktot,2))      !Nktot=# tot kpoints, 2= 2D
  call TB_build_kgrid([Nkx,Nky],kgrid)
  allocate(wmats(L))
  wmats = pi/beta*(2*arange(1,L)-1)



  !+ ReSigma(iw_n)[L] + ImSIgma(iw_n)[L] + Tz[1] + <|dTz|**2>[1] 
  Nparams = 2 + 2*L
  allocate( params(Nparams), params_prev(Nparams), pwork(Nparams))

  !Start from MF solution
  params = [dble(zeros(L)),dble(zeros(L)),Tz,abs(dTz0)]
  inquire(file="params.restart",exist=iexist)
  if(iexist)call read_array("params.restart",params)
  call save_array("params.init",params)

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
  write(unit,*)gt,gn,params(2*L+1),params(2*L+2)
  close(unit)
  write(*,*) "Tz,dTz=",params(2*L+1),params(2*L+2)


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
     write(unit,*) 2d0*pi*m/beta, Chi_qv(params,[0d0,0d0],m)
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
    f = gt/2d0*(Tz**2) - 2d0*integral/pi2/pi2
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
    Fk  = sqrt( (mh - gt*Tz/2d0 + ek)**2 + (x2+y2) )
  end function fmf_bhz
  !#################################








  subroutine solve_eqs(p) 
    real(8),dimension(:),intent(inout) :: p ![2+2L]
    real(8)                            :: Tz,dTz,TzTmp,dTzTmp
    real(8),dimension(L)               :: ReSigma,ImSigma
    real(8)                            :: kvec(2),qvec(2)
    real(8)                            :: wn
    real(8)                            :: ReS(L),ImS(L)
    real(8)                            :: ReG(L),ImG(L),denSigma(L)
    real(8)                            :: Meff,ek,xk,yk,SimEk
    real(8)                            :: Den,ChiTmp
    complex(8)                         :: N_el 
    complex(8)                         :: Smats(Nspin,Nspin,Norb,Norb,L)
    complex(8)                         :: Gkmats(Nso,Nso,L)
    real(8)                            :: n_k(Nktot,Nso)
    integer                            :: ik,n,m


    !split params as required:
    ReSigma = p(1:L)
    ImSigma = p(L+1:2*L)
    Tz      = p(2*L+1)
    dTz     = p(2*L+2)
    !
    ReG     = 0d0
    ImG     = 0d0
    ReS     = 0d0
    ImS     = 0d0
    TzTmp   = 0d0
    dTzTmp  = 0d0
    !
    do ik=1,Nktot
       kvec = Kgrid(ik,:)
       qvec = Kgrid(ik,:)
       ek   = -1d0*(cos(kvec(1))+cos(kvec(2)))
       xk   = lambda*sin(kvec(1))
       yk   = lambda*sin(kvec(2))
       !
       do n=1,L
          wn     = pi/beta*(2*n-1)-ImSigma(n) 		
          Meff   = Mh - Tz*gt/2d0 + ReSigma(n)
          simEk  = (Meff + ek)**2 + xk**2 + yk**2
          Den    = wn**2d0 + simEk
          ! ReG(n) = ReG(n) + (Meff + ek)/Den  
          ! ImG(n) = ImG(n) + wn/Den
          ReS(n) = ReS(n) + (Meff + ek)/Den
          ImS(n) = ImS(n) + wn/Den
          TzTmp  = TzTmp  + (Meff + ek)/Den  	
       enddo
       !
       do m=1,Lb
          ChiTmp = Chi_qv(p,qvec,m)
          dTzTmp = dTzTmp + ChiTmp/(1d0-gt*ChiTmp)                   
       enddo
    enddo
    !
    ReG      = -ReG/Nktot
    ImG      = -ImG/Nktot
    ! denSigma = (1d0-gt*Tz*ReG)**2d0 + (gt*Tz*ImG)**2d0
    ! ReSigma  = dTz*(ReG - gt*Tz*(ReG**2d0+ImG**2d0))/denSigma
    ! ImSigma  = dTz*ImG/denSigma    
    ReSigma = -ReS*dTz/Nktot
    ImSigma = -ImS*dTz/Nktot
    Tz      = -4d0*TzTmp/beta/Nktot 
    dTz     =  2d0*dTzTmp*gt**2d0/beta/Nktot
    !
    !Update params:
    p(1:L)     = ReSigma
    p(L+1:2*L) = ImSigma
    p(2*L+1)   = Tz
    p(2*L+2)   = dTz


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

    write(*,*)iter,Tz,dTz,ReSigma(1),ImSigma(1), N_el
    call splot("Sigma_iw_iter"//str(iter,3)//".dat",wmats,dcmplx(ReSigma(:),ImSigma(:)))
    return
  end subroutine solve_eqs






  function Chi_qv(p,qvec,m) result(chi)
    real(8),dimension(:),intent(in) :: p
    real(8),dimension(:),intent(in) :: qvec
    integer,intent(in)              :: m
    real(8)		   	    :: chi,chi_tmp
    !Push required parameters to integrand
    pwork     = p
    qvec_work = qvec
    m_work    = m
    call gauss_quad(chi_qv_k, [0d0,0d0],[pi2,pi2],ans=chi_tmp)
    chi = chi_tmp/pi2/pi2
  end function Chi_qv



  function chi_qv_k(kvec) result(chi)
    real(8),dimension(:) :: kvec
    real(8)              :: chi
    real(8)              :: Tz,dTz
    real(8),dimension(L) :: ReSigma,ImSigma
    real(8)              :: iwn,iwnm
    real(8)              :: wn,wn_plus_m
    real(8)              :: Mk,Mk_plus_q
    real(8)              :: ek,ek_plus_q
    real(8)              :: xk,xk_plus_q
    real(8)              :: yk,yk_plus_q
    real(8)              :: SimEk,SimEk_plus_q
    real(8)              :: Dk,Dk_plus_q
    real(8)              :: num,den
    real(8)              :: num_t,den_t, Mk_t, Mkq_t, Ek2, Ekq2    
    real(8)              :: tail_num, tail_den, tail, nu       
    integer              :: ik,n
    real(8)              :: kx,ky,qx,qy,vkq
    !
    ReSigma   = pwork(1:L)
    ImSigma   = pwork(L+1:2*L)
    Tz        = pwork(2*L+1)
    dTz       = pwork(2*L+2)
    !
    kx        = kvec(1)
    ky        = kvec(2)
    qx        = qvec_work(1)
    qy        = qvec_work(2)
    nu        = 2*m_work*pi/beta
    !
    ek        = -1d0*(cos(kx)+cos(ky))
    ek_plus_q = -1d0*(cos(kx+qx)+cos(ky+qy))
    xk        = lambda*sin(kx)
    yk        = lambda*sin(ky)
    xk_plus_q = lambda*sin(kx+qx)
    yk_plus_q = lambda*sin(ky+qy)
    vkq       = xk*xk_plus_q + yk*yk_plus_q
    !
    Mk_t      = Mh - Tz*gt/2d0 + ek
    Mkq_t     = Mh - Tz*gt/2d0 + ek_plus_q          
    !
    Ek2       = Mk_t**2 + xk**2 + yk**2
    Ekq2      = Mkq_t**2 + xk_plus_q**2 + yk_plus_q**2          
    !       
    do n=1,Lf
       iwn          = pi/beta*(2*n-1)
       iwnm         = pi/beta*(2*(n+m)-1)
       !
       wn           = iwn-ImSigma(n) 			
       wn_plus_m    = iwnm-ImSigma(n+m) 	 
       !
       Mk           = Mh - Tz*gt/2d0 + ReSigma(n) + ek
       Mk_plus_q    = Mh - Tz*gt/2d0 + ReSigma(n+m) + ek_plus_q
       !
       simEk        = Mk**2 + xk**2 + yk**2
       simEk_plus_q = Mk_plus_q**2 + xk_plus_q**2 + yk_plus_q**2
       !
       Dk           = wn**2d0 + simEk
       Dk_plus_q    = wn_plus_m**2 + simEk_plus_q
       !
       num          = wn*wn_plus_m - Mk*Mk_plus_q + vkq
       den          = Dk*Dk_plus_q
       !
       num_t        = iwn*iwnm - Mk_t*Mkq_t + vkq
       den_t        = (iwnm**2d0 + Ekq2)*(iwn**2d0 + Ek2)
       !
       Chi          = Chi + 2d0*(num/den-num_t/den_t)/beta/Nktot
    enddo
    tail_num = (Ekq2-Ek2)*(vkq-Mk_t*Mkq_t-Ek2) + (vkq-Mk_t*Mkq_t+Ek2)*nu**2d0
    tail_den = sqrt(Ek2)*((Ekq2-Ek2)**2d0 + 2d0*(Ekq2+Ek2)*nu**2d0 + nu**4d0)
    Chi      = Chi + tanh(0.5d0*beta*sqrt(Ek2))*tail_num/tail_den/Nktot                  
  end function chi_qv_k






  subroutine build_self_energy(p,Sigma)
    real(8),dimension(:)                          :: p
    complex(8),dimension(Nspin,Nspin,Norb,Norb,L) :: Sigma
    complex(8),dimension(Nspin*Norb,Nspin*Norb,L) :: Self
    do i=1,L
       Self(:,:,i)      = -gt*p(2*L+1)/2d0*Gamma5 + p(i)*Gamma5 + xi*p(L+i)*GammaN
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


end program bhz_2d










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
!      Mk_t      = Mh - Tz*gt/2d0 + ek
!      Mkq_t     = Mh - Tz*gt/2d0 + ek_plus_q          
!      !
!      Ek2       = Mk_t**2 + xk**2 + yk**2
!      Ekq2      = Mkq_t**2 + xk_plus_q**2 + yk_plus_q**2          
!      !       
!      do n=1,Lf
!         wn            = pi/beta*(2*n-1)-ImSigma(n) 			
!         wn_plus_m     = pi/beta*(2*(n+m)-1)-ImSigma(n+m) 	 
!         !
!         Mk            = Mh - Tz*gt/2d0 + ReSigma(n) + ek
!         Mk_plus_q     = Mh - Tz*gt/2d0 + ReSigma(n+m) + ek_plus_q
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
