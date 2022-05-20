program bhz_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none


  integer,parameter                       :: Norb=2,Nspin=2,Nso=Nspin*Norb
  integer				  :: Nparams
  integer                                 :: Nk,Nktot,Nkx,Npts,Lf,Lb,L
  integer                                 :: Nky
  integer                                 :: i,j,k,ik
  integer                                 :: info,unit
  integer 				  :: Iter,MaxIter,Nsuccess=2  
  real(8)                                 :: Uloc,Jh,Jhratio,Sz,Tz
  real(8)                                 :: mh,lambda
  real(8)                                 :: xmu,beta,eps
  real(8)                                 :: wmix,it_error,tz0,dtz0
  real(8)                                 :: tol,sigma
  real(8)                                 :: gint,x(1),dx(1)
  real(8)                                 :: integral
  real(8)                                 :: ky_g
  logical                                 :: withgf
  logical                                 :: iexist,converged
  character(len=20)                       :: Finput
  complex(8)                              :: Hloc(Nso,Nso),Hmf_glob(Nso,Nso)
  !
  real(8),dimension(:,:),allocatable      :: kgrid,kpath
  complex(8),dimension(:,:,:),allocatable :: Hk
  complex(8),dimension(Nso,Nso)           :: Gamma1,Gamma2,Gamma5,GammaN
  real(8),dimension(:),allocatable        :: params,params_prev,wmats





  call parse_cmd_variable(Finput,"FINPUT",default="input.conf")
  call parse_input_variable(Lf,"LF",Finput,default=256,comment="# of fermionic Mats frequencies, L=Lf+Lb")
  call parse_input_variable(Lb,"LB",Finput,default=64,comment="# of bosonix Mats frequencies, L=Lf+Lb")
  call parse_input_variable(tz0,"tz0",Finput,default=-0.1d0,comment="Guess for MF search of Tz (tz0<0)")
  call parse_input_variable(dtz0,"dtz0",Finput,default=0.1d0,comment="Guess for dTz fluctuations (dtz0>0)")
  call parse_input_variable(nkx,"NKX",Finput,default=10)
  call parse_input_variable(Uloc,"ULOC",Finput,default=0d0)
  call parse_input_variable(Jhratio,"JHRATIO",Finput,default=0d0)
  call parse_input_variable(Mh,"MH",Finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",Finput,default=0.3d0)
  call parse_input_variable(xmu,"XMU",Finput,default=0.d0)  
  call parse_input_variable(beta,"BETA",Finput,default=1000.d0)
  call parse_input_variable(it_error,"IT_ERROR",Finput,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",Finput,default=100)    
  call parse_input_variable(eps,"EPS",Finput,default=4.d-2)
  call parse_input_variable(wmix,"WMIX",Finput,default=1d0)

  call parse_input_variable(withgf,"WITHGF",Finput,default=.false.)
  !
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
  Jh     = Jhratio*Uloc
  gint   = Uloc*(1d0-5d0*Jhratio)
  L      = Lf+Lb
  write(*,*)"Using L freq.=",L


  !>SOLVE MF PROBLEM 1st: >>ACTHUNG<< This solution does not use BZ basis defined later!!
  call start_timer()
  x(1)=-abs(tz0)
  dx(1)=0.1d0
  call fmin(bhz_f,x,lambda=dx)
  tz=x(1)
  open(free_unit(unit),file="mf_tzVSuloc.dat")
  write(unit,*)uloc,tz
  close(unit)
  write(*,*) "Tz=",tz
  call stop_timer(" Mean-Field")



  !> SOLVE FLUCTUATIONS:
  !Setup the k-space lattice basis:
  call TB_set_bk([pi2,0d0],[0d0,pi2])


  allocate(kgrid(Nktot,2))      !Nktot=# tot kpoints, 2= 2D
  call TB_build_kgrid([Nkx,Nky],kgrid)
  allocate(wmats(L))
  wmats = pi/beta*(2*arange(1,L)-1)



  !Tz[1] + <|dTz|**2>[1] + ReSigma(iw_n)[L] + ImSIgma(iw_n)[L]
  Nparams = 2 + 2*L
  allocate( params(Nparams), params_prev(Nparams))

  !Start from MF solution
  params = [dble(zeros(L)),dble(zeros(L)),Tz,abs(dTz0)]

  inquire(file="params.restart",exist=iexist)
  if(iexist)call read_array("params.restart",params)
  call save_array("params.init",params) 	!FP: verifying the initial params

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
  open(free_unit(unit),file="tz_dtzVSuloc.dat")
  write(unit,*)uloc,params(2*L+1),params(2*L+2)
  close(unit)
  write(*,*) "Tz,dTz=",params(2*L+1),params(2*L+2)



contains





  !--------------------------------------------------------------------!
  !BHZ HAMILTONIAN:
  !--------------------------------------------------------------------!
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



  subroutine solve_eqs(p) 
    real(8),dimension(:),intent(inout) :: p ![2+2L]
    real(8)                         :: Tz,dTz,TzTmp,dTzTmp
    real(8),dimension(L)            :: ReSigma,ImSigma
    real(8)                         :: kvec(2),qvec(2)
    real(8)                         :: wn
    real(8)                         :: ReS(L),ImS(L)
    real(8)                         :: Meff,ek,xk,yk,SimEk
    real(8)                         :: Den,ChiTmp
    integer                         :: ik,n,m

    !split params as required:
    ReSigma = p(1:L)
    ImSigma = p(L+1:2*L)
    Tz      = p(2*L+1)
    dTz     = p(2*L+2)
    ReS     = 0d0
    ImS     = 0d0
    TzTmp   = 0d0
    dTzTmp  = 0d0
    do ik=1,Nktot
       kvec = Kgrid(ik,:)
       qvec = Kgrid(ik,:)
       ek   = -1d0*(cos(kvec(1))+cos(kvec(2)))
       xk   = lambda*sin(kvec(1))
       yk   = lambda*sin(kvec(2))
       !
       do n=1,L
          wn     = pi/beta*(2*n-1)-ImSigma(n)
          Meff   = Mh - Tz*gint/2d0 + ReSigma(n)
          simEk  = (Meff + ek)**2 + xk**2 + yk**2
          Den    = wn**2d0 + simEk
          ReS(n) = ReS(n)  + (Meff + ek)/Den  
          ImS(n) = ImS(n)  + wn/Den
          TzTmp  = TzTmp   + (Meff + ek)/Den  
       enddo
       !
       do m=1,Lb
          ChiTmp = Chi_qv(p,qvec,m)
          dTzTmp = dTzTmp + ChiTmp/(1d0-gint*ChiTmp)
       enddo
    enddo
    ReSigma = -ReS*dTz/Nktot
    ImSigma = -ImS*dTz/Nktot
    Tz      = -TzTmp*4d0*gint/beta/Nktot
    dTz     =  2d0*dTzTmp*gint**2d0/beta/Nktot
    !
    !Update params:
    p(1:L)     = ReSigma
    p(L+1:2*L) = ImSigma
    p(2*L+1)   = Tz
    p(2*L+2)   = dTz

    write(*,*)iter,Tz,dTz,ReSigma(1),ImSigma(1)
    call splot("Sigma_iw_iter"//str(iter,3)//".dat",wmats,dcmplx(ReSigma(:),ImSigma(:)))
    return
  end subroutine solve_eqs




  function Chi_qv(p,qvec,m) result(chi)
    real(8),dimension(:),intent(in) :: p
    real(8),dimension(:),intent(in) :: qvec
    integer,intent(in)              :: m
    real(8)		   	    :: chi
    real(8)                         :: Tz,dTz
    real(8),dimension(L)            :: ReSigma,ImSigma
    real(8)                         :: kvec(2)
    real(8)                         :: wn,wn_plus_m
    real(8)                         :: Mk,Mk_plus_q
    real(8)                         :: ek,ek_plus_q
    real(8)                         :: xk,xk_plus_q
    real(8)                         :: yk,yk_plus_q
    real(8)                         :: SimEk,SimEk_plus_q
    real(8)                         :: Dk,Dk_plus_q
    real(8)                         :: num,den
    integer                         :: ik,n
    real(8)                         :: kx,ky,qx,qy,vkq
    ReSigma = p(1:L)
    ImSigma = p(L+1:2*L)
    Tz      = p(2*L+1)
    dTz     = p(2*L+2)
    Chi     = 0d0
    qx = qvec(1)
    qy = qvec(2)
    do ik=1,Nktot
       kx        = Kgrid(ik,1)
       ky        = Kgrid(ik,2)     
       ek        = -1d0*(cos(kx)+cos(ky))
       ek_plus_q = -1d0*(cos(kx+qx)+cos(ky+qy))
       xk        = lambda*sin(kx)
       yk        = lambda*sin(ky)
       xk_plus_q = lambda*sin(kx+qx)
       yk_plus_q = lambda*sin(ky+qy)
       vkq       = xk*xk_plus_q + yk*yk_plus_q
       !
       do n=1,Lf
          wn            = pi/beta*(2*n-1)-ImSigma(n)
          wn_plus_m    = pi/beta*(2*(n+m)-1)-ImSigma(n+m)
          !
          Mk            = Mh - Tz*gint/2d0 + ReSigma(n)
          Mk_plus_q     = Mh - Tz*gint/2d0 + ReSigma(n+m)
          !
          simEk         = (Mk + ek)**2 + xk**2 + yk**2
          simEk_plus_q  = (Mk_plus_q + ek)**2 + xk_plus_q**2 + yk_plus_q**2
          !
          Dk            = wn**2d0 + simEk
          Dk_plus_q     = wn_plus_m**2 + simEk_plus_q
          !
          num = wn*wn_plus_m - Mk*Mk_plus_q + vkq
          den = Dk*Dk_plus_q
          Chi = Chi + num/den/beta/Nktot
       enddo
    enddo
  end function Chi_qv







  !For MF calculations:
  function bhz_f(a) result(f)
    real(8),dimension(:) :: a
    real(8)              :: f
    real(8)              :: integral
    tz = a(1)
    call quad2d(Nky,integral)    
    f = gint/2d0*(Tz**2) - 2d0*integral
  end function bhz_f

  subroutine quad2d(Ly,int)
    integer,optional                 :: Ly
    integer                          :: i,Ly_
    real(8)                          :: int,inty
    real(8),dimension(:),allocatable :: Yarray
    Ly_=50;if(present(Ly))Ly_=Ly
    !
    allocate(Yarray(Ly_))
    Yarray = linspace(0d0,pi2,Ly_)
    int=0d0
    do i=1,Ly_
       ky_g = Yarray(i)
       call quad(hk_x,0d0,pi2,result=inty)
       int=int+inty/Ly_/pi2
    enddo
  end subroutine quad2d

  function hk_x(kx) result(fx)
    real(8) :: kx
    real(8) :: fx
    real(8) :: ek,x2,y2
    ek = -1d0*(cos(kx)+cos(ky_g))
    x2  =  lambda*sin(kx)  ;x2=x2**2
    y2  =  lambda*sin(ky_g);y2=y2**2
    fx  = sqrt( (mh - gint*Tz/2d0 + ek)**2 + (x2+y2) )
  end function hk_x

end program bhz_2d


