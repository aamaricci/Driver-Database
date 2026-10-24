! Exact single-site BHZ solver using nssDMRG local Hamiltonian and operators.
! Build: make (pkg-config dmrg; ARITH=cmplx or dble must match the library)
program atomic_bhz
  use iso_fortran_env, only: dp=>real64
  use SF_PARSE_INPUT, only: parse_cmd_variable, parse_input_variable
  use SF_LINALG, only: eigh
  use SITES, only: site, electron_site, assignment(=)
  use MATRIX_SPARSE, only: sparse_matrix
  use INPUT_VARS, only: Norb, Nspin, Uloc, Ust, local_jh=>Jh, Jx, Jp, &
       local_hfmode=>HFmode, Dmrg_mode, PBCdmrg, LOGfile
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: nelec=2,nu=41,nm=81,k,a,ispin,iorb,d,g,iu,im,nspec,nground
  real(dp) :: umin=0,umax=8,mmin=0,mmax=4,jh=0.5_dp,jratio=-1
  real(dp) :: jx_factor=0,jp_factor=0,u,m,hund,gap,sz,orb,doub,s2val
  logical :: hfmode=.true.
  character(512) :: input
  type(site) :: dot
  type(sparse_matrix) :: op_sparse
  character(:),allocatable :: key
  real(dp) :: c(16,16,4),num(16,16,4),ntot(16,16),sp(16,16),s2(16,16)
  real(dp) :: tz(16,16),docc(16,16),szop(16,16),h(16,16),v(16,16),e(16)
  integer :: basis(16)

  call parse_cmd_variable(input,'FINPUT',default='inputATOM.conf')
  call parse_input_variable(nelec,'NELEC',input,default=2)
  call parse_input_variable(nu,'NU',input,default=41)
  call parse_input_variable(nm,'NM',input,default=81)
  call parse_input_variable(umin,'UMIN',input,default=0._dp)
  call parse_input_variable(umax,'UMAX',input,default=8._dp)
  call parse_input_variable(mmin,'MMIN',input,default=0._dp)
  call parse_input_variable(mmax,'MMAX',input,default=4._dp)
  call parse_input_variable(jh,'JH',input,default=0.5_dp)
  call parse_input_variable(jratio,'JRATIO',input,default=-1._dp)
  call parse_input_variable(jx_factor,'JX_FACTOR',input,default=0._dp)
  call parse_input_variable(jp_factor,'JP_FACTOR',input,default=0._dp)
  call parse_input_variable(hfmode,'HFMODE',input,default=.true.)
  if(.not.all(ieee_is_finite([umin,umax,mmin,mmax,jh,jratio,jx_factor,jp_factor]))) &
       error stop 'Parameters must be finite'
  if(nelec<0.or.nelec>4.or.nu<1.or.nm<1)error stop 'Require 0<=nelec<=4 and nu,nm>=1'
  if(umin<0.or.umax<umin.or.mmax<mmin.or.jh<0.or.jx_factor<0.or.jp_factor<0) &
       error stop 'Invalid parameter range or negative Hund factors'

  ! Initialize only the library state needed by electron_site; no DMRG run.
  Norb=2; Nspin=2; Dmrg_mode='normal'; PBCdmrg=.false.
  allocate(Uloc(Norb))
  open(newunit=LOGfile,status='scratch',action='readwrite')
  call make_dot(0._dp,0._dp,0._dp,0._dp,0._dp,dot)
  do ispin=1,Nspin
     do iorb=1,Norb
        a=iorb+(ispin-1)*Norb
        key=dot%name()//dot%okey(iorb,ispin,ilink='n')
        op_sparse=dot%operators%op(key=key)
        c(:,:,a)=real_operator(op_sparse)
     enddo
  enddo
  call op_sparse%free()
  call dot%free()
  do a=1,4
     num(:,:,a)=matmul(transpose(c(:,:,a)),c(:,:,a))
  enddo
  szop=(num(:,:,1)+num(:,:,2)-num(:,:,3)-num(:,:,4))/2
  sp=matmul(transpose(c(:,:,1)),c(:,:,3))+matmul(transpose(c(:,:,2)),c(:,:,4))
  s2=matmul(szop,szop)+(matmul(sp,transpose(sp))+matmul(transpose(sp),sp))/2
  tz=(num(:,:,1)+num(:,:,3)-num(:,:,2)-num(:,:,4))/2
  docc=matmul(num(:,:,1),num(:,:,3))+matmul(num(:,:,2),num(:,:,4))
  ntot=sum(num,dim=3)
  call select_sector(nelec)
  call self_test()
  open(newunit=nground,file='atomic_ground.dat',status='replace')
  open(newunit=nspec,file='atomic_spectrum.dat',status='replace')
  write(nground,'(a)') '# Equal-weight average over the degenerate ground manifold; Tz=(n1-n2)/2'
  write(nground,'(a)') '# U M J E0 degeneracy gap_above_manifold S2 Sz Tz D n1 n2'
  write(nspec,'(a)') '# U M J level energy'
  do iu=1,nu
     u=umin
     if(nu>1)u=umin+(umax-umin)*real(iu-1,dp)/real(nu-1,dp)
     hund=jh
     if(jratio>=0)hund=jratio*u
     do im=1,nm
        m=mmin
        if(nm>1)m=mmin+(mmax-mmin)*real(im-1,dp)/real(nm-1,dp)
        call build(u,m,hund,jx_factor*hund,jp_factor*hund,h)
        call diagonalize(h,e,v)
        g=count(abs(e(1:d)-e(1))<1.e-10_dp*max(1._dp,maxval(abs(e(1:d)))))
        gap=0
        if(g<d)gap=e(g+1)-e(1)
        s2val=average(s2,v,g); sz=average(szop,v,g)
        orb=average(tz,v,g); doub=average(docc,v,g)
        write(nground,'(4(es24.15,1x),i3,1x,7(es24.15,1x))') &
             u,m,hund,e(1),g,gap,s2val,sz,orb,doub,nelec/2._dp+orb,nelec/2._dp-orb
        do k=1,d
           write(nspec,'(3(es24.15,1x),i3,1x,es24.15)')u,m,hund,k,e(k)
        enddo
     enddo
     write(nground,*)
  enddo
  close(nground); close(nspec); close(LOGfile)
  deallocate(Uloc)
  print '(a)', 'Analytic spectrum checks passed. Wrote atomic_ground.dat and atomic_spectrum.dat'
contains
  ! electron_site owns all interaction signs, couplings and HFMODE shifts.
  subroutine make_dot(u,m,hund,exchange,pair,dot)
    real(dp),intent(in) :: u,m,hund,exchange,pair
    type(site),intent(out) :: dot
#ifdef _CMPLX
    complex(dp) :: hloc(4,4)
#else
    real(dp) :: hloc(4,4)
#endif
    integer :: a
    Uloc=u; Ust=u-2*hund; local_jh=hund; Jx=exchange; Jp=pair
    local_hfmode=hfmode
    hloc=0
    do a=1,4
       hloc(a,a)=m*(-1)**(mod(a-1,2))
    enddo
    rewind(LOGfile)
    dot=electron_site(hloc)
  end subroutine

  subroutine build(u,m,hund,exchange,pair,h)
    real(dp),intent(in) :: u,m,hund,exchange,pair
    real(dp),intent(out) :: h(16,16)
    type(site) :: dot
    type(sparse_matrix) :: h_sparse
    call make_dot(u,m,hund,exchange,pair,dot)
    h_sparse=dot%operators%op(key='H')
    h=real_operator(h_sparse)
    call h_sparse%free()
    call dot%free()
  end subroutine

  ! This atomic Hamiltonian and its C operators are real in both builds.
  function real_operator(op) result(matrix)
    type(sparse_matrix),intent(in) :: op
    real(dp) :: matrix(16,16)
    complex(dp),allocatable :: dense(:,:)
    dense=cmplx(op%as_matrix(),kind=dp)
    if(any(shape(dense)/=[16,16]))error stop 'Expected a two-orbital local operator'
    if(maxval(abs(aimag(dense)))>1.e-12_dp)error stop 'Unexpected complex local operator'
    matrix=real(dense,dp)
  end function

  subroutine select_sector(n)
    integer,intent(in) :: n
    real(dp) :: diagonal(16,16)
    integer :: i
    diagonal=0; d=0
    do i=1,16
       diagonal(i,i)=ntot(i,i)
       if(abs(ntot(i,i)-n)>1.e-12_dp)cycle
       d=d+1; basis(d)=i
    enddo
    if(maxval(abs(ntot-diagonal))>1.e-12_dp)error stop 'N must be diagonal in the local basis'
    if(d==0)error stop 'Empty particle sector'
  end subroutine

  subroutine diagonalize(h,e,v)
    real(dp),intent(in) :: h(16,16)
    real(dp),intent(out) :: e(16),v(16,16)
    real(dp) :: block(d,d)
    integer :: i,j
    if(maxval(abs(h-transpose(h)))>1.e-12_dp)error stop 'Non-Hermitian Hamiltonian'
    if(maxval(abs(matmul(h,ntot)-matmul(ntot,h)))>1.e-10_dp) &
         error stop 'Local Hamiltonian does not conserve N'
    block=0; e=0; v=0
    do j=1,d
       do i=1,d
          block(i,j)=h(basis(i),basis(j))
       enddo
    enddo
    call eigh(block,e(1:d))
    do i=1,d
       v(basis(i),1:d)=block(i,1:d)
    enddo
  end subroutine

  function average(op,v,g)result(value)
    real(dp),intent(in) :: op(16,16),v(16,16)
    integer,intent(in) :: g
    real(dp) :: value
    integer :: k
    value=0
    do k=1,g
       value=value+dot_product(v(:,k),matmul(op,v(:,k)))
    enddo
    value=value/g
  end function

  subroutine self_test()
    real(dp) :: test_h(16,16),test_e(16),test_v(16,16),expected(6),tmp
    integer :: saved_d,saved_basis(16),i,j,l
    logical :: saved_hfmode
    saved_d=d; saved_basis=basis
    saved_hfmode=hfmode; hfmode=.false.
    call select_sector(2)
    do j=0,1
       call build(4._dp,0.7_dp,0.6_dp,0.6_dp*j,0.6_dp*j,test_h)
       call diagonalize(test_h,test_e,test_v)
       expected=[2.2_dp,2.2_dp,2.8_dp-0.6_dp*j,2.8_dp+0.6_dp*j, &
            4-sqrt(1.4_dp**2+(0.6_dp*j)**2),4+sqrt(1.4_dp**2+(0.6_dp*j)**2)]
       do i=2,6
          kloop: do l=i,2,-1
             if(expected(l)>=expected(l-1))exit kloop
             tmp=expected(l); expected(l)=expected(l-1); expected(l-1)=tmp
          enddo kloop
       enddo
       if(maxval(abs(test_e(1:6)-expected))>1.e-12_dp)error stop 'Analytic spectrum check failed'
       if(abs(average(s2,test_v,1)-2)>1.e-12_dp)error stop 'Triplet spin check failed'
    enddo
    d=saved_d; basis=saved_basis; hfmode=saved_hfmode
  end subroutine
end program
