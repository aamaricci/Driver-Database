program ss_bethe
  USE SLAVE_SPINS
  !
  USE SCIFOR
  USE DMFT_TOOLS
  !
  implicit none
  complex(8),allocatable,dimension(:,:,:) :: Hk
  complex(8),allocatable,dimension(:,:)   :: Hloc
  real(8)                                 :: ts(5),v
  real(8)                                 :: Mh(5)
  integer                                 :: Nkx,Nktot,Nlso,ilat,Npts,Nkpath,ik
  real(8),dimension(:,:),allocatable      :: kpath
  real(8),dimension(:),allocatable        :: Zeta,Self

  call parse_input_variable(ts,"ts","inputSS.conf",default=[1d0,1d0,1d0,1d0,1d0])
  call parse_input_variable(v,"v","inputSS.conf",default=0d0)
  call parse_input_variable(Mh,"Mh","inputSS.conf",default=[0d0,0d0,0d0,0d0,0d0])
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

  if(Nspin/=1.OR.Norb>5.OR.Nlat>5)stop "Wrong setup from input file: Nspin!=1 OR Norb>5 OR Nlat>5"
  Nlso=Nlat*Nspin*Norb

  Nktot = Nkx**2
  ! !solve along a path in the 3D BZ.
  Npts = 4
  allocate(kpath(Npts,2))
  kpath(1,:)=[0d0,0d0]
  kpath(2,:)=[1d0,1d0]
  kpath(3,:)=[1d0,0d0]
  kpath(4,:)=[0d0,0d0]
  kpath = kpath*pi


  allocate(Hk(Nlso,Nlso,Nktot),Hloc(Nlso,Nlso))
  call TB_set_bk([pi2,0d0],[0d0,pi2])
  call TB_build_model(Hk,hk_model,Nlso,[Nkx,Nkx])
  Hloc = sum(Hk,dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero 
  call TB_Solve_model(hk_model,Nlso,kpath,Nkpath,&
       colors_name=[black,red,blue,green,magenta],&
       points_name=[character(len=40) ::'G', 'M', 'X', 'G'],&
       file="Bands_2d",iproject=.false.)


  call ss_solve(Hk,ineq_sites=(/(1,ilat=1,Nlat)/) )


  !Solve for the renormalized bands:
  call TB_Solve_model(ss_hk_model,Nlso,kpath,Nkpath,&
       colors_name=[black,red,blue,green,magenta],&
       points_name=[character(len=40) ::'G', 'M', 'X', 'G'],&
       file="zBands_2d",iproject=.false.)


contains

  function hk_model(kvec,N) result(hk)
    real(8),dimension(:)      :: kvec
    integer                   :: N,i
    complex(8),dimension(N,N) :: diagZ
    complex(8),dimension(N,N) :: Hk,H
    real(8)                   :: kx,ky,cx,cy
    kx = kvec(1)
    ky = kvec(2)
    cx = cos(kx)
    cy = cos(ky)
    hk = zero
    do i=1,N
       hk(i,i) = Mh(i)-2d0*ts(i)*(cx+cy)
    enddo
    if(Norb==2)then
       hk(1,2) = v
       hk(2,1) = v
    endif
  end function hk_model


  function ss_Hk_model(kvec,N) result(Hk)
    real(8),dimension(:)      :: kvec
    integer                   :: N
    complex(8),dimension(N,N) :: h,hk
    h = Hk_model(kvec,N)
    call ss_get_ssHk(h,Hloc)
    hk = h
  end function ss_Hk_model

end program ss_bethe
