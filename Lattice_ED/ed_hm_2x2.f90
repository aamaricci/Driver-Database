program hubbard_plaquette
  USE EDLAT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none


  character(len=16)                :: finput
  real(8)                          :: ts,eloc
  integer                          :: Nso,Nk
  real(8),dimension(:),allocatable :: wr
  complex(8),allocatable           :: Gij(:,:,:,:)
  complex(8),allocatable           :: G0ij(:,:,:,:)
  complex(8),allocatable           :: Gk(:,:,:),G0k(:,:,:),Sigk(:,:,:)
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(ts,"TS",finput,default=-0.25d0,comment="hopping parameter")
  call parse_input_variable(eloc,"ELOC",finput,default=0d0,comment="local energies")
  call ed_read_input(str(finput))
  !
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  if(Nspin/=1.OR.Norb/=1.OR.Nsites(1)/=4)stop "This test driver is 2x2 square plaquette only"


  !2d square with OBC and staggered energies
  call ed_Hij_init(Nsites)
  !> ionic potential
  call ed_Hij_add_link(1,1,1,1,1, one*eloc)
  call ed_Hij_add_link(2,2,1,1,1,-one*eloc)
  call ed_Hij_add_link(3,3,1,1,1, one*eloc)
  call ed_Hij_add_link(4,4,1,1,1,-one*eloc)
  !> hoppings
  call ed_Hij_add_link(1,2,1,1,1,one*ts)
  call ed_Hij_add_link(1,4,1,1,1,one*ts)
  call ed_Hij_add_link(2,1,1,1,1,one*ts)
  call ed_Hij_add_link(2,3,1,1,1,one*ts)
  call ed_Hij_add_link(3,2,1,1,1,one*ts)
  call ed_Hij_add_link(3,4,1,1,1,one*ts)
  call ed_Hij_add_link(4,1,1,1,1,one*ts)
  call ed_Hij_add_link(4,3,1,1,1,one*ts)
  !
  call ed_Hij_info()
  call ed_Hij_write()


  call ed_init_solver()
  call ed_solve()

  allocate(wr(Lreal))
  wr = linspace(wini,wfin,Lreal)

  Nso = Nsites(1)
  allocate(Gij(Nspin,Nso,Nso,Lreal))
  call ed_get_gimp_realaxis(Gij)

  allocate(G0ij, mold=Gij)
  call ed_Hij_get_g0func(dcmplx(wr,eps),G0ij)


  Nk = Nso
  allocate(Gk(Nspin,Nk,Lreal))
  allocate(G0k, mold=Gk)
  allocate(Sigk, mold=Gk)
  call get_zeros()

contains




  include 'get_zero.f90'


  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Hk model for the 2d square lattice
  !-------------------------------------------------------------------------------------------
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:) :: kpoint
    integer              :: N,ih
    real(8)              :: kx,ky
    complex(8)           :: hk(N,N)
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = zero
    do ih=1,N
       Hk(ih,ih) = -one*2d0*ts*(cos(kx)+cos(ky))
    enddo
  end function hk_model



  !   !---------------------------------------------------------------------
  !   !PURPOSE: GET ZEROS ON THE REAL AXIS
  !   !---------------------------------------------------------------------
  !   subroutine get_zeros()
  !     integer                                       :: i,j,ik,ix,iy,Nso,Nktot,Npts
  !     integer                                       :: iorb,jorb
  !     integer                                       :: isporb,jsporb
  !     integer                                       :: ispin,jspin
  !     integer                                       :: iso,unit
  !     real(8),dimension(Nspin*Norb)                 :: dzeta
  !     complex(8),allocatable,dimension(:,:,:)       :: Hk_bare,Hk_topo
  !     complex(8),dimension(Nspin*Norb,Nspin*Norb)   :: zeta,z_adj,fg,gdelta,fgk
  !     complex(8),dimension(:,:,:,:,:),allocatable   :: gloc
  !     complex(8),dimension(:,:,:,:,:,:),allocatable :: Sigmareal,Sigmamats
  !     complex(8),dimension(:,:,:,:),allocatable     :: gk,gfoo,ReSmat
  !     complex(8)                                    :: iw
  !     complex(8),dimension(:,:),allocatable         :: detGiw
  !     real(8),dimension(Lreal)                      :: Den
  !     real(8),dimension(:),allocatable              :: Ipoles,Xcsign,Iweight
  !     real(8),dimension(:,:),allocatable            :: Ipoles3d,kpoints
  !     real(8),dimension(:,:),allocatable            :: Mpoles,Mweight
  !     real(8),dimension(:,:,:),allocatable          :: Mpoles3d
  !     integer                                       :: Linterval
  !     integer                                       :: count,Ninterval,maxNinterval,int
  !     real(8)                                       :: sign,sign_old
  !     real(8),dimension(:,:),allocatable :: kpath
  !     !
  !     Nso=Nspin*Norb
  !     !
  !     allocate(kpath(3,2))
  !     kpath(1,:)=[0.0,0.0]!G-e<-R
  !     kpath(2,:)=[1.0,0.0]!G
  !     kpath(3,:)=[2.0,0.0]!G+e->R
  !     !kpath(4,:)=[1.0-0.35,0.0]!G-e<-R
  !     !kpath(5,:)=[1.0,0.0]!G
  !     !kpath(6,:)=[1.0+0.35,0.0]!G+e->R
  !     kpath=kpath*pi
  !     Npts  = size(kpath,1)
  !     Nktot = (Npts-1)*Nkpath
  !     !
  !     if(allocated(Hk_bare))deallocate(Hk_bare)
  !     allocate(Hk_bare(Nspin*Norb,Nspin*Norb,Nktot));Hk_bare=zero
  !     call TB_build_model(hk_bare,hk_periodized,Nspin*Norb,kpath,Nkpath)
  !     !
  !     !
  !     allocate(kpoints(Nktot,2))
  !     call TB_build_kgrid(kpath,Nkpath,kpoints)
  !     if(allocated(Sigmamats))deallocate(Sigmamats)
  !     if(allocated(Sigmareal))deallocate(Sigmareal)
  !     allocate(Sigmamats(Nktot,Nspin,Nspin,Norb,Norb,Lmats))
  !     allocate(Sigmareal(Nktot,Nspin,Nspin,Norb,Norb,Lreal))
  !     do ik=1,Nktot
  !        Sigmamats(ik,:,:,:,:,:)=periodize_sigma_mscheme_mats(kpoints(ik,:),wprint=.false.)
  !        Sigmareal(ik,:,:,:,:,:)=periodize_sigma_mscheme_real(kpoints(ik,:),wprint=.false.)
  !     enddo
  !     !
  !     Linterval = 50000 !Maximum number of allowed intervals to look for zeros&poles
  !     !
  !     allocate(Xcsign(0:Linterval))
  !     allocate(Ipoles(Nktot),Iweight(Nktot))
  !     allocate(Mpoles(Nktot,Linterval),Mweight(Nktot,Linterval))
  !     !
  !     !FINDING THE POLES:
  !     !assume \eps=0.d0 ==> the ImSigma(poles)=0 this condition should be automatically
  !     !verified at the pole from definition of the pole (the ImSigma at the pole is just
  !     !an artificial broadening of an otherwise delta function, whose position should be 
  !     !determined by ReSigma only.
  !     Ipoles=0.d0   
  !     Mpoles=0.d0
  !     write(LOGfile,*)"Solving for the zeros..."
  !     maxNinterval=-1
  !     do ik=1,Nktot
  !        do i=1,Lreal
  !           zeta(:,:) = (wr(i)+xmu)*eye(Nso) - Hk_bare(:,:,ik) - nn2so(Sigmareal(ik,:,:,:,:,i))
  !           call inv(zeta)
  !           Den(i) = dreal(zeta(1,1))*dreal(zeta(2,2)) - dreal(zeta(1,2)*zeta(2,1))
  !        enddo
  !        Xcsign(0)=0.d0
  !        count=0
  !        sign_old=sgn(Den(Lreal/2+1))
  !        do i=Lreal/2+1,Lreal
  !           sign=sgn(Den(i))
  !           if(sign*sign_old<1)then
  !              count=count+1
  !              if(count>Linterval)stop "Allocate Xcsign to a larger array."
  !              Xcsign(count)=wr(i)
  !           endif
  !           sign_old=sign
  !        enddo
  !        Ninterval=count
  !        if(count>maxNinterval)maxNinterval=count
  ! l       call init_finter(finter_func,wr,Den,3)
  !        do int=1,Ninterval
  !           Mpoles(ik,int) = brentq(det_poles,Xcsign(int-1),Xcsign(int))
  !           Mweight(ik,int)= get_weight(hk_bare(:,:,ik)-nn2so(Sigmamats(ik,:,:,:,:,1)))
  !        enddo
  !        ipoles(ik) = brentq(det_poles,0.d0,wr(Lreal))
  !        iweight(ik)= get_weight(hk_bare(:,:,ik)-nn2so(Sigmamats(ik,:,:,:,:,1)))
  !        call delete_finter(finter_func)
  !     enddo
  !     call splot("BHZzeros.ed",ipoles,iweight)
  !     do int=1,maxNinterval
  !        unit=free_unit()
  !        open(unit,file="BHZzeros_int"//reg(txtfy(int))//".ed")
  !        if(any((Mpoles(:,int)/=0.d0)))then
  !           do ik=1,Nktot
  !              if(Mpoles(ik,int)/=0.d0)write(unit,*)ik-1,Mpoles(ik,int),Mweight(ik,int)
  !           enddo
  !        endif
  !        close(unit)
  !     enddo
  !     !
  !   end subroutine get_zeros





end program hubbard_plaquette






