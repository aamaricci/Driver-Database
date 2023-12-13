subroutine get_zeros
  integer                                 :: ik,ir,jr,Lk,Lr
  real(8)                                 :: arg
  real(8),dimension(:,:),allocatable      :: kgrid,rgrid
  real(8),dimension(2)                    :: Kvec,Ri,Rj
  real(8),dimension(:,:),allocatable      :: Akw,A0kw
  real(8),dimension(:),allocatable        :: Aloc,A0loc
  complex(8),dimension(:,:,:),allocatable :: Hk
  ! integer                                       :: i,j,ik,ix,iy,Nso,Nktot,Npts
  complex(8),dimension(1,1)   :: zeta
  real(8),dimension(Lreal)                      :: Den
  real(8),dimension(:),allocatable              :: Ipoles,Xcsign,Iweight
  real(8),dimension(:,:),allocatable            :: Mpoles,Mweight
  integer                                       :: Linterval
  integer                                       :: count,Ninterval,maxNinterval,int
  real(8)                                       :: sign,sign_old
  type(finter_type)                                                      :: finter_func


  Lk = Nk
  Lr = Nk
  call TB_set_ei(eix=[1d0,0d0],eiy=[0d0,1d0])
  call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])


  allocate(kgrid(Lk,2))
  call TB_build_kgrid([2,2],Kgrid)
  ! Kgrid(1,:)=[0,0]*pi
  ! Kgrid(2,:)=[1,0]*pi
  ! Kgrid(3,:)=[1,1]*pi
  ! Kgrid(4,:)=[0,1]*pi  
  do ik=1,Lk
     print*,Kgrid(ik,:)
  enddo

  
  allocate(rgrid(Lr,2))
  call TB_build_rgrid([2,2],Rgrid);Rgrid=Rgrid-1d0
  Rgrid(1,:)=[0,0]
  Rgrid(2,:)=[1,0]
  Rgrid(3,:)=[1,1]
  Rgrid(4,:)=[0,1]
  do ir=1,Lr
     print*,Rgrid(ir,:)     
  enddo



  !Build Gk and sum_k Gk
  allocate(Akw(Nk,Lreal),Aloc(Lreal));Aloc=0d0
  allocate(A0kw(Nk,Lreal),A0loc(Lreal));A0loc=0d0
  Gk = zero
  G0k= zero
  do ik=1,Lk
     Kvec = Kgrid(ik,:)
     do ir=1,Lr
        Ri = Rgrid(ir,:)
        do jr=1,Lr
           Rj = Rgrid(jr,:)
           arg= dot_product(Kvec,Ri-Rj)
           Gk(:,ik,:) = Gk(:,ik,:)  + exp(-xi*arg)*Gij(:,ir,jr,:)/Lr
           G0k(:,ik,:)= G0k(:,ik,:) + exp(-xi*arg)*G0ij(:,ir,jr,:)/Lr
        enddo
     enddo
     Sigk(:,ik,:) = one/G0k(:,ik,:) - one/Gk(:,ik,:)
     Akw(ik,:)  = -dimag(Gk(1,ik,:))/pi
     A0kw(ik,:) = -dimag(G0k(1,ik,:))/pi
     Aloc       = Aloc + Akw(ik,:)/Nk
     A0loc      = A0loc + A0kw(ik,:)/Nk
     call splot("Gk"//str(ik)//"_realw.dat",wr,Gk(1,ik,:))
     call splot("A0k"//str(ik)//"_realw.dat",wr,A0kw)
     call splot("Sigma_k"//str(ik)//"_realw.dat",wr,Sigk(1,ik,:))
  enddo
  call splot("Gloc_realw.dat",wr,Aloc)
  call splot("G0loc_realw.dat",wr,A0loc)


  allocate(Hk(1,1,Nk));Hk=zero
  call TB_build_model(hk,hk_model,1,Kgrid)


  Linterval = 50000 !Maximum number of allowed intervals to look for zeros&poles

  allocate(Xcsign(0:Linterval))
  allocate(Ipoles(Lk),Iweight(Lk))
  allocate(Mpoles(Lk,Linterval),Mweight(Lk,Linterval))
  ! !     !



end subroutine get_zeros
