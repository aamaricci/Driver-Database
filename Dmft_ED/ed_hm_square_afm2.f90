!   Solve the Hubbard model with AFM 2 atoms in the basis 
program ed_hm_square_afm2
  USE EDIPACK2
  USE EDIPACK2INEQ
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none

  integer                                       :: ip,iloop,ilat,ineq,Lk,Nso,Nlso,ispin,iorb
  logical                                       :: converged
  integer                                       :: Nineq,Nlat
  !Bath:
  integer                                       :: Nb
  real(8),allocatable                           :: Bath_ineq(:,:),Bath_prev(:,:)

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal,Sreal_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal,Greal_ineq
  !Hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk ![Nlat*Nspin*Norb,Nlat*Nspin*Norb,Nk]
  complex(8),allocatable,dimension(:,:)         :: modelHloc ![Nlat*Nspin*Norb,Nlat*Nspin*Norb]
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc_ineq
  !variables for the model:
  character(len=16)                             :: finput
  real(8)                                       :: ts,wmixing
  integer                                       :: Nktot,Nkx,Nkpath
  logical                                       :: spinsym,neelsym

  integer                                       :: comm,rank
  logical                                       :: master


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput   , "FINPUT" , default='inputED.conf')
  call parse_input_variable(ts     , "TS"     , finput, default=1d0)
  call parse_input_variable(nkx    , "NKX"    , finput, default=25)
  call parse_input_variable(nkpath , "NKPATH" , finput, default=500)
  call parse_input_variable(wmixing, "WMIXING", finput, default=0.75d0)
  call parse_input_variable(spinsym, "SPINSYM", finput, default=.false.)
  call parse_input_variable(neelsym, "NEELSYM", finput, default=.true.)

  call ed_read_input(trim(finput))


  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")


  Nlat=2
  Nineq=2
  Nso=Nspin*Norb
  Nlso=Nlat*Nso
  if(Norb/=1)stop  "Norb != 1"
  if(Nspin/=2)stop "Nspin != 2"

  if(neelsym)then
     Nineq=1
     write(*,*)"Using Neel symmetry to refold BZ"
     write(*,*)"Using Nineq sites=",Nineq
     write(*,*)"Symmetries used are:"
     write(*,*)"(site=2,l,s)=(site=1,l,-s)"
  endif

  if(spinsym)sb_field=0.d0

  !Allocate Weiss Field:
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  !
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  !
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb))


  call build_hk()
  Hloc = lso2nnn_reshape(modelHloc,Nlat,Nspin,Norb)
  do ip=1,Nineq
     Hloc_ineq(ip,:,:,:,:) = Hloc(ip,:,:,:,:)
  enddo


  !Setup solver
  Nb=ed_get_bath_dimension()
  allocate(Bath_ineq(Nineq,Nb))
  allocate(Bath_prev(Nineq,Nb))
  call ed_init_solver(Bath_ineq)
  call ed_break_symmetry_bath(Bath_ineq,sb_field, (/( (-1d0)**(ip+1), ip=1,Nineq)/) )


  call ed_set_hloc(Hloc_ineq,Nineq)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !solve the impurity problem:
     call ed_solve(Bath_ineq)
     !
     !retrieve inequivalent self-energies:
     call ed_get_sigma(Smats_ineq,Nineq,axis='m')
     call ed_get_sigma(Sreal_ineq,Nineq,axis='r')
     !
     !extend them to the lattice using symmetry if this applies
     do ip=1,Nineq
        Smats(ip,:,:,:,:,:) = Smats_ineq(ip,:,:,:,:,:)
        Sreal(ip,:,:,:,:,:) = Sreal_ineq(ip,:,:,:,:,:)
     enddo
     if(neelsym)then
        do ispin=1,2
           Smats(2,ispin,ispin,:,:,:)=Smats(1,3-ispin,3-ispin,:,:,:)
           Sreal(2,ispin,ispin,:,:,:)=Sreal(1,3-ispin,3-ispin,:,:,:)
        enddo
     endif
     !
     !
     ! compute the local gf:
     call dmft_get_gloc(Hk,Gmats,Smats,axis='m')
     call dmft_write_gf(Gmats,"Gloc",axis='m',iprint=4)
     !
     !fold to the inequivalent sites
     do ip=1,Nineq
        Gmats_ineq(ip,:,:,:,:,:) = Gmats(ip,:,:,:,:,:)
     enddo
     !
     !
     !
     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='delta')then
        call dmft_self_consistency(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq)
     else
        call dmft_self_consistency(Gmats_ineq,Smats_ineq,Weiss_ineq)
     endif
     !
     !
     ! fit baths and mix result with old baths
     call ed_chi2_fitgf(Weiss_ineq,Bath_ineq,ispin=1)
     if(spinsym)then
        call ed_spin_symmetrize_bath(Bath_ineq,save=.true.)
     else
        call ed_chi2_fitgf(Weiss_ineq,Bath_ineq,ispin=2)
     endif
     !
     !
     ! Mixing:
     if(iloop>1)Bath_ineq = wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath_ineq

     ! Convergence
     if(master)converged = check_convergence(Weiss_ineq(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
     call Bcast_MPI(Comm,converged)
     !
     if(master)call end_loop
  enddo
  call ed_finalize_solver()

  call dmft_get_gloc(Hk,Greal,Sreal,axis='r')
  call dmft_write_gf(Greal,"Gloc",axis='r',iprint=4)


  call dmft_kinetic_energy(Hk,Smats)

  call finalize_MPI()



contains



  !--------------------------------------------------------------------!
  !PURPOSE: BUILD THE H(k) FOR THE BHZ-AFM MODEL.
  !--------------------------------------------------------------------!
  subroutine build_hk()
    integer                                 :: Npts
    integer                                 :: i,j,k,ik,iorb,jorb
    integer                                 :: ix,iy,iz
    real(8)                                 :: kx,ky,kz
    real(8),dimension(:),allocatable        :: kxgrid,kygrid
    real(8),dimension(:,:),allocatable      :: kpath
    real(8),dimension(2)                    :: bk1,bk2,kvec
    real(8)                                 :: n(Nlso)
    complex(8)                              :: w
    complex(8)                              :: Gmats(Nlso,Nlso,Lmats),Greal(Nlso,Nlso,Lreal)
    complex(8)                              :: Smats(Nlso,Nlso,Lmats),Sreal(Nlso,Nlso,Lreal)
    !
    Nktot=Nkx*Nkx
    write(LOGfile,*)"Using Nk_total="//txtfy(Nktot)
    !
    !
    !>Reciprocal lattice basis vector  
    bk1=  pi*[ 1d0, -1d0 ]
    bk2=2*pi*[ 0d0,  1d0 ]
    call TB_set_bk(bk1,bk2)
    !
    !
    !>Get TB Hamiltonian matrix
    allocate(Hk(Nlso,Nlso,Nktot))

    call TB_build_model(Hk,hk_model,Nlso,[Nkx,Nkx])

    allocate(modelHloc(Nlso,Nlso))
    modelHloc = sum(Hk(:,:,:),dim=3)/Nktot
    where(abs(dreal(modelHloc))<1.d-9)modelHloc=0.d0
    call TB_write_Hloc(modelHloc)
    !
    !
    !solve along the standard path in the 2D BZ.
    Npts=4
    allocate(kpath(Npts,3))
    kpath(1,:)=kpoint_Gamma
    kpath(2,:)=kpoint_X1
    kpath(3,:)=kpoint_M1
    kpath(4,:)=kpoint_Gamma
    call TB_solve_model(hk_model,Nlso,kpath,Nkpath,&
         colors_name=[red1,blue1, red1,blue1],&
         points_name=[character(len=20) :: "{\Symbol G}","X","M","{\Symbol G}"],&
         file="Eigenbands_afm2.ed")
  end subroutine build_hk




  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    complex(8),dimension(N,N)     :: hk
    complex(8),dimension(N,N)     :: h0,tk
    !
    if(N/=Nlso)stop "hk_model error: N != Nlso" 
    kx = kpoint(1)
    ky = kpoint(2)
    !
    ! Hk =  -t * | 0                  1 + e^ikx(e^ikx + e^iky) |
    !            | 1 + e^-ikx(e^-ikx + e^-iky)   0             |
    !
    hk=zero
    hk(1,2) = -ts*(one+exp(xi*2*kx)+exp(xi*(kx+ky))+exp(xi*(kx-ky)))
    hk(3,4) = -ts*(one+exp(xi*2*kx)+exp(xi*(kx+ky))+exp(xi*(kx-ky)))
    hk(2,1) = -ts*(one+exp(-xi*2*kx)+exp(-xi*(kx+ky))+exp(-xi*(kx-ky)))
    hk(4,3) = -ts*(one+exp(-xi*2*kx)+exp(-xi*(kx+ky))+exp(-xi*(kx-ky)))
  end function hk_model







  !--------------------------------------------------------------------!
  !Reshaping functions:                                                !
  !--------------------------------------------------------------------!

  function nnn2lso_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Fin
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
       Fout(is,js) = Fin(ilat,ispin,jspin,iorb,jorb)
    enddo
  end function nnn2lso_reshape

  function lso2nnn_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Fin
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
       Fout(ilat,ispin,jspin,iorb,jorb) = Fin(is,js)
    enddo
  end function lso2nnn_reshape




end program ed_hm_square_afm2



