!   Solve the Hubbard model with AFM 2 atoms in the basis 
program ed_hm_square_afm2
  USE DMFT_ED
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
  real(8),allocatable,dimension(:)              :: Wtk
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
  call parse_input_variable(ts     , "TS"     , finput, default=1.d0)
  call parse_input_variable(nkx    , "NKX"    , finput, default=25)
  call parse_input_variable(nkpath , "NKPATH" , finput, default=500)
  call parse_input_variable(wmixing, "WMIXING", finput, default=0.75d0)
  call parse_input_variable(spinsym, "SPINSYM", finput, default=.false.)
  call parse_input_variable(neelsym, "NEELSYM", finput, default=.true.)

  call ed_read_input(trim(finput),comm)


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
  call ed_init_solver(comm,Bath_ineq)
  do ip=1,Nineq
     call ed_break_symmetry_bath(Bath_ineq(ip,:),sb_field,(-1d0)**(ip+1))
  enddo


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !solve the impurity problem:
     call ed_solve(Comm,Bath_ineq,Hloc_ineq)
     !
     !retrieve inequivalent self-energies:
     do ilat=1,Nineq
      call ed_get_sigma_matsubara(Smats_ineq(ilat,:,:,:,:,:),Nineq)
     enddo
     !
     if(neelsym)then
        do ispin=1,2
           Smats(2,ispin,ispin,:,:,:)=Smats(1,3-ispin,3-ispin,:,:,:)
           Sreal(2,ispin,ispin,:,:,:)=Sreal(1,3-ispin,3-ispin,:,:,:)
        enddo
     endif
     !
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
     !
     !fold to the inequivalent sites
     do ip=1,Nineq
        Gmats_ineq(ip,:,:,:,:,:) = Gmats(ip,:,:,:,:,:)
     enddo
     !
     !
     !
     ! compute the Weiss field (only the Nineq ones)
     call dmft_self_consistency(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,cg_scheme)
     !
     !
     ! fit baths and mix result with old baths
     call ed_chi2_fitgf(Comm,Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=1)
     if(spinsym)then
        call ed_spin_symmetrize_bath(Bath_ineq,save=.true.)
     else
        call ed_chi2_fitgf(Comm,Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=2)
     endif
     !
     !
     ! Mixing:
     if(iloop>1)Bath_ineq = wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath_ineq

     ! Convergence
     converged = check_convergence(Weiss_ineq(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
     call Bcast_MPI(Comm,converged)
     !
     call end_loop
  enddo


  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)




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
    allocate(Wtk(Nktot))
    allocate(modelHloc(Nlso,Nlso))
    call TB_build_model(Hk,hk_model,Nlso,[Nkx,Nkx])
    Wtk=1.d0/dble(Nktot)
    modelHloc = sum(Hk(:,:,:),dim=3)/Nktot
    where(abs((modelHloc))<1.d-9)modelHloc=0.d0
    !
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
    hk=zero
    hk(1,2) = -ts*(one+exp(xi*2*kx)+exp(xi*(kx+ky))+exp(xi*(kx-ky)))
    hk(2,1) = -ts*(one+exp(-xi*2*kx)+exp(-xi*(kx+ky))+exp(-xi*(kx-ky)))
  end function hk_model
  

  function lso2nnn_reshape(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function lso2nnn_reshape
  
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


end program ed_hm_square_afm2



