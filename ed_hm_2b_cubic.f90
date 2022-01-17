program ed_hm_2b_cubic
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none  
  integer                                     :: iloop,Nb,Nk,Nso,Nktot
  logical                                     :: converged
  real(8)                                     :: wband,ts,de
  real(8),allocatable                         :: dens(:)
  !
  !Bath:
  real(8),allocatable,dimension(:)            :: Bath,Bath_Prev
  !
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:) :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:) :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Weiss,Weiss_
  real(8)                                     :: crystal_field,var,wmixing
  !
  complex(8),allocatable,dimension(:,:,:)     :: Hk
  complex(8),dimension(:,:,:,:),allocatable   :: Hloc




  
  !Local, user defined variables. Parse_input (in Scifor) add them to the list
  call parse_input_variable(wband,"WBAND","inputED.conf",default=1.d0)
  call parse_input_variable(Nk,"Nk","inputED.conf",default=10)
  call parse_input_variable(crystal_field,"CRYSTAL_FIELD","inputED.conf",default=0.d0)
  call parse_input_variable(wmixing,"WMIXING","inputED.conf",default=1.d0)
  !ED read input
  call ed_read_input("inputED.conf")
  !Add ctrl variables to DMFT_TOOLS memory pool
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  !                               
  !
  Nso = Nspin*Norb
  !
  !
  ts=wband/6.d0

  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));Gmats=0.d0
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));Smats=0.d0
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats));Weiss=0.d0
  allocate(Weiss_(Nspin,Nspin,Norb,Norb,Lmats));Weiss_=0.d0
  !
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal));Greal=0.d0
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal));Sreal=0.d0
  !
  allocate(dens(Norb));

  !+- the non-interacting hamiltonian
  !Here Hk and Hloc are allocated and built
  call build_Hk
  
  !+- Setup the Solver -+!
  !Get the dimension of the bath, user side.
  !this is evaluated internally by code given all the inputs
  Nb=ed_get_bath_dimension()
  !
  allocate(Bath(Nb))            !actual bath
  allocate(Bath_prev(Nb))       !prev bath, for mixing
  call ed_init_solver(bath)


  !+- DMFT loop -+!
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath,Hloc)
     !

     !Retrieve impurity self-energies 
     call ed_get_Sigma_matsubara(Smats)
     call ed_get_Sigma_realaxis(Sreal)

     !
     !Compute the local greens function of the lattice
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     !

     !+- get the new Weiss field by imposing that G_{lattice} = G_{impurity}
     call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,trim(cg_scheme))
     !
     !+- printout the local GF and the Weiss field
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)
     !
     !

     !+- get the new bath parameters
     call ed_chi2_fitgf(Weiss,bath,ispin=1)

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath

     converged = check_convergence(Weiss(1,1,1,1,:)+Weiss(1,1,2,2,:),dmft_error,nsuccess,nloop)

     !
     !+- retrieve the density and check for chemical potential adjustments
     call ed_get_dens(dens)
     if(nread/=0.d0) call ed_search_variable(xmu,dens(1)+dens(2),converged)

  enddo
  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)
  !
  call dmft_kinetic_energy(Hk,Sreal)
  !

contains


  !---------------------------------------------------------------------
  !GET 2bands cubic lattice HAMILTONIAN (from the NonInteracting code)
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j
    !
    !Set the reciprocal lattice basis vector
    call TB_set_bk(bkx=[pi2,0d0,0d0],bky=[0d0,pi2,0d0],bkz=[0d0,0d0,pi2])
    !
    write(LOGfile,*)"Build H(k):"
    Nktot=Nk**3
    write(*,*)"# of k-points     :",Nktot
    write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    if(allocated(Hloc))deallocate(Hloc)
    allocate(Hk(Nso,Nso,Nktot)) ;Hk=zero
    !
    call TB_build_model(Hk,hk_2bands,Nso,[Nk,Nk,Nk])
    ! if(present(file))then
    call TB_write_hk(Hk,"Hk.dat",&
         Nlat=1,&
         Nspin=1,&
         Norb=Norb,&
         Nkvec=[Nk,Nk,Nk])
    ! endif
    allocate(Hloc(Nspin,Nspin,Norb,Norb));  Hloc = 0d0
    Hloc = j2so(sum(Hk,dim=3))/Nktot !send a 4dim array into a 2dim one
    where(abs(dreal(Hloc))<1d-6)Hloc=zero
    ! call TB_write_Hloc(so2j(Hloc))
  end subroutine build_hk


  !--------------------------------------------------------------------!
  !2bands cubic lattice HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_2bands(kvec,N) result(hk)
    integer                   :: N
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: ek,kx,ky,kz
    if(N/=2)stop "hk_2bands error: N != Nspin*Norb == 2"
    kx=kvec(1)
    ky=kvec(2)
    kz=kvec(3)
    ek = -2.d0*ts*(cos(kx)+cos(ky)+cos(kz))
    Hk = ek*eye(N)              !build the non-local part of the H(k)
    Hk(1,1) = Hk(1,1) + crystal_field*0.5d0
    Hk(2,2) = Hk(2,2) - crystal_field*0.5d0
  end function hk_2bands





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


end program ed_hm_2b_cubic
