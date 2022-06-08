program ed_kanemele
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none

  integer                                       :: iloop,Lk,Nso,Nlso,Nlat,Nineq
  logical                                       :: converged
  integer                                       :: ispin,ilat,i!,j
  !
  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:,:)            :: Bath,Bath_prev
  !
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Greal
  complex(8),allocatable,dimension(:,:,:,:)     :: Sekin
  !
  !Hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk
  complex(8),allocatable,dimension(:,:)         :: kmHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  !
  integer,allocatable,dimension(:)              :: ik2ix,ik2iy
  real(8),dimension(2)                          :: e1,e2   !real-space lattice basis
  real(8),dimension(2)                          :: bk1,bk2 !reciprocal space lattice basis
  real(8),dimension(2)                          :: d1,d2,d3
  real(8),dimension(2)                          :: a1,a2,a3
  real(8),dimension(2)                          :: bklen
  !
  !Variables for the model:
  integer                                       :: Nk,Nkpath
  real(8)                                       :: t1,t2,phi,Mh,wmixing
  character(len=32)                             :: finput
  character(len=32)                             :: hkfile
  real(8),allocatable,dimension(:)              :: dens
  !
  !Flags and options
  character(len=32)                             :: bathspins
  logical                                       :: neelsym,xkick,ykick,getbands
  !
  !Replica Hamiltonian
  real(8),dimension(:,:),allocatable            :: lambdasym_vector ![Nlat,:]
  complex(8),dimension(:,:,:,:,:),allocatable   :: Hsym_basis
  !
  !MPI
  integer                                       :: comm,rank
  logical                                       :: master


  !MPI INIT:
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  !Parse additional variables && read Input && read H(k)^2x2
  call parse_cmd_variable(finput,"FINPUT",default='inputKANEMELE.conf')
  !
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in",&
       comment='Hk will be written here')
  call parse_input_variable(nk,"NK",finput,default=100,&
       comment='Number of kpoints per direction')
  call parse_input_variable(nkpath,"NKPATH",finput,default=500,&
       comment='Number of kpoints per interval on kpath. Relevant only if GETBANDS=T.')
  call parse_input_variable(t1,"T1",finput,default=1d0,&
       comment='NN hopping, fixes noninteracting bandwidth')
  call parse_input_variable(t2,"T2",finput,default=0.1d0,&
       comment='Haldane-like NNN hopping-strenght, corresponds to lambda_SO in KM notation')
  call parse_input_variable(phi,"PHI",finput,default=pi/2d0,&
       comment='Haldane-like flux for the SOI term, KM model corresponds to a pi/2 flux')
  call parse_input_variable(mh,"MH",finput,default=0d0,&
       comment='On-site staggering, aka Semenoff-Mass term')
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0,&
       comment='Mixing parameter: 0 means 100% of the old bath (no update at all), 1 means 100% of the new bath (pure update)')
  call parse_input_variable(bathspins,"BathSpins",finput,default="x",&
       comment='x; xy; xz; xyz. Meaning the replica bath will have Sx; Sx,Sy; Sx,Sz; Sx,Sy,Sz components.')
  call parse_input_variable(neelsym,"NEELSYM",finput,default=.true.,&
       comment='If T AFM(xy) symmetry is enforced on the self energies at each loop')
  call parse_input_variable(xkick,"xKICK",finput,default=.true.,&
       comment='If T the bath spins get an initial AFM(x) distortion')
  call parse_input_variable(ykick,"yKICK",finput,default=.false.,&
       comment='If T the bath spins get an initial AFM(y) distortion')
  call parse_input_variable(getbands,"GETBANDS",finput,default=.false.,&
       comment='If T the noninteracting model is solved and the bandstructure stored')
  !
  call ed_read_input(trim(finput),comm)
  !
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")


  !SOME PRELIMINARY CHECKS FOR THIS DRIVER:
  if(.not.(bath_type=="replica".AND.ed_mode=='nonsu2'))&
       stop "Wrong setup from input file: AFMxy requires NONSU2-mode and REPLICA-bath"
  if(BathSpins=="xz".AND.yKICK)&
       stop "Wrong setup from input file: AFM(y) sb-field is not allowed with a xz-only bath"
  if(BathSpins=="x".AND.yKICK)&
       stop "Wrong setup from input file: AFM(y) sb-field is not allowed with a x-only bath"
  if(Norb/=1.OR.Nspin/=2)&
       stop "Wrong setup from input file: Norb=1 AND Nspin=2 is the correct configuration."
  Nlat=2
  Nso=Nspin*Norb
  Nlso=Nlat*Nso                 !=4 = 2(ineq sites)*2(spin)*1(orb)


  !SETUP LATTICE AND H(k)
  !Lattice basis (a=1; a0=sqrt3*a) is:
  !e_1 = a0 [ sqrt3/2 , 1/2 ] = 3/2a[1, 1/sqrt3]
  !e_2 = a0 [ sqrt3/2 ,-1/2 ] = 3/2a[1,-1/sqrt3]
  e1 = 3d0/2d0*[1d0, 1d0/sqrt(3d0)]
  e2 = 3d0/2d0*[1d0,-1d0/sqrt(3d0)]

  !lattice basis: nearest neighbor: A-->B, B-->A
  d1= [  1d0/2d0 , sqrt(3d0)/2d0 ]
  d2= [  1d0/2d0 ,-sqrt(3d0)/2d0 ]
  d3= [ -1d0     , 0d0           ]

  !next nearest-neighbor displacements: A-->A, B-->B, cell basis
  a1 = d1-d3                    !3/2*a[1,1/sqrt3]
  a2 = d2-d3                    !3/2*a[1,-1/sqrt3]
  a3 = d1-d2

  !reciprocal lattice vectors:
  bklen=2d0*pi/3d0
  bk1=bklen*[ 1d0, sqrt(3d0)] 
  bk2=bklen*[ 1d0,-sqrt(3d0)]
  call TB_set_bk(bkx=bk1,bky=bk2)

  !Build the Hamiltonian on a grid or on path
  call build_hk(trim(hkfile),getbands)
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb));Hloc=zero
  Hloc = lso2nnn_reshape(kmHloc,Nlat,Nspin,Norb)


  !ALLOCATE LOCAL FIELDS:
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal=zero


  !SETUP HREPLICA SYMMETRIES: 
  select case(trim(BathSpins))

  case default
     stop "BathSpins not in [x; xy; xz; xyz]"

  case("x")  !Only X spin component in the bath
     allocate(lambdasym_vector(Nlat,2))
     allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,2))
     Hsym_basis(:,:,:,:,1)=so2nn_reshape(pauli_sigma_0,Nspin,Norb)
     Hsym_basis(:,:,:,:,2)=so2nn_reshape(pauli_sigma_x,Nspin,Norb)
     lambdasym_vector(1,:)=[0d0, 0d0]
     lambdasym_vector(2,:)=[0d0, 0d0]

  case("xy")  !Only XY spin components in the bath
     allocate(lambdasym_vector(Nlat,3))
     allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
     Hsym_basis(:,:,:,:,1)=so2nn_reshape(pauli_sigma_0,Nspin,Norb)
     Hsym_basis(:,:,:,:,2)=so2nn_reshape(pauli_sigma_x,Nspin,Norb)
     Hsym_basis(:,:,:,:,3)=so2nn_reshape(pauli_sigma_y,Nspin,Norb)
     lambdasym_vector(1,:)=[0d0, 0d0, 0d0]
     lambdasym_vector(2,:)=[0d0, 0d0, 0d0]

  case("xz")  !Only XZ spin components in the bath
     allocate(lambdasym_vector(Nlat,3))
     allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
     Hsym_basis(:,:,:,:,1)=so2nn_reshape(pauli_sigma_0,Nspin,Norb)
     Hsym_basis(:,:,:,:,2)=so2nn_reshape(pauli_sigma_x,Nspin,Norb)
     Hsym_basis(:,:,:,:,3)=so2nn_reshape(pauli_sigma_z,Nspin,Norb)
     lambdasym_vector(1,:)=[0d0, 0d0, 0d0]
     lambdasym_vector(2,:)=[0d0, 0d0, 0d0]

  case("xyz") !Full XYZ spin freedom in the bath
     allocate(lambdasym_vector(Nlat,4))
     allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,4))
     Hsym_basis(:,:,:,:,1)=so2nn_reshape(pauli_sigma_0,Nspin,Norb)
     Hsym_basis(:,:,:,:,2)=so2nn_reshape(pauli_sigma_x,Nspin,Norb)
     Hsym_basis(:,:,:,:,3)=so2nn_reshape(pauli_sigma_y,Nspin,Norb)
     Hsym_basis(:,:,:,:,4)=so2nn_reshape(pauli_sigma_z,Nspin,Norb);
     lambdasym_vector(1,:)=[0d0, 0d0, 0d0, 0d0]
     lambdasym_vector(2,:)=[0d0, 0d0, 0d0, 0d0]

  end select

  !SETUP SYMMETRY BREAKING KICKS
  if(xKICK)then
     lambdasym_vector(1,2)= +sb_field
     lambdasym_vector(2,2)= -sb_field
     !For the log file
     if(master)write(*,*) "*************************************************"
     if(master)write(*,*) "*                                               *"
     if(master)write(*,*) "*  !Applying an AFMx kick to the initial bath!  *"
     if(master)write(*,*) "*                                               *"
     if(master)write(*,*) "*************************************************"
  endif
  if(yKICK)then  !Safe: look at the preliminary checks
     lambdasym_vector(1,3)= +sb_field
     lambdasym_vector(2,3)= -sb_field
     !For the log file
     if(master)write(*,*) "*************************************************"
     if(master)write(*,*) "*                                               *"
     if(master)write(*,*) "*  !Applying an AFMy kick to the initial bath!  *"
     if(master)write(*,*) "*                                               *"
     if(master)write(*,*) "*************************************************"
  endif

  !SETUP H_replica
  call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
  !this is now elevated to RDMFT: ineq sites (1,2) for the lambdas

  !SETUP SOLVER
  Nb=ed_get_bath_dimension(Hsym_basis)
  allocate(Bath(Nlat,Nb))
  allocate(Bath_prev(Nlat,Nb))
  call ed_init_solver(comm,Bath)



  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     if(neelsym)then
        !solve just one sublattice and get the other by Neel symmetry (xy version)
        call ed_solve(comm,Bath(1,:),Hloc(1,:,:,:,:))
        call ed_get_sigma_matsubara(Smats(1,:,:,:,:,:))
        call ed_get_sigma_realaxis(Sreal(1,:,:,:,:,:))
        Smats(2,1,1,:,:,:) = Smats(1,1,1,:,:,:)   !S(iw)_{B,up,up} = S(iw)_{A,up,up}
        Smats(2,2,2,:,:,:) = Smats(1,2,2,:,:,:)   !S(iw)_{B,dw,dw} = S(iw)_{A,dw,dw}
        Smats(2,1,2,:,:,:) = -Smats(1,1,2,:,:,:)  !S(iw)_{B,up,dw} = -S(iw)_{A,up,dw}
        Smats(2,2,1,:,:,:) = -Smats(1,2,1,:,:,:)  !S(iw)_{B,dw,up} = -S(iw)_{A,dw,up}
        Sreal(2,1,1,:,:,:) = Sreal(1,1,1,:,:,:)   !S(w)_{B,up,up}  = S(w)_{A,up,up}
        Sreal(2,2,2,:,:,:) = Sreal(1,2,2,:,:,:)   !S(w)_{B,dw,dw}  = S(w)_{A,dw,dw}
        Sreal(2,1,2,:,:,:) = -Sreal(1,1,2,:,:,:)  !S(w)_{B,up,dw}  = -S(w)_{A,up,dw}
        Sreal(2,2,1,:,:,:) = -Sreal(1,2,1,:,:,:)  !S(w)_{B,dw,up}  = -S(w)_{A,dw,up}
        if(master)write(*,*) "***********************************"
        if(master)write(*,*) "*                                 *"
        if(master)write(*,*) "*  !Enforcing NEEL(xy) symmetry!  *"
        if(master)write(*,*) "*                                 *"
        if(master)write(*,*) "***********************************"
     else
        !solve both sublattices independently with the RDMFT wrapper: 
        !mpi_lanc=T => MPI lanczos, mpi_lanc=F => MPI for ineq sites
        !Hloc is now mandatory here
        call ed_solve(comm,Bath,Hloc,mpi_lanc=.true.)
        !retrieve all self-energies:
        call ed_get_sigma_matsubara(Smats,Nlat)
        call ed_get_sigma_realaxis(Sreal,Nlat)
        !
     endif
     !
     !COMPUTE THE LOCAL GF:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
     !
     !COMPUTE THE WEISS FIELD (only the Nineq ones)
     call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,cg_scheme)
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=4)
     !
     !FIT THE NEW BATH, starting from the old bath + the supplied delta
     !IF(NONSU2): Sz-conservation is broken -> Allows for magXY
     call ed_chi2_fitgf(comm,Bath,Weiss,Hloc) !Hloc mandatory here, it sets impHloc
     !
     !MIXING:
     if(iloop>1)Bath=wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath
     !
     !CHECK CONVERGENCE. This is now entirely MPI-aware:
     converged = check_convergence(Weiss(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
     !
     call end_loop
  enddo

  !Extract and print retarded self-energy and Green's function 
  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Greal",iprint=4)

  !Compute Kinetic Energy
  call dmft_kinetic_energy(Hk,Smats)


  allocate(Sekin(Nlat,Nspin*Norb,Nspin*Norb,Lmats))
  do i=1,Lmats
     Sekin(:,:,:,i) = nnn2lnn_reshape(Smats(:,:,:,:,:,i),Nlat,Nspin,Norb)
  enddo
  call dmft_kinetic_energy_exp(Hk,Sekin)

  if(master) then
     write(*,*) "!***************************!"
     write(*,*) "!*                         *!"
     write(*,*) "!*   !!!  FINISHED  !!!    *!"
     write(*,*) "!*                         *!"
     write(*,*) "!***************************!"
  endif

  call finalize_MPI()


contains



  !---------------------------------------------------------------------
  !PURPOSE: Get Kane Mele Model Hamiltonian
  !---------------------------------------------------------------------
  subroutine build_hk(file,getbands)
    character(len=*),optional          :: file
    integer                            :: i,j,ik
    integer                            :: ix,iy
    real(8)                            :: kx,ky
    real(8),dimension(2)               :: pointK,pointKp
    integer                            :: iorb,jorb
    integer                            :: isporb,jsporb
    integer                            :: ispin,jspin
    integer                            :: unit
    real(8),dimension(:,:),allocatable :: KPath
    logical                            :: getbands
    !
    Lk= Nk*Nk
    write(LOGfile,*)"Build H(k) Kane-Mele:",Lk
    write(LOGfile,*)"# of SO-bands     :",Nlso
    !
    if(allocated(Hk))deallocate(Hk)
    !
    allocate(Hk(Nlso,Nlso,Lk));Hk=zero
    !
    !
    call TB_build_model(Hk,hk_kanemele_model,Nlso,[Nk,Nk],wdos=.false.)
    !
    !
    allocate(kmHloc(Nlso,Nlso))
    kmHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(kmHloc))<1.d-4)kmHloc=0d0
    if(master)call TB_write_Hloc(kmHloc)
    if(master)call TB_write_Hloc(kmHloc,'Hloc.txt')
    !
    !
    if(getbands)then
       pointK = [2*pi/3, 2*pi/3/sqrt(3d0)]
       pointKp= [2*pi/3,-2*pi/3/sqrt(3d0)]
       if(master)write(*,*) "***************************************"
       if(master)write(*,*) "*                                     *"
       if(master)write(*,*) "*  !Solving noninteracting TB model!  *"
       if(master)write(*,*) "*                                     *"
       if(master)write(*,*) "***************************************"
       if(master)then
          allocate(Kpath(4,2))
          KPath(1,:)=[0,0]
          KPath(2,:)=pointK
          Kpath(3,:)=pointKp
          KPath(4,:)=[0d0,0d0]
          call TB_Solve_model(hk_kanemele_model,Nlso,KPath,Nkpath,&
               colors_name=[red1,blue1,red1,blue1],&
               points_name=[character(len=10) :: "G","K","K`","G"],&
               file="Eigenbands.nint",iproject=.false.)
       endif
    endif
    !
  end subroutine build_hk



  !--------------------------------------------------------------------!
  !Kane-Mele HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_kanemele_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nlso
    complex(8),dimension(2,2)       :: hk11,hk22
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                         :: h0,hx,hy,hz
    !----- Old Version ---------
    !real(8)                         :: kdotd(3),kdota(3)
    !!(k.d_j)
    !kdotd(1) = dot_product(kpoint,d1)
    !kdotd(2) = dot_product(kpoint,d2)
    !kdotd(3) = dot_product(kpoint,d3)
    !!(k.a_j)
    !kdota(1) = dot_product(kpoint,a1)
    !kdota(2) = dot_product(kpoint,a2)
    !kdota(3) = dot_product(kpoint,a3)
    !
    !h0 = 2*t2*cos(phi)*sum( cos(kdota(:)) )
    !hx =-t1*sum( cos(kdotd(:)) )
    !hy =-t1*sum( sin(kdotd(:)) )
    !hz = 2*t2*sin(phi)*sum( sin(kdota(:)) )
    !
    !------- New Version [Consistent with Bernevig&Hughes book] -----------
    ! This way the Hamiltonian will have a Bloch structure so we can obtain 
    ! local (i.e. in the unit-cell) quantities by averaging over all k-points. 
    real(8)           :: kdote1, kdote2
    !
    kdote1 = dot_product(kpoint,e1)
    kdote2 = dot_product(kpoint,e2)
    !
    h0 = 2*t2*cos(phi)*( cos(kdote1) + cos(kdote2) + cos(kdote1-kdote2) )
    hx = t1*( cos(kdote1) + cos(kdote2) + 1)
    hy = t1*( sin(kdote1) + sin(kdote2) )
    hz = 2*t2*sin(phi)*( sin(kdote1) - sin(kdote2) - sin(kdote1-kdote2) )
    !
    hk11 = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z + Mh*pauli_z
    !
    hk22 = h0*pauli_0 + hx*pauli_x + hy*pauli_y - hz*pauli_z + Mh*pauli_z
    !
    hk          = zero
    !hk(1:2,1:2) = hk11
    !hk(3:4,3:4) = hk22
    !create hk such that it works well with lso2nnn_reshape(kmHloc,Nlat,Nspin,Norb)
    hk(1:3:2,1:3:2) = hk11
    hk(2:4:2,2:4:2) = hk22
    !
  end function hk_kanemele_model





  !--------------------------------------------------------------------!
  !Reshaping functions:                                                !
  !--------------------------------------------------------------------!
  function nnn2lnn_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
    integer                                          :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Fin
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb) :: Fout
    integer                                          :: iorb,ispin,ilat,is
    integer                                          :: jorb,jspin,js
    Fout=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb
                   js = jorb + (jspin-1)*Norb
                   Fout(ilat,is,js) = Fin(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function nnn2lnn_reshape

  function lnn2nnn_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb) :: Fin
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb
                   js = jorb + (jspin-1)*Norb
                   Fout(ilat,ispin,jspin,iorb,jorb) = Fin(ilat,is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function lnn2nnn_reshape



  function nnn2lso_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Fin
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Fout(is,js) = Fin(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function nnn2lso_reshape

  function lso2nnn_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Fin
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Fout(ilat,ispin,jspin,iorb,jorb) = Fin(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function lso2nnn_reshape

  function so2nn_reshape(Fin,Nspin,Norb) result(Fout)
    integer                                               :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb)           :: Fin
    complex(8),dimension(Nspin,Nspin,Norb,Norb)           :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Fout(ispin,jspin,iorb,jorb) = Fin(is,js)
             enddo
          enddo
       enddo
    enddo
  end function so2nn_reshape

  function nn2so_reshape(Fin,Nspin,Norb) result(Fout)
    integer                                               :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb)           :: Fin
    complex(8),dimension(Nspin*Norb,Nspin*Norb)           :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Fout(is,js) = Fin(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function nn2so_reshape











  subroutine dmft_kinetic_energy_exp(Hk,Sigma)
    complex(8),dimension(:,:,:)                                     :: Hk ! [Nlat*Nspin*Norb][Nlat*Nspin*Norb][Lk]
    complex(8),dimension(:,:,:,:)                                   :: Sigma ! [Nlat][Nspin*Norb][Nspin*Norb][L]
    !
    integer                                                         :: Lk,Nlso,Liw,Nso,Nlat
    integer                                                         :: ik
    integer                                                         :: i,iorb,ilat,ispin,io,is
    integer                                                         :: j,jorb,jlat,jspin,jo,js
    !
    integer                                                         :: mpi_ierr
    integer                                                         :: mpi_rank
    integer                                                         :: mpi_size
    logical                                                         :: mpi_master
    !
    ! integer                                                       :: Norb,Nporb
    ! integer                                                       :: Nspin  
    ! real(8)                                                       :: beta
    ! real(8)                                                       :: xmu
    real(8),dimension(:),allocatable                                :: wm
    !
    complex(8),dimension(size(Sigma,1),size(Sigma,2),size(Sigma,3)) :: Sigma_HF ![Nlat][Nso][Nso]
    complex(8),dimension(size(Hk,1),size(Hk,2))                     :: Ak,Bk,Ck,Hloc,Hloc_tmp
    complex(8),dimension(size(Hk,1),size(Hk,2))                     :: Gk,Tk
    complex(8),dimension(size(Hk,1),size(Hk,2))                     :: GkTmp,HkTmp
    real(8),dimension(size(Hk,1))                                   :: Nk
    complex(8),dimension(size(Hk,1),size(Hk,2))                     :: Evec
    real(8),dimension(size(Hk,1))                                   :: Eval,Coef
    real(8)                                                         :: spin_degeneracy
    !
    real(8),dimension(size(Hk,1))                                   :: H0,Hl
    real(8),dimension(size(Hk,1))                                   :: H0free,Hlfree
    real(8),dimension(size(Hk,1))                                   :: H0tmp,Hltmp
    real(8),dimension(size(Hk,1))                                   :: Ekin_,Eloc_

    !
    !MPI setup:
#ifdef _MPI    
    if(check_MPI())then
       mpi_size  = get_size_MPI()
       mpi_rank =  get_rank_MPI()
       mpi_master= get_master_MPI()
    else
       mpi_size=1
       mpi_rank=0
       mpi_master=.true.
    endif
#else
    mpi_size=1
    mpi_rank=0
    mpi_master=.true.
#endif
    !
    !Retrieve parameters:
    ! call get_ctrl_var(Norb,"NORB")
    ! call get_ctrl_var(Nspin,"NSPIN")
    ! call get_ctrl_var(beta,"BETA")
    ! call get_ctrl_var(xmu,"XMU")
    !Get generalized Lattice-Spin-Orbital index
    Nlso = size(Hk,1)
    Lk   = size(Hk,3)
    Nlat = size(Sigma,1)
    Nso  = size(Sigma,2)
    Liw  = size(Sigma,4)
    !Testing:
    if(Nso/=Norb*Nspin)stop "dmft_kinetic_energy_exp: Nso != Norb*Nspin [from Sigma]"
    call assert_shape(Hk,[Nlat*Nso,Nlso,Lk],"dmft_kinetic_energy_exp","Hk") !implcitly test that Nlat*Nso=Nlso
    call assert_shape(Sigma,[Nlat,Nso,Nso,Liw],"dmft_kinetic_energy_exp","Sigma")
    !
    !Allocate and setup the Matsubara freq.
    if(allocated(wm))deallocate(wm);allocate(wm(Liw))
    wm = pi/beta*(2*arange(1,Liw)-1)
    !
    !Get HF part of the self-energy
    Sigma_HF = dreal(Sigma(:,:,:,Liw))![Nlat,Nso,Nso]
    !
    ! Get the local Hamiltonian, i.e. the block diagonal part of the full Hk summed over k
    Hloc_tmp=sum(Hk(:,:,:),dim=3)/dble(Lk)
    Hloc=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Nspin*Norb
                   Hloc(is,js)=Hloc_tmp(is,js) 
                enddo
             enddo
          enddo
       enddo
    enddo
    where(abs(dreal(Hloc))<1.d-6)Hloc=0d0
    ! if(size(Hloc,1)<16)then
    !    if(mpi_master)call print_hloc(Hloc)
    ! else
    !    if(mpi_master)call print_hloc(Hloc,"dmft_ekin_Hloc.dat")
    ! endif
    !
    !Start the timer:
    if(mpi_master)write(*,"(A)") "Kinetic energy computation"
    if(mpi_master)call start_timer
    H0    = 0d0
    Hl    = 0d0
    H0tmp = 0d0
    Hltmp = 0d0
    !Get principal part: Tr[ Hk.(Gk-Tk) ]
    do ik=1,Lk
       Ak    = Hk(:,:,ik) - Hloc(:,:)
       do i=1+mpi_rank,Liw,mpi_size
          GkTmp(:,:) = (xi*wm(i) + xmu)*eye(Nlso) - blocks_to_matrix(Sigma(:,:,:,i),Nlat,Nso)  - Hk(:,:,ik)        
          call inv(GkTmp)
          Gk = GkTmp
          !
          GkTmp=zero
          GkTmp(:,:) = (xi*wm(i) + xmu)*eye(Nlso) - blocks_to_matrix(Sigma_HF,Nlat,Nso)  - Hk(:,:,ik)
          call inv(GkTmp)
          Tk = GkTmp
          !
          Bk = matmul(Ak, Gk - Tk)
          Ck = matmul(Hloc, Gk - Tk)
          do is=1,Nlso
             H0tmp(is) = H0tmp(is) + Bk(is,is)/dble(Lk)
             Hltmp(is) = Hltmp(is) + Ck(is,is)/dble(Lk)
          enddo
       enddo
       if(mpi_master)call eta(ik,Lk)
    enddo
#ifdef _MPI
    if(check_MPI())then
       call AllReduce_MPI(MPI_COMM_WORLD,H0tmp,H0)
       call AllReduce_MPI(MPI_COMM_WORLD,Hltmp,Hl)
    else
       H0=H0tmp
       Hl=Hltmp
    endif
#else
    H0=H0tmp
    Hl=Hltmp
#endif
    if(mpi_master)call stop_timer
    spin_degeneracy=3d0-Nspin !2 if Nspin=1, 1 if Nspin=2
    H0 = H0/beta*2d0*spin_degeneracy
    Hl = Hl/beta*2d0*spin_degeneracy
    !
    !
    !get tail subtracted contribution: Tr[ Hk.Tk ]
    H0free=0d0
    Hlfree=0d0
    do ik=1,Lk
       Evec= Hk(:,:,ik) +  blocks_to_matrix(Sigma_HF,Nlat,Nso)     
       call eigh(Evec,Eval)
       Nk  = fermi(Eval,beta)
       Gk  = matmul(Evec,matmul(diag(Nk),conjg(transpose(Evec))))
       Ak  = matmul(Hk(:,:,ik)-Hloc, Gk)
       do is=1,Nlso
          H0free(is) = H0free(is) + Ak(is,is)/dble(Lk)
       enddo
       Ak  = matmul(Hloc, Gk)
       do is=1,Nlso
          Hlfree(is) = Hlfree(is) + Ak(is,is)/dble(Lk)
       enddo
    enddo
    H0free=spin_degeneracy*H0free
    Hlfree=spin_degeneracy*Hlfree
    !
    Ekin_=H0+H0free
    Eloc_=Hl+Hlfree
    !
    if(mpi_master)call write_kinetic_value(Ekin_,Eloc_,Nlat,Nso)
    !
    deallocate(wm)
  end subroutine dmft_kinetic_energy_exp

  subroutine write_kinetic_value(Ekin,Eloc,Nlat,Nso)
    real(8),dimension(:)               :: Ekin
    real(8),dimension(size(Ekin))      :: Eloc
    real(8),dimension(:,:),allocatable :: Ekin_,Eloc_
    integer,optional                   :: Nlat,Nso
    integer                            :: Nlso
    integer                            :: i,iso,ilat,unit
    if(.not.present(Nlat))then
       !
       Nlso = size(Ekin)
       !
       unit = free_unit()
       open(unit,file="my_kinetic_energy.info")
       write(unit,"(A1,90(A14,1X))")"#",&
            str(1)//"<K>",str(2)//"<Eloc>",&
            (str(2+i)//"<K"//str(i)//">",i=1,Nlso),&
            (str(2+Nlso+i)//"<Eloc"//str(i)//">",i=1,Nlso)
       close(unit)
       !
       unit = free_unit()
       open(unit,file="my_kinetic_energy.dat")
       write(unit,"(90F15.9)")sum(Ekin),sum(Eloc),(Ekin(i),i=1,Nlso),(Eloc(i),i=1,Nlso)
       close(unit)
       !
    else
       !
       if(.not.present(Nso))stop "ERROR write_kinetic_value: Nlat present but Nso not present."
       !
       Nlso = size(Ekin)
       if(Nlso /= Nlat*Nso)stop "Error write_kinetic_value: Nlso != Nlat*Nso" 
       !
       unit = free_unit()
       open(unit,file="my_kinetic_energy.info")
       write(unit,"(A1,90(A14,1X))")"#",&
            str(1)//"<K>",str(2)//"<Eloc>",&
            (str(2+i)//"<K"//str(i)//">",i=1,Nso),&
            (str(2+Nso+i)//"<Eloc"//str(i)//">",i=1,Nso)
       close(unit)
       !
       !
       allocate(Ekin_(Nlat,Nso))
       allocate(Eloc_(Nlat,Nso))
       do ilat=1,Nlat
          do iso=1,Nso
             i = iso + (ilat-1)*Nso
             Ekin_(ilat,iso) = Ekin(i)
             Eloc_(ilat,iso) = Eloc(i)
          enddo
       enddo
       !
       unit = free_unit()
       open(unit,file="my_kinetic_energy.dat")       
       write(unit,"(90F15.9)")sum(Ekin_)/Nlat,sum(Eloc_)/Nlat
       do ilat=1,Nlat
          write(unit,"(100000F15.9)")sum(Ekin_(ilat,:)),sum(Eloc_(ilat,:)),&
               (Ekin_(ilat,i),i=1,Nso),&
               (Eloc_(ilat,i),i=1,Nso)
       enddo
       close(unit)
    endif
  end subroutine write_kinetic_value

  !+-----------------------------------------------------------------------------+!
  !PURPOSE:
  ! Bcast/Reduce a vector of Blocks [Nlat][Nso][Nso] onto a matrix [Nlat*Nso][Nlat*Nso]
  !+-----------------------------------------------------------------------------+!
  function blocks_to_matrix(Vblocks,Nlat,Nso) result(Matrix)
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: i,j,ip
    Matrix=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Matrix(i:j,i:j) =  Vblocks(ip,:,:)
    enddo
  end function blocks_to_matrix

  function matrix_to_blocks(Matrix,Nlat,Nso) result(Vblocks)
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: i,j,ip
    Vblocks=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Vblocks(ip,:,:) = Matrix(i:j,i:j)
    enddo
  end function matrix_to_blocks




end program ed_kanemele
