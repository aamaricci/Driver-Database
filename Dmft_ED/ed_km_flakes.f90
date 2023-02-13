program ed_kanemele_flakes
   USE DMFT_ED    !0.6.0
   USE SCIFOR     !4.9.6
   USE DMFT_TOOLS !2.4.3
   USE HONEYTOOLS !0.2.1
   USE HONEYPLOTS !0.2.1
   USE MPI

   implicit none

   integer                                       :: iloop,Nso,Nlso,Nlat
   logical                                       :: converged
   integer                                       :: ilat
   !
   !Bath:
   integer                                       :: Nb
   real(8),allocatable,dimension(:,:)            :: Bath,Bath_prev
   !
   !The local hybridization function:
   complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
   complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Sreal
   complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Greal
   !
   !Hamiltonian input:
   complex(8),allocatable,dimension(:,:,:)       :: Hij
   complex(8),allocatable,dimension(:,:)         :: kmHij
   complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
   complex(8),allocatable,dimension(:,:)         :: Hloc_lso
   !
   real(8),dimension(2)                          :: e1,e2   !real-space lattice basis
   real(8),dimension(2)                          :: d1,d2,d3
   real(8),dimension(2)                          :: a1,a2,a3
   !
   !Variables for the model:
   integer,parameter                             :: Lk=1 ! just one k-point
   integer                                       :: radius
   real(8)                                       :: t1,t2,nphi,phi,Mh,Bz,wmixing
   character(len=32)                             :: finput
   character(len=32)                             :: HijFILE
   !
   !Flags and options
   character(len=32)                             :: bathspins
   logical                                       :: neelsym,xkick,ykick,getbands
   !
   !Replica Hamiltonian
   real(8),dimension(:,:,:),allocatable          :: lambdasym_vectors ![Nlat,Nbath,Nsym]
   complex(8),dimension(:,:,:,:,:),allocatable   :: Hsym_basis   ![size(Hloc),Nsym]
   real(8),dimension(:,:),allocatable            :: onsite_band  !temporary [Nlat,Nbath]
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


   !Parse additional variables && read Input
   call parse_cmd_variable(finput,"FINPUT",default='inputKANEMELE.conf')
   !
   call parse_input_variable(HijFILE,"HijFILE",finput,default="Hij.in",&
      comment='Hk will be written here')
   call parse_input_variable(radius,"RADIUS",finput,default=1,&
      comment='Integer radius of the flake, in hexagonal units')
   call parse_input_variable(t1,"T1",finput,default=1d0,&
      comment='NN hopping, fixes noninteracting bandwidth')
   call parse_input_variable(t2,"T2",finput,default=0.1d0,&
      comment='Haldane-like NNN hopping-strenght, corresponds to lambda_SO in KM notation')
   call parse_input_variable(nphi,"NPHI",finput,default=0.5d0,&
      comment='Haldane-like flux for the SO term, in units of PI')
   call parse_input_variable(mh,"MH",finput,default=0d0,&
      comment='On-site staggering, aka Semenoff-Mass term')
   call parse_input_variable(Bz,"Bz",Finput,default=0d0,&
      comment='External AFMz Zeeman field: Hij = Hij - Bz * tau_z ⊗ sigma_z')
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
   call add_ctrl_var(nbath,"nbath")
   call add_ctrl_var(ed_hw_bath,"ed_hw_bath")


   !INPUT VALIDATION
   !
   if(.not.(bath_type=="replica".AND.ed_mode=='nonsu2'))&
      stop "Wrong setup from input file: AFMxy requires NONSU2-mode and REPLICA-bath"
   if(BathSpins=="xz".AND.yKICK)&
      stop "Wrong setup from input file: AFM(y) sb-field is not allowed with a xz-only bath"
   if(BathSpins=="x".AND.yKICK)&
      stop "Wrong setup from input file: AFM(y) sb-field is not allowed with a x-only bath"
   if(Norb/=1.OR.Nspin/=2)&
      stop "Wrong setup from input file: Norb=1 AND Nspin=2 is the correct configuration."
   !
   Nso=Nspin*Norb ! >> Nlat undefined until nanoflake is built <<

   !SETUP THE HONEYCOMB LATTICE
   !
   !Lattice basis is:
   !e₁ = a₀ [ sqrt3/2 , 1/2 ] = 3/2a[1, 1/sqrt3]
   !e₂ = a₀ [ sqrt3/2 ,-1/2 ] = 3/2a[1,-1/sqrt3]
   e1 = 3d0/2d0*[1d0, 1d0/sqrt(3d0)]
   e2 = 3d0/2d0*[1d0,-1d0/sqrt(3d0)]
   !
   !Unit-cell displacements: nearest neighbor A-->B, B-->A
   d1= [  1d0/2d0 , sqrt(3d0)/2d0 ]
   d2= [  1d0/2d0 ,-sqrt(3d0)/2d0 ]
   d3= [ -1d0     , 0d0           ]
   !
   !Cell displacements: next nearest-neighbor A-->A, B-->B
   a1 = d1-d3                    !3/2*a[1,1/sqrt3]
   a2 = d2-d3                    !3/2*a[1,-1/sqrt3]
   a3 = d1-d2                    !sqrt3[0,1]

   !BUILD THE REAL SPACE HAMILTONIAN
   call build_Hij(trim(HijFILE)) ! {Hij} -> defines Nlat
   allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
   Hloc = lso2nnn_reshape(kmHij,Nlat,Nspin,Norb)
   ! Note that kmHij is not diagonal in Nlat
   Hloc_lso = nnn2lso_reshape(Hloc,Nlat,Nspin,Norb)
   if(master)call TB_write_Hloc(Hloc_lso,'Hloc.txt')

   !ALLOCATE LOCAL FIELDS
   allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
   allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
   allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
   allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
   allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal=zero

   !SETUP HREPLICA SYMMETRIES:
   call build_replica_band(onsite_band,ed_hw_bath,Nbath)
   select case(trim(BathSpins))

    case default
      stop "BathSpins not in [x; xy; xz; xyz]"

    case("x")  !Only X spin component in the bath
      allocate(lambdasym_vectors(Nlat,Nbath,2))
      allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,2))
      Hsym_basis(:,:,:,:,1)=so2nn_reshape(pauli_sigma_0,Nspin,Norb)
      Hsym_basis(:,:,:,:,2)=so2nn_reshape(pauli_sigma_x,Nspin,Norb)
      lambdasym_vectors(:,:,1)=onsite_band
      lambdasym_vectors(:,:,2)=0d0 !unbroken symmetry by default

    case("xy")  !Only XY spin components in the bath
      allocate(lambdasym_vectors(Nlat,Nbath,3))
      allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
      Hsym_basis(:,:,:,:,1)=so2nn_reshape(pauli_sigma_0,Nspin,Norb)
      Hsym_basis(:,:,:,:,2)=so2nn_reshape(pauli_sigma_x,Nspin,Norb)
      Hsym_basis(:,:,:,:,3)=so2nn_reshape(pauli_sigma_y,Nspin,Norb)
      lambdasym_vectors(:,:,1)=onsite_band
      lambdasym_vectors(:,:,2)=0d0 !unbroken magnetic
      lambdasym_vectors(:,:,3)=0d0 !symmetry by default

    case("xz")  !Only XZ spin components in the bath
      allocate(lambdasym_vectors(Nlat,Nbath,3))
      allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
      Hsym_basis(:,:,:,:,1)=so2nn_reshape(pauli_sigma_0,Nspin,Norb)
      Hsym_basis(:,:,:,:,2)=so2nn_reshape(pauli_sigma_x,Nspin,Norb)
      Hsym_basis(:,:,:,:,3)=so2nn_reshape(pauli_sigma_z,Nspin,Norb)
      lambdasym_vectors(:,:,1)=onsite_band
      lambdasym_vectors(:,:,2)=0d0 !unbroken magnetic
      lambdasym_vectors(:,:,3)=0d0 !symmetry by default

    case("xyz") !Full XYZ spin freedom in the bath
      allocate(lambdasym_vectors(Nlat,Nbath,4))
      allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,4))
      Hsym_basis(:,:,:,:,1)=so2nn_reshape(pauli_sigma_0,Nspin,Norb)
      Hsym_basis(:,:,:,:,2)=so2nn_reshape(pauli_sigma_x,Nspin,Norb)
      Hsym_basis(:,:,:,:,3)=so2nn_reshape(pauli_sigma_y,Nspin,Norb)
      Hsym_basis(:,:,:,:,4)=so2nn_reshape(pauli_sigma_z,Nspin,Norb);
      lambdasym_vectors(:,:,1)=onsite_band
      lambdasym_vectors(:,:,2)=0d0 !unbroken magnetic
      lambdasym_vectors(:,:,3)=0d0 !symmetry on all
      lambdasym_vectors(:,:,4)=0d0 !axes by default

   end select

   !SETUP SYMMETRY BREAKING KICKS (AFM order)
   if(xKICK)then
      lambdasym_vectors(1,:,2) = +sb_field
      lambdasym_vectors(2,:,2) = -sb_field
      !For the log file
      if(master)write(*,*) "*************************************************"
      if(master)write(*,*) "*                                               *"
      if(master)write(*,*) "*  !Applying an AFMx kick to the initial bath!  *"
      if(master)write(*,*) "*                                               *"
      if(master)write(*,*) "*************************************************"
   endif
   if(yKICK)then  !Safe: look at the preliminary checks
      lambdasym_vectors(1,:,3) = +sb_field
      lambdasym_vectors(2,:,3) = -sb_field
      !For the log file
      if(master)write(*,*) "*************************************************"
      if(master)write(*,*) "*                                               *"
      if(master)write(*,*) "*  !Applying an AFMy kick to the initial bath!  *"
      if(master)write(*,*) "*                                               *"
      if(master)write(*,*) "*************************************************"
   endif

   !SETUP H_replica
   call ed_set_Hreplica(Hsym_basis,lambdasym_vectors)
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
      !COMPUTE THE LOCAL GF
      call dmft_gloc_matsubara(Hij,Gmats,Smats)
      call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
      !
      !COMPUTE THE WEISS FIELD (only the Nineq ones)
      call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,cg_scheme)
      call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=4)
      !
      !FIT THE NEW BATH, starting from the old bath + the supplied delta
      !IF(NONSU2): Sz-conservation is broken, so a unique fit for both spins is needed
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
   call dmft_gloc_realaxis(Hij,Greal,Sreal)
   call dmft_print_gf_realaxis(Greal,"Greal",iprint=4)

   !Compute Kinetic Energy
   call dmft_kinetic_energy(Hij,Smats)

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
   !PURPOSE: Get real-space Kane Mele Model Hamiltonian
   !---------------------------------------------------------------------
   subroutine build_Hij(filename)
      character(len=*),optional     :: filename
      !
      if(master)write(LOGfile,*)"Build H(ij) for a Kane-Mele flake"
      if(master)write(LOGfile,*)"# of SO-bands:",Nso
      !
      !Get Hlso (automatic allocation)
      kmHij = Hij_kanemele_flake(radius)
      !
      if(allocated(Hij))deallocate(Hij)
      allocate(Hij(Nlso,Nlso,Lk))
      Hij(:,:,Lk) = kmHij
      !
      call TB_write_Hloc(kmHij,filename)
      !
   end subroutine build_Hij


   !--------------------------------------------------------------------!
   ! > Kane-Mele HAMILTONIAN in real-space, through HoneyTools library
   !--------------------------------------------------------------------!
   function Hij_kanemele_flake(radius) result(Hlso)
      integer,intent(in)               :: radius
      complex(8),allocatable           :: Hup(:,:)  ![Nlat,Nlat]
      complex(8),allocatable           :: Hdw(:,:)  ![Nlat,Nlat]
      complex(8),allocatable           :: Hlso(:,:) ![Nlso,Nlso]
      complex(8),allocatable           :: Hnnnn(:,:,:,:,:,:)
      complex(8),allocatable           :: UPstates(:,:)
      real(8),allocatable              :: UPlevels(:)
      complex(8),allocatable           :: DWstates(:,:)
      real(8),allocatable              :: DWlevels(:)
      character(32)                    :: fig_name
      type(unit_cell)                  :: km_basis
      type(xy_lattice)                 :: km_flake
      type(xy_lattice)                 :: subflake
      type(xy_site)                    :: site
      integer,allocatable              :: indices(:)
      logical,allocatable              :: t1_mask(:,:)
      logical,allocatable              :: t2_mask(:,:)
      integer                          :: unit
      integer                          :: i,j,k,l
      !
      km_basis = unit_cell(hex_orientation(e1,e2,angle=0))
      km_flake = get_flake(radius,layout=km_basis)
      call xy_next_nearest_neighbors(lattice=km_flake,nn_mask=t1_mask,nnn_mask=t2_mask)
      !
      !Determine Nlat and Nlso
      Nlat = km_flake%size
      Nlso = Nlat*Nso
      if(master)write(LOGfile,*)"# of sites:",Nlat
      !
      !Build spin hamiltonians ([Nlat,Nlat])
      allocate(Hup(Nlat,Nlat),Hdw(Nlat,Nlat))
      Hup = zero
      Hdw = zero
      !HOPPING AMPLITUDES: EASY!
      where(t1_mask)
         Hup = -t1
         Hdw = -t1
      end where
      !SUBLATTICE TERMS: EASY!
      subflake = get_sublattice(km_flake,"A")
      indices = subflake%site%key
      do i = 1,size(indices)
         Hup(indices(i),indices(i)) = + Mh - Bz
         Hdw(indices(i),indices(i)) = + Mh + Bz
      enddo
      subflake = get_sublattice(km_flake,"B")
      indices = subflake%site%key
      do i = 1,size(indices)
         Hup(indices(i),indices(i)) = - Mh + Bz
         Hdw(indices(i),indices(i)) = - Mh - Bz
      enddo
      !NOW THE PAINFUL SOC PHASES <'TT_TT'>
      phi = nphi * pi
      l = 0 !counter
      do i = 1,km_flake%size
         do j = i+1,km_flake%size
            do k = 1,6
               site = km_flake%site(i)
               site = xy_nnn_hop(km_basis,site,k)
               if(site==km_flake%site(j))then
                  if(mod(k,2)==1)then
                     if(km_flake%site(i)%label=="A")then
                        Hup(i,j) = +t2 * exp(+xi*phi)
                        Hdw(i,j) = -t2 * exp(+xi*phi)
                     else
                        Hup(i,j) = +t2 * exp(-xi*phi)
                        Hdw(i,j) = -t2 * exp(-xi*phi)
                     endif
                  else
                     if(km_flake%site(i)%label=="A")then
                        Hup(i,j) = +t2 * exp(-xi*phi)
                        Hdw(i,j) = -t2 * exp(-xi*phi)
                     else
                        Hup(i,j) = +t2 * exp(+xi*phi)
                        Hdw(i,j) = -t2 * exp(+xi*phi)
                     endif
                  endif
                  l = l + 2
                  Hup(j,i) = conjg(Hup(i,j))
                  Hdw(j,i) = conjg(Hdw(i,j))
               endif
            enddo
         enddo
      enddo
      if(l/=count(t2_mask)) error stop "WRONG NN COUNT"
      !
      !Print to file Hup and Hdw
      call TB_write_Hloc(Hup,"Hup.txt")
      call TB_write_Hloc(Hdw,"Hdw.txt")
      !
      if(getbands)then
         if(master)write(*,*) "***************************************"
         if(master)write(*,*) "*                                     *"
         if(master)write(*,*) "*  !Solving noninteracting TB model!  *"
         if(master)write(*,*) "*                                     *"
         if(master)write(*,*) "***************************************"
         UPstates = Hup
         DWstates = Hdw
         allocate(UPlevels(Nlat))
         allocate(DWlevels(Nlat))
         call eigh(UPstates,UPlevels,jobz='V',uplo='U')
         call eigh(DWstates,DWlevels,jobz='V',uplo='U')
         call TB_write_Hloc(UPstates,'km_up_states.txt')
         call TB_write_Hloc(DWstates,'km_dw_states.txt')
         if(master)call save_array("km_up_levels.txt",UPlevels)
         if(master)call save_array("km_dw_levels.txt",DWlevels)
      endif
      !
      !FILL IN THE STANDARD QCMPLAB HIGHER-RANK FORMAT
      allocate(Hnnnn(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
      Hnnnn = zero
      do i = 1,Nlat
         Hnnnn(i,i,1,1,1,1) = Hup(i,i)
         Hnnnn(i,i,2,2,1,1) = Hdw(i,i)
         do j = 1,Nlat
            Hnnnn(i,j,1,1,1,1) = Hup(i,j)
            Hnnnn(j,i,1,1,1,1) = Hup(j,i)
            Hnnnn(i,j,2,2,1,1) = Hdw(i,j)
            Hnnnn(j,i,2,2,1,1) = Hdw(j,i)
         enddo
      enddo
      !
      !RESHAPE TO NLSO FORMAT...
      Hlso = nnnn2lso_reshape(Hnnnn,Nlat,Nspin,Norb)
      !
      !Plotting calls and lattice I/O
      call plot(km_flake,backend='gnuplot',set_terminal='dumb')
      fig_name = trim('flake'//str(radius)//'.svg')
      call plot(km_flake,t1_mask,figure_name=fig_name)
      unit = free_unit()
      open(unit,file="flake.txt",action="write",position="rewind")
      call xy_print(km_flake,unit,quiet=.true.)
      close(unit)
      !
   end function Hij_kanemele_flake





   !--------------------------------------------------------------------!
   !Reshaping functions:                                                !
   !--------------------------------------------------------------------!

   function nnnn2lso_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
      !! 🚨 THIS VERSION IS NOT DIAGONAL IN NLAT 🚨 !!
      integer                                               :: Nlat,Nspin,Norb
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Fin
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Fout
      integer                                               :: iorb,ispin,ilat,is
      integer                                               :: jorb,jspin,jlat,js
      Fout=zero
      do ilat=1,Nlat
         do jlat=1,Nlat
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        !lattice-spin-orbit stride
                        is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                        js = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin
                        Fout(is,js) = Fin(ilat,jlat,ispin,jspin,iorb,jorb)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function nnnn2lso_reshape

   function lso2nnnn_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
      !! 🚨 THIS VERSION IS NOT DIAGONAL IN NLAT 🚨 !!
      integer                                               :: Nlat,Nspin,Norb
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Fin
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Fout
      integer                                               :: iorb,ispin,ilat,is
      integer                                               :: jorb,jspin,jlat,js
      Fout=zero
      do ilat=1,Nlat
         do jlat=1,Nlat
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        !lattice-spin-orbit stride
                        is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                        js = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin
                        Fout(ilat,jlat,ispin,jspin,iorb,jorb) = Fin(is,js)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function lso2nnnn_reshape

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
                     !lattice-spin-orbit stride
                     is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                     js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
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
                     !lattice-spin-orbit stride
                     is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
                     js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin
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
      integer                                               :: iorb,ispin,is
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
      integer                                               :: iorb,ispin,is
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


   subroutine build_replica_band(lambdas,bandwidth,Nreplica)
      real(8),allocatable,dimension(:,:),intent(out)   :: lambdas
      real(8),intent(in)                               :: bandwidth
      integer,intent(in)                               :: Nreplica
      integer                                          :: ireplica
      real(8),dimension(Nreplica)                      :: tempvec
      real(8)                                          :: tempval
      !
      allocate(lambdas(Nlat,Nbath))
      !
      do ireplica=1,Nreplica
         tempval = ireplica - 1 - (Nreplica-1)/2d0       ![-(N-1)/2 : (N-1)/2]
         tempval = tempval * 2 * bandwidth/(Nreplica-1)  !-bandwidth:bandwidth
         tempvec(ireplica) = tempval
      enddo
      !
      if(mod(Nreplica,2)==0)then
         tempvec(Nreplica/2) = -1d-1  !Much needed small energies around the
         tempvec(Nreplica/2+1) = 1d-1 !Fermi level (if EF itself not present)
      endif
      !
      do ilat=1,Nlat
         lambdas(ilat,:) = tempvec
      enddo
      !
   end subroutine build_replica_band

end program ed_kanemele_flakes
