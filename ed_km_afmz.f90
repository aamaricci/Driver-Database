program ed_kanemele
   USE DMFT_ED    !0.6.0
   USE SCIFOR     !4.9.4
   USE DMFT_TOOLS !2.3.8
   USE MPI
   implicit none

   integer                                       :: iloop,Lk,Nso,Nlso,Nlat,Nineq
   logical                                       :: converged
   integer                                       :: ispin,ilat!,i,j
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
   real(8),allocatable,dimension(:)              :: dens
   !
   !Variables for the model:
   integer                                       :: Nk,Nkpath
   real(8)                                       :: t1,t2,phi,Mh,Bz,wmixing
   character(len=32)                             :: finput
   character(len=32)                             :: hkfile
   !
   !Flags and options
   logical                                       :: neelsym,afmkick,getbands
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
   call parse_input_variable(Bz,"Bz",Finput,default=0d0,&
      comment='External AFMz Zeeman field: Hk = Hk - Bz * tau_z ⊗ sigma_z')
   call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0,&
      comment='Mixing parameter: 0 means 100% of the old bath (no update at all), 1 means 100% of the new bath (pure update)')
   call parse_input_variable(neelsym,"NEELSYM",finput,default=.true.,&
      comment='If T AFM(z) symmetry is enforced on the self energies at each loop')
   call parse_input_variable(afmkick,"AFMKICK",finput,default=.true.,&
      comment='If T the bath spins get an initial AFM(z) distortion')
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
   if(ed_mode=="superc")stop "Wrong setup from input file: NORMAL or NONSU2 ed-mode here"
   !
   if(bath_type=="hybrid")stop "Wrong setup from input file: NORMAL or REPLICA bath here"
   !
   if(Norb/=1.OR.Nspin/=2)stop "Wrong setup from input file: Norb=1 AND Nspin=2 is the correct configuration for the model"
   Nlat=2
   Nso=Nspin*Norb
   Nlso=Nlat*Nso                 !4 = 2(ineq sites)*2(spin)*1(orb)

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

   !RECIPROCAL LATTICE VECTORS
   bklen=2d0*pi/3d0
   bk1=bklen*[ 1d0, sqrt(3d0)]
   bk2=bklen*[ 1d0,-sqrt(3d0)]
   call TB_set_bk(bkx=bk1,bky=bk2)


   !BUILD THE RECIPROCAL SPACE HAMILTONIAN
   call build_hk(trim(hkfile),getbands)
   allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb));Hloc=zero
   Hloc = lso2nnn_reshape(kmHloc,Nlat,Nspin,Norb)

   !ALLOCATE LOCAL FIELDS
   allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
   allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
   allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
   allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
   allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal=zero

   !SETUP SOLVER
   if(bath_type=="replica")then
      !SETUP REPLICA BATH
      !setup symmetry-basis:
      allocate(lambdasym_vectors(Nlat,Nbath,2))
      allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,2))
      Hsym_basis(:,:,:,:,1)=so2nn_reshape(pauli_sigma_0,Nspin,Norb)
      Hsym_basis(:,:,:,:,2)=so2nn_reshape(pauli_sigma_z,Nspin,Norb)
      call build_replica_band(onsite_band,ed_hw_bath,Nbath)
      lambdasym_vectors(:,:,1)=onsite_band
      lambdasym_vectors(:,:,2)=0d0 !unbroken
      !setup symmetry-breaking AFM kick, if requested in the inputfile
      if(afmkick)then
         lambdasym_vectors(1,:,2) = +sb_field
         lambdasym_vectors(2,:,2) = -sb_field
      endif
      !setup H_replica
      call ed_set_Hreplica(Hsym_basis,lambdasym_vectors)
      !this is now elevated to RDMFT: ineq sites (1,2) for the lambdas
      Nb=ed_get_bath_dimension(Hsym_basis)
      allocate(Bath(Nlat,Nb))
      allocate(Bath_prev(Nlat,Nb))
      call ed_init_solver(comm,Bath)
   elseif(bath_type=="normal")then
      !SETUP NORMAL BATH
      !Usual bath initialization (no H_replica)
      Nb=ed_get_bath_dimension()
      allocate(Bath(Nlat,Nb))
      allocate(Bath_prev(Nlat,Nb))
      call ed_init_solver(comm,Bath)
      !setup symmetry-breaking AFM kick, if requested in the inputfile
      if(afmkick)then
         do ilat=1,Nlat
            call ed_break_symmetry_bath(Bath(ilat,:),sb_field,(-1d0)**(ilat+1))
         enddo
      endif
   endif
   if(afmkick)then !For the log file
      if(master)write(*,*) "*************************************************"
      if(master)write(*,*) "*                                               *"
      if(master)write(*,*) "*  !Applying an AFMz kick to the initial bath!  *"
      if(master)write(*,*) "*                                               *"
      if(master)write(*,*) "*************************************************"
   endif

   !DMFT loop
   iloop=0;converged=.false.
   do while(.not.converged.AND.iloop<nloop)
      iloop=iloop+1
      call start_loop(iloop,nloop,"DMFT-loop")
      !
      !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
      if(neelsym)then
         !solve just one sublattice and get the other by Neel symmetry (z version)
         call ed_solve(comm,Bath(1,:),Hloc(1,:,:,:,:))
         call ed_get_sigma_matsubara(Smats(1,:,:,:,:,:))
         call ed_get_sigma_realaxis(Sreal(1,:,:,:,:,:))
         Smats(2,2,2,:,:,:) = Smats(1,1,1,:,:,:) !S(iw)_{B,dw,dw} = S(iw)_{A,up,up}
         Smats(2,1,1,:,:,:) = Smats(1,2,2,:,:,:) !S(iw)_{B,up,up} = S(iw)_{A,dw,dw}
         Sreal(2,2,2,:,:,:) = Sreal(1,1,1,:,:,:) !S(w)_{B,dw,dw}  = S(w)_{A,up,up}
         Sreal(2,1,1,:,:,:) = Sreal(1,2,2,:,:,:) !S(w)_{B,up,up}  = S(w)_{A,dw,dw}
         if(master)write(*,*) "************************************"
         if(master)write(*,*) "*                                  *"
         if(master)write(*,*) "*   !Enforcing NEEL(z) symmetry!   *"
         if(master)write(*,*) "*                                  *"
         if(master)write(*,*) "************************************"
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
      call dmft_gloc_matsubara(Hk,Gmats,Smats)
      call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4)
      !
      !COMPUTE THE WEISS FIELD (only the Nineq ones)
      call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,cg_scheme)
      call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=4)
      !
      !Fit the new bath, starting from the old bath + the supplied delta
      !Behaves differently depending on the ed_mode input:
      !IF(NORMAL): normal/magZ phase is solved, so we can fit spin1 and spin2 INDIPENDENTLY
      !IF(NONSU2): Sz-conservation is broken, so a unique fit for both spins is needed
      select case(ed_mode)
       case default
         stop "ed_mode has to be normal or nonsu2"
       case("normal")
         call ed_chi2_fitgf(comm,Bath,Weiss,Hloc,ispin=1)
         call ed_chi2_fitgf(comm,Bath,Weiss,Hloc,ispin=2)
       case("nonsu2")
         call ed_chi2_fitgf(comm,Bath,Weiss,Hloc)
      end select
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
               colors_name=[red,blue,tomato,aquamarine],& !\psi=[A_up,A_dw,B_up,B_dw]
               points_name=[character(len=10) :: "{/Symbol G}","K","K`","{/Symbol G}"],&
               file="EigenbandsKMH.nint",iproject=.false.)
         endif
      endif
      !
   end subroutine build_hk



   !--------------------------------------------------------------------!
   !Kane-Mele HAMILTONIAN
   !--------------------------------------------------------------------!
   function hk_kanemele_model(kpoint,Nlso) result(hk)
      real(8),dimension(:)            :: kpoint
      integer                         :: Nlso
      complex(8),dimension(2,2)       :: HkUP,HkDW
      complex(8),dimension(Nlso,Nlso) :: Hk
      real(8)                         :: h0,hx,hy,hz
      real(8)                         :: kdote1, kdote2
      !
      kdote1 = dot_product(kpoint,e1)
      kdote2 = dot_product(kpoint,e2)
      !
      h0 = 2*t2*cos(phi)*( cos(kdote1) + cos(kdote2) + cos(kdote1-kdote2) )
      hx = t1*( cos(kdote1) + cos(kdote2) + 1)
      hy = t1*( sin(kdote1) + sin(kdote2) )
      hz = 2*t2*sin(phi)*( sin(kdote1) - sin(kdote2) - sin(kdote1-kdote2) )
      !
      HkUP = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z + Mh*pauli_z
      HkDW = h0*pauli_0 + hx*pauli_x + hy*pauli_y - hz*pauli_z + Mh*pauli_z
      !
      Hk = zero
      !
      Hk(1:3:2,1:3:2) = HkUP !Recall that the spinors are written in the
      Hk(2:4:2,2:4:2) = HkDW !\psi=[A_up,A_dw,B_up,B_dw] convention.
      !
      !EXTERNAL ANTIFERROMAGNETIC ZEEMAN FIELD (z-axis)
      Hk = Hk - Bz * kron_pauli(pauli_z,pauli_z)
      !
   end function hk_kanemele_model





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




end program ed_kanemele
