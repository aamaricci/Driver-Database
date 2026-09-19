program BHZ_1d
  USE SCIFOR
  USE DMRG
#ifdef _MPI
  USE MPI
#endif
  implicit none

  integer                                        :: Nso
  character(len=64)                              :: finput
  integer                                        :: i,j,unit,iorb,ispin,Nsites,N3d
  real(8)                                        :: eh,mh,lambda,K
  real(8)                                        :: Eloc,Etot,Product
  type(site)                                     :: Dot
  complex(8),dimension(4,4)                      :: GammaZ,GammaX
  complex(8),dimension(:,:),allocatable          :: Hloc,Hlr
  type(sparse_matrix),dimension(:,:),allocatable :: N,C
  type(sparse_matrix),dimension(:),allocatable   :: dens,docc,sz,s2z,Mvec
  type(sparse_matrix)                            :: Tz,Kij,Hi,Tx,O4f
  complex(8)                                     :: Corr
  integer                                        :: irank,comm,rank,ierr
  logical                                        :: master,imeasure,irun,get3Dcorr,getAVcorr
  real(8),dimension(:,:),allocatable             :: dqs

#ifdef _MPI  
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
#endif


  call parse_cmd_variable(finput,"FINPUT",default='DMRG.conf')
  call parse_input_variable(irun,"irun",finput,default=.true.,comment="Bool to run DMRG. F for post-processing")
  call parse_input_variable(imeasure,"imeasure",finput,default=.true.,comment="Bool to perform measurements. T for post-processing.")
  call parse_input_variable(mh,"MH",finput,default=0.5d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_input_variable(getAVcorr,"GETAVCORR",finput,default=.false.)
  call parse_input_variable(get3dcorr,"GET3DCORR",finput,default=.false.)
  call read_input(finput)


  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"

  if(Imeasure)then
    save_block=.true.
    save_umat=.true.
  endif       

  Nsites=2*Ldmrg
  Nso=Nspin*Norb
  !
  gammaX=kron( pauli_sigma_z, pauli_tau_x)
  gammaZ=kron( pauli_sigma_0, pauli_tau_z)

  !>Local Hamiltonian:
  allocate(Hloc(Nso,Nso))
  Hloc = Mh*GammaZ
  Dot  = electron_site(Hloc)

  !>Hopping Hamiltonian (i->i+1, right hop direction)
  if(allocated(Hlr))deallocate(Hlr)
  allocate(Hlr(Nso,Nso))
  Hlr = -0.5d0*GammaZ + 0.5d0*xi*lambda*GammaX
  if(master)then
    call print_matrix(Hloc,"Hloc.dmrg")
    call print_matrix(Hlr,"Hlr.dmrg")
  endif  

  !Init DMRG
  call init_dmrg(Hlr,ModelDot=[Dot])

  !Run DMRG algorithm
  if(Irun)call run_DMRG()

  if(Imeasure)then
     !Post-processing and measure quantities:
     allocate(C(Norb,Nspin),N(Norb,Nspin))
     do ispin=1,Nspin
        do iorb=1,Norb
           C(iorb,ispin) = dot%operators%op(key="C"//dot%okey(iorb,ispin))
           N(iorb,ispin) = matmul(C(iorb,ispin)%dgr(),C(iorb,ispin))
        enddo
     enddo
     allocate(Mvec(3*Norb),sz(Norb))
     do iorb=1,Norb
        sz(iorb)          = n(iorb,1)-n(iorb,2)
        Mvec(iorb)        = n(iorb,1)+n(iorb,2)
        Mvec(iorb+Norb)   = matmul(n(iorb,1),n(iorb,2))
        Mvec(iorb+2*Norb) = matmul(sz(iorb),sz(iorb))
     enddo
     Tz = (n(2,1)+n(2,2)-n(1,1)-n(1,2))/2d0
     !
     !Measure <Tz>,<N>,<D>,<Sz^2> (write to file integrated)
     call Measure_DMRG([Tz,Mvec],file="tz_n_d_s2zVSj", pos=arange(1,Nsites))
     !

     !Measure <Tz(1).Tz(j)>-<Tz(1)><Tz(j)>
     if(master)print*,"<Tz_i.Tz_j>-<Tz_i><Tz_j>"     
     call get_correlations("tz.tz",Tz,[0d0,0d0],Tz,[0d0,0d0])
     !Measure <Sz(i).Sz(j)>
     if(master)print*,"<Sz_i.Sz_j>"
     call get_correlations("sz.sz",Sz(1)+Sz(2),[0d0,0d0],Sz(1)+Sz(2),[0d0,0d0])
     !Measure <Tx(i).Tx(j)>
     !Tx = 1/2\sum_\sigma c^+_{1i\sigma}.c_{2i\sigma} + c^+_{2i\sigma}.c_{1i\sigma}
     if(master)print*,"<Tx_i.Tx_j>"
     Tx = 0.5d0*(matmul(C(1,1)%dgr(),C(2,1))+matmul(C(2,1)%dgr(),C(1,1)))
     Tx = Tx + 0.5d0*(matmul(C(1,2)%dgr(),C(2,2))+matmul(C(2,2)%dgr(),C(1,2)))
     call get_correlations("tx.tx",Tx,[0d0,0d0],Tx,[0d0,0d0])
     !
     if(master)print*,"<Tz_i.Tz_{i+1}.Tz_j.Tz_{j+1}>"
     allocate(dqs(2,4));dqs=0d0
     if(master)unit=fopen("O4fz.O4fz_ij"//str(label_DMRG('u')),append=.true.)
     call start_timer()
     do i=1,Nsites-1
      do j=i,Nsites-1
        product = dreal(Measure_Product_DMRG([Tz,Tz,Tz,Tz],dqs,["none","none","none","none"],[i,i+1,j,j+1]))
        if(master)write(unit,*)i,j,product
      enddo
      if(master)write(unit,*)""
      if(master)call eta(i,Nsites-1)
     enddo
     call stop_timer()
     
     !Measure energies: <K>,<Hloc>
     if(master)print*,"measure energies"
     if(master)unit=fopen("Ekin_Eloc_Etot"//str(label_DMRG('u')),append=.true.)
     call Measure_Energy_DMRG(Hlr,K,Eloc,Etot,Kij)
     if(master)write(unit,*)K,Eloc,Etot
     if(master)close(unit)


     call End_Measure_DMRG()
     do ispin=1,Nspin
      do iorb=1,Norb
        call C(iorb,ispin)%free()
        call N(iorb,ispin)%free()
      enddo
    enddo
  endif


  !Finalize DMRG
  call finalize_dmrg()
#ifdef _MPI
  call finalize_MPI()
#endif



contains


  subroutine get_correlations(label,OpA,dqA,OpB,dqB)
    character(len=*)                   :: label
    type(sparse_matrix)                :: OpA,OpB
    real(8),dimension(:)               :: dqA,dqB
    integer                            :: i,j,r,ic,count,cut
    real(8),dimension(:,:),allocatable :: Cij
    real(8),dimension(:),allocatable   :: Cavg
    real(8)                            :: sum_corr
    !
    allocate(Cij(Nsites,Nsites))
    Cij=zero
    !Catch them all..
    
    !Plot <O_1.O_j>
    if(master)call start_timer("<O1.Oj>")
    do j=1,Nsites
      Cij(1,j) = dreal(Measure_Corr_DMRG(OpA,dqA,OpB,dqB,1,j,connected=.true.))
      if(master)call eta(j,Nsites)
    enddo
    if(master)call stop_timer()
    if(master) call splot(str(label)//"_1j"//str(label_DMRG('u')), 1d0*arange(1,Nsites), Cij(1,:))
    !
    !
    !Plot <O_i.O_j> with |j-i| as Manhattan distance
    if(master) unit = fopen(str(label)//"_r"//str(label_DMRG('u')), append=.false.)
    if(master)call start_timer("<Oi.Oj>_|i-j|")
    ic = Nsites / 2      
    do r = 1, ic - 1
      i = ic - r / 2
      j = ic + (r + 1) / 2
      if(Cij(i,j)==zero)Cij(i,j) = dreal(Measure_Corr_DMRG(OpA,dqA,OpB,dqB,i,j,connected=.true.))
      if(master) write(unit,*) abs(i-j), Cij(i,j)
      if(master)call eta(r,ic)
    enddo
    if(master)call stop_timer()
    if(master) close(unit)
    !
    !
    !Plot <av{O_i.O_j}> average over radius R as a function of R
    if(getAVcorr)then
      if(master) unit = fopen(str(label)//"_Rav"//str(label_DMRG('u')), append=.false.)   
      if(master)call start_timer("<av{Oi.Oj}>_R")  
      cut = 2   !max(2, Nsites / 2) 
      allocate(Cavg(Nsites - 2*cut))
      Cavg = zero
      do r = 1, Nsites - 2*cut
        sum_corr = zero
        count    = 0
        do i = 1+cut, Nsites-cut-r
          j = i + r
          if(Cij(i,j)==zero)Cij(i,j) = dreal(Measure_Corr_DMRG(OpA,dqA,OpB,dqB,i,j,connected=.true.))
          sum_corr = sum_corr + Cij(i,j)
          count    = count + 1
        enddo
        if (count > 0) Cavg(r) = sum_corr /dble(count)
        if(master) write(unit,*) r, Cavg(r)
        if(master)call eta(r,Nsites-2*cut)
      enddo
      if(master)call stop_timer()
      if(master) close(unit)
      deallocate(Cavg)
    endif
    !
    !
    if(get3dcorr)then
      if(master)call start_timer()
      do i=1,Nsites
        if(Cij(i,j)==zero)Cij(i,i) = dreal(Measure_Corr_DMRG(OpA,dqA,OpB,dqB,i,i,connected=.true.))
        do j=i+1,Nsites
          if(Cij(i,j)==zero)Cij(i,j) = dreal(Measure_Corr_DMRG(OpA,dqA,OpB,dqB,i,j,connected=.true.))
          Cij(j,i) = Cij(i,j)
        enddo
        if(master) call eta(i,Nsites)         
      enddo
      if(master)call stop_timer()
      if(master) call splot3d(str(label)//"_ij"//str(label_DMRG('u')), 1d0*arange(1,Nsites), 1d0*arange(1,Nsites), Cij)
    endif
    !
  end subroutine get_correlations

end program BHZ_1d



