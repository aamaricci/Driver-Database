program ed_bilayer
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer :: ip,im,ipp,iw
  real(8) :: wx
  real(8) :: Vex
  real(8) :: Ef
  real(8) :: h0
  !
  real(8) :: ecut,pcut
  integer :: lp,lphi  
  integer :: unit_io
  !
  real(8),dimension(:),allocatable :: wr,modp,phip  
  complex(8),dimension(:,:),allocatable :: Pi_irreducible,Pi_irreducible_tmp
  !
  real(8),dimension(:),allocatable :: Im_Sigma,Im_Sigma_tmp
  complex(8),dimension(:),allocatable :: Sigma,DX

  !
  complex(8),dimension(:,:),allocatable :: Lambda_vertex
  !
  complex(8),dimension(:),allocatable :: int_tmp,int_tmp_inn

  real(8),dimension(:),allocatable :: int_tmp_im

  real(8),dimension(:),allocatable :: epp_tmp
  real(8)   :: epsilonp
  complex(8) :: wcmplx
  real(8) :: wp,mx,mel,bwp,wtmp,ImL
  !  
  character(len=16)                           :: finput
  !MPI Vars:
  integer                                     :: comm,rank,mpierr,mpiSize
  logical                                     :: master
  !
  !
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  mpiSize = get_Size_MPI(comm)
  !
  !

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='input_Xpol.in')  
  call parse_input_variable(wx,"wx",finput,default=1.d0,comment='q=0  exciton frequency ')
  call parse_input_variable(Vex,"Vex",finput,default=1.d0,comment='e-X interaction')
  call parse_input_variable(ef,"ef",finput,default=0.d0,comment='fermi energy of the 2d electron gas')
  
  call parse_input_variable(ecut,"ecut",finput,default=1.d0,comment='cut off energy for momentum integration [eV]')  
  call parse_input_variable(lp,"lp",finput,default=100,comment='number of points for momenutm integration (modulus)')
  call parse_input_variable(lphi,"lphi",finput,default=100,comment='number of points for momenutm integration (phase)')


  call parse_input_variable(mX,"mX",finput,default=1.d0,comment='mass enhancement exciton')
  call parse_input_variable(mel,"mel",finput,default=1.d0,comment='mass enhancement electron')

  !
  call ed_read_input(trim(finput),comm)
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  !
  !
  h0=7.62d-02 ! hbar^2/m [eV x nm^2]
  pcut=sqrt(2.d0*ecut/h0)
  !
  ! 
  allocate(wr(lreal));  wr=linspace(wini,wfin,lreal)
  allocate(modp(lp));   modp=linspace(0.d0,pcut,lp)
  allocate(phip(lphi)); phip=linspace(0.d0,2.d0*pi,lphi)
  !
  allocate(Pi_irreducible(lreal,lp));     Pi_irreducible = 0.d0
  allocate(Pi_irreducible_tmp(lreal,lp)); Pi_irreducible_tmp = 0.d0
  !  
  allocate(int_tmp(lp)); int_tmp=0.d0
  allocate(int_tmp_inn(lphi)); int_tmp_inn=0.d0

  allocate(epp_tmp(lphi)); epp_tmp=0.d0

  do iw=1+rank,lreal,mpiSize     
     wcmplx=wr(iw) + xi*eps
     
     if(master) write(*,*) iw
     
     do ip=1,lp
        !
        Pi_irreducible_tmp(iw,ip) = 0.d0
        !
        int_tmp=0d0
        do ipp=1,lp
           !
           int_tmp_inn = 0.d0
           epp_tmp(1:lphi) = h0*0.5*(modp(ip)**2.d0+modp(ipp)**2.d0-2.d0*modp(ipp)*modp(ip)*dcos(phip(1:lphi)))/mel-Ef
           wp = wx + h0*0.5*0.5*modp(ipp)**2.d0/mx        
           bwp = 1.d0/(exp(-beta*wp)-1.d0)
           int_tmp_inn(1:lphi) = (fermi(epp_tmp,beta) + bwp)/(wcmplx - epp_tmp - wp)*modp(ipp)/(2d0*pi)/(2d0*pi)
           !
           int_tmp(ipp) = trapz(int_tmp_inn,phip(1),phip(lphi))           
           !
        end do
        Pi_irreducible_tmp(iw,ip) = 1.d0*trapz(int_tmp,modp(1),modp(lp))
        !
     end do
  end do
  !
  call mpi_allreduce(Pi_irreducible_tmp,Pi_irreducible,lreal*lp,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  
  allocate(Lambda_vertex(lreal,lp));     Lambda_vertex = 0.d0

  Lambda_vertex=Vex*Vex*Pi_irreducible/(1.0+Vex*Pi_irreducible)

  if(master) then     
     unit_io=free_unit()
     open(unit=unit_io,file="pi_irreducible.out")
     do iw=1,lreal        
        do ip=1,lp
           write(unit_io,'(5F18.10)') wr(iw),modp(ip),Pi_irreducible(iw,ip)
        end do
        write(unit_io,*)
     end do
     close(unit_io)
     !
     open(unit=unit_io,file="Lambda_vertex.out")
     do iw=1,lreal        
        do ip=1,lp
           write(unit_io,'(5F18.10)') wr(iw),modp(ip),lambda_vertex(iw,ip)
        end do
        write(unit_io,*)
     end do
     close(unit_io)

     open(unit=unit_io,file="q0_pi_irreducible.out")
     do iw=1,lreal        
        write(unit_io,'(5F18.10)') wr(iw),Pi_irreducible(iw,1)
     end do     
     close(unit_io)

     open(unit=unit_io,file="q0_lambda_vertex.out")
     do iw=1,lreal        
        write(unit_io,'(5F18.10)') wr(iw),Lambda_vertex(iw,1)
     end do     
     close(unit_io)

  end if
  !
  allocate(Im_Sigma(lreal));     Im_Sigma = 0.d0
  allocate(Im_Sigma_tmp(lreal)); Im_Sigma_tmp = 0.d0


  allocate(Sigma(lreal));     Sigma = 0.d0

  allocate(int_tmp_im(lp));int_tmp_im = 0.d0
  !
  do iw=1+rank,lreal,mpiSize
     !
     Im_Sigma_tmp(iw) = 0.d0
     !
     int_tmp_im=0.d0
     do ip=1,lp

        !find Im \Lambda_p(\Omega+\epsilon_p)
        epsilonp = h0*0.5*modp(ip)**2.d0
        wtmp = wr(iw) + epsilonp 
        if(wtmp.lt.wini.or.wtmp.gt.wfin) then
           ImL=0.d0
        else
           !+- interpolate
           call linear_spline(wr(:),dimag(Lambda_vertex(:,ip)),wtmp,ImL)
        end if
        !
        int_tmp_im(ip) = (fermi(epsilonp-ef,beta) - fermi(wtmp-ef,beta) )*ImL*modp(ip)/(2.d0*pi)
        !
     end do
     !
     Im_sigma_tmp(iw) = -1.d0*trapz(int_tmp_im,modp(1),modp(lp))           
     !
  end do
  !
  call mpi_allreduce(Im_sigma_tmp,Im_sigma,lreal,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  call get_kkt(Im_sigma,Im_sigma_tmp,wr,'IR')
  !
  Sigma = Im_sigma_tmp + xi*Im_sigma
  !
  
  allocate(Dx(lreal));Dx=0.d0
  do iw=1,Lreal
     wcmplx = wr(iw) + xi*eps
     Dx(iw) = 1.d0/(wcmplx-wx-Sigma(iw))
  end do
  if(rank==0) then
     unit_io=free_unit()
     open(unit=unit_io,file="Sigma.out")
     do iw=1,lreal        
        write(unit_io,'(5F18.10)') wr(iw),Sigma(iw)
     end do
     close(unit_io)

     open(unit=unit_io,file="DX.out")
     do iw=1,lreal        
        write(unit_io,'(5F18.10)') wr(iw),Dx(iw)
     end do
     close(unit_io)

  end if

  !
  call finalize_MPI()

contains
  


  subroutine get_KKT(ReS,ImS,wreal,mode_)
    real(8),dimension(:) :: ReS,ImS
    real(8),dimension(:) :: wreal
    character(len=2),optional :: mode_
    character(len=2) :: mode
    real(8) :: A,B
    integer :: iv,iw,Lw
    real(8),dimension(:),allocatable :: ImS_tmp
    !
    Lw = size(wreal)
    if(size(ReS).ne.Lw) then
       if(rank==0) write(*,*) 'size(ReS).ne.Lw'
       CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)       
       stop
    end if
    !
    if(size(ImS).ne.Lw) then
       if(rank==0) write(*,*) 'size(ImS).ne.Lw'
       CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)       
       stop
    end if
    !
    mode='RI'
    if(present(mode_)) mode=mode_
    if(mode.ne.'IR'.and.mode.ne.'RI') stop "wrong mode KKT"
    !
    allocate(ImS_tmp(Lw))
    ImS_tmp=0.d0
    ImS=0.d0
    do iw=1+rank,Lw,mpiSize
       do iv=1,Lw-1          
          A = ReS(iv) -wr(iv)*(ReS(iv)-ReS(iv+1))/(wr(iv)-wr(iv+1))
          B = (ReS(iv)-ReS(iv+1))/(wr(iv)-wr(iv+1))          
          ImS_tmp(iw) =ImS_tmp(iw)  -B*(wr(iv+1)-wr(iv))          
          if(iv+1.ne.iw) ImS_tmp(iw) = ImS_tmp(iw) - (A+B*wr(iw))*log(abs(wr(iw)-wr(iv+1)))
          if(iv.ne.iw)   ImS_tmp(iw) = ImS_tmp(iw) + (A+B*wr(iw))*log(abs(wr(iw)-wr(iv)))
       end do
       ImS_tmp(iw) = ImS_tmp(iw)/pi
       if(mode.eq.'IR') ImS_tmp(iw)=-ImS_tmp(iw)
    end do
    CALL MPI_ALLREDUCE(ImS_tmp,ImS,Lw,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    !
  end subroutine get_KKT



end program
