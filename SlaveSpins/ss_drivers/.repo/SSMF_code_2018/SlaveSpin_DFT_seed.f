CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This program uses Slave-Spin mean-field for calculating the renormalization 
C     factors (mass enhancements) and band shifts for a multi orbital Hubbard model 
C     hamiltonian with semi-circular 
C     circular Density of States (DOS) or for a user-provided tight-binding hamiltonian
C     INPUT:
C
C     OUTPUT:
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      program Slave_spin_DFT

      implicit double precision (a-h, o-z)
      real*8 J,mixlam,mixh,mixZ
      logical Hund_Ising,zeroTemp, usedos, no_ins_seed, DFTplusSS,fixpop
      logical calclambda0,moveon,calckHam,verb
      include 'param.dat'
      real*8 nm(M), x(M)
      real*8 epsm(M),Z(M,M),Zold(M,M),tm(M),hh(M)
      real*8 ffh(M,M),nf(M,M)
      real*8 cm(M),Sx(M),Sz(M),SxSx(M,M),S(M),SdS(M,M),SzSz(M,M,2)
      real*8 hm(M,M),hmold(M,M),lam(M),lamold(M),lamold2(M),lam0(M)
      real*8 testlam(M),testl(M),testZ(M),testdens
      logical check
      character*80 xxx,prefix
      common /init/ tm,epsm,U,J,Beta,xmu,Vhyb,Vnl
      common /newton/ nm,Z
      common / whichorb /Lorb

      open (1,file='input.dat',form='formatted',status='old')
      open (2,file='options.dat', form='formatted', status='old')

      ! reads in the file "options"
      read(2,*)xxx 
      read(2,*)calclambda0, moveon      
      read(2,*)xxx
      read(2,*)prefix, calckHam
      read(2,*)xxx
      read(2,*)verb
      
      if (DFTplusSS) then       ! it uses the user-provided realistic hamiltonian in file [prefix].rHam 
         if (calckHam) then
            call calc_kHAM(trim(prefix)) ! uses [prefix].rHam to generate the k-space Hamiltonian
         else
            call read_kHam(trim(prefix)) ! uses instead the previously calculated [prefix].kHam
         endif
      endif                     

c     Note that the epsm(L) are used also in the ab-initio case (are overwritten in kHAM) 
c     Ab-initio hamiltonian with local hybridization is not implemented
c     Analytical hamiltonian with basis is not implemented

C*********Guess*********************************************************
      open (11,file='seed.dat',form='formatted',status='unknown')
      
      read(11,*,err=10,end=10)xxx
      do L1=1,M
         read (11,*,err=10,end=10)(Z(L1,L2),L2=L1,M)
         if(no_ins_seed) then
            do L2=L1,M
               if (Z(L1,L2).lt.Zlimit) then
                  print*,'reset Z',L1,L2, 'to',Zreset
                  Z(L1,L2)=Zreset
               endif
            enddo
         endif
         do L2=L1+1,M
            Z(L2,L1)=Z(L1,L2)
         enddo
      enddo 

      read(11,*)xxx
      do L1=1,M
         read(11,*,err=10, end=10)lam(L1),lam0(L1)
      enddo

      read(11,*)xxx ! blank lines are jumped
      read(11,*)xmu, dndmu
      xmuold=xmu

      goto 11

 10   print*,"seed absent or broken"
      do L1=1,M    ! this kicks in if guess is absent or broken
         xmu=0.d0
         dndmu=1.0
         lam(L1)=-xmu
         lam0(L1)=0.d0
         do L2=1,M
            Z(L1,L2)=1.d0
         enddo
      enddo
      
 11   continue

C********* calculation of the lambda0's *******************************
C  calculates the proper shifts, at the given filling, 
C  if moveon=.false. it simply writes the right lambda0's and chemical potential xmu (for U=0) for the requested filling
C  if moveon=.true. it discards the xmu(U=0) and moves on to the Slave spin calculation for the given U,J with the calculated lambda0's
      if (calclambda0) then
         print*,'---------------------------------------',
     &'----------------------------'
         print*," calculating lambda0's"
         print*,'---------------------------------------',
     &'----------------------------'

         read(1,*)xxx           ! dummy reading here, except for densfin
         read(1,*)U, ratioJU 
         read(1,*)xxx           !densfin is used in calculation of lambda0s
         read(1,*)Beta,densfin
         if (.not.DFTplusSS) then
            read(1,*)xxx
            read(1,*)(tm(L),L=1,M)
            read(1,*)(epsm(L),L=1,M)
            read(1,*)Vhyb,Vnl
         endif
         
         call calclambdazero(densfin,lam,lam0) 
         
         if(.not.moveon) then
            call writeseed(Z,lam,lam0,xmu,dndmu) 
            print*,"end of lambda0 calculation"
            stop
         else
            xmu=xmuold
         endif
      rewind(1)
      endif

C*********  starts solving the problem  *******************************
         print*,' '
         print*,'---------------------------------------',
     &'----------------------------'
         print*,'***************************************',
     &'****************************'

      ! reads in the file "input"
      read(1,*)xxx ! Correlations
      read(1,*)U, ratioJU
      read(1,*)xxx !Inverse Temp, total population (if fixpop=.true.)
      read(1,*)Beta,densfin
      if (.not.DFTplusSS) then ! it uses an analytic hamiltonian or DOS
         read(1,*)xxx
         read(1,*)(tm(L),L=1,M)
         read(1,*)(epsm(L),L=1,M)
         read(1,*)Vhyb,Vnl

         if(Mb.ne.1) stop "Mb≠1, not implemented for analytic ham"
         print*,' '
         print*,'---------------------------------------',
     &        '----------------------------'
         
         if (usedos) then !uses an analytical DOS (provided in subroutine fermionic) and energy sum
            print*,"using analytic DOS and energy sums..."
         else! it uses an analytical Hamiltonian (provided in subroutine analytic_kHam) and k sum
            print*,"using analytic k-space hamiltonian..."
            call analytic_kHam()
         endif                     

      endif ! not DFTplusSS

      J=ratioJU*U
      
      ! these are the output files
      open (41,file='Z_vs_U_mu', form='formatted', 
     &     status='unknown', access='append')
      open (42,file='corr_vs_U_mu_n', form='formatted', 
     &     status='unknown', access='append')

      
         print*,'---------------------------------------',
     &'----------------------------'
         print*,' '
         print*,'U= ',U,'   J= ',J
         print*,'ntot= ',densfin
         if(zeroTemp)then
            print*,"zero Temperature"
            else
               print*,'Beta= ',Beta
         endif

         ! here we plot the convergence test as a function of the iterations
         open(33,file='test',form='formatted',status='unknown')
      
      do L=1,M
         lamold(L)=lam(L)
         hh(L)=lam(L)-lam0(L)+xmu
         S(L)=sqrt(Z(L,L))
      enddo

C**********Starts the external loop **************************************
      nit=0
      
      print*,'iterations'
      print111,' '
      
C*****Calculation of hm's (fermionic problem) from seed for nit=1 ***************
 100  nit=nit+1
      print111,']'
      
      do L1=1,M
         lamold2(L1)=lam(L1) 
         do L2=1,M
            Zold(L1,L2)=Z(L1,L2) 
         enddo
      enddo

      call Fermionic(Z,nf,ffH,hh,emin,emax)
      
      !!diagonal effective fields
      do L=1,M
         hm(L,L)=0.d0
         do LL=1,M
            hm(L,L)=hm(L,L)+S(LL)*ffH(L,LL)
        enddo
      if (verb) then
         print*,L,"lam",lam(L),"  S",S(L),"  hm",hm(L,L),"nf",nf(L,L)
      endif
      enddo
      !!off-diagonal effective fields
      do L1=1,M
          do L2=1,M
             if (L1.ne.L2) hm(L1,L2)=nf(L1,L2) 
         enddo
      enddo

******** mixing of the hm ************
      if (verb) print*,"mixh=",mixh
      do L1=1,M
         if (verb) print*,"L,hm,hmold",L1,hm(L1,L1),hmold(L1,L1)
          do L2=1,M
            if (nit.gt.1) then 
               hm(L1,L2)=hm(L1,L2)*mixh+hmold(L1,L2)*(1.d0-mixh)
            endif
            hmold(L1,L2)=hm(L1,L2)
         enddo
      enddo

      do L=1,M ! this gauge of the slave-spin mean-field depends on the orbital density
         cm(L)=1/sqrt(nf(L,L)-nf(L,L)**2)-1
      enddo
               
ccccccccccccccccc internal iterations only on lambda's
         nn=0
 200     if (verb) print*,""
         do L=1,M
            lamold(L)=lam(L) 
            if (verb) print*,L,"hm",hm(L,L),"lam",lam(L),"cm",cm(L)
         enddo
         nn=nn+1

C*****Diagonalization of the Spin hamiltonian****************************
         call Diagonalization(hm,lam,cm)
         
C***********************************************************************
C*****Calculation of the average values of Sx(L) and Sz(L)  ************

         call AverageValues(Sz,Sx,SxSx,S,SdS,SzSz,cm)

         do L1=1,M
            do L2=1,M
               if (L1.eq.L2) then
                  Z(L1,L2)=S(L1)**2
c                  print*,"S=",S(L1)
               else
                  Z(L1,L2)=Sds(L1,L2)
               endif
            enddo
         enddo
         
         do L=1,M
            nm(L)=0.5d0+Sz(L)
         enddo

C***** Fermionic problem - diagonalization and find the         ********
C***** of Lambdas that satisfy the constraint equations         ********
      if ((.not.DFTplusSS).and.((abs(Vhyb)+abs(Vnl)).lt.tiny)) then 
cc!     separated conjugated search for each orbital
         do L=1,M
            Lorb=L
            x(1)=hh(L)
         call broydn(x,1,check)
c         if (check.eqv..true.) stop "problem in separate broyd"
            hh(L)=x(1)
         enddo
      else ! common conjugated gradient search
         do L=1,M
            x(L)=hh(L)
         enddo
         call broydn(x,M,check)
c         if (check.eqv..true.) stop "problem in common broyd"
         do L=1,M
            hh(L)=x(L)
         enddo
      endif

      do L=1,M
         lam(L)=hh(L)+lam0(L)-xmu
      enddo

      do L=1,M
         testl(L)=abs(lam(L)-lamold(L))
         if (verb) print*,"Lam inside",L,lam(L),lamold(L)
         lam(L)=lam(L)*mixlam+lamold(L)*(1.d0-mixlam) ! **mixing of the lam's
      enddo
c      write(34,*)(testl(L),L=1,M)
c      call flush(34)
      do L=1,M
         if ((testl(L).gt.test_tol)
     &        .and.(nn.lt.nlam)) then
            goto 200
         endif
      enddo      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc adjustment of the chemical potential to get the desired final population cc
      if (fixpop) then
         dens=0.d0              ! total density, based on the spin hamiltonian
         do L=1,M
            dens=dens+2*nm(L)
         enddo
         if (nit.gt.2) then 
c            print*,"calcolo di dndmu: dens, densold", dens, densold
c            print*,"                 xmu, xmuold",xmu,xmuold
            dndmu=(dens-densold)/(xmu-xmuold)
         endif
         if (dndmu.lt.compmin) then 
            print*,"compressibility is too low - put by hand to",compmin
            dndmu=dndmu/abs(dndmu)*compmin
         endif
c         print*, "old mu, compress ",xmu,dndmu
         densold=dens; xmuold=xmu
         xmu=xmu-(dens-densfin)/dndmu
c         print*, "new mu = ",xmu
      do L=1,M
         lam(L)=hh(L)+lam0(L)-xmu
         lam(L)=mixlam*lam(L)+(1.d0-mixlam)*lamold2(L) ! mixing
         do LL=1,M
            Z(L,LL)=mixZ*Z(L,LL)+(1.d0-mixZ)*Zold(L,LL)
         enddo
      enddo
      endif

c      call writeseed(Z,lam,lam0,xmu,dndmu)
c      call flush (11)

C***********************************************************************
C*********************test of self-consistency**************************
      do L=1,M
         testlam(L)=abs(lam(L)-lamold2(L))
         testZ(L)=abs(Z(L,L)-Zold(L,L))
         testdens=abs(dens-densfin)
      enddo
      write(33,*)(testlam(L),testZ(L),L=1,M),testdens
      call flush(33)
      do L=1,M
         if (testlam(L).gt.test_tol) goto 100
         if (testZ(L).gt.test_tol) goto 100
      enddo
      if ((fixpop).and.
     &     (testdens.gt.dens_tol))then
         goto 100
      endif
      ! discard the last change in mu
      if (fixpop) then
         do L=1,M
            lam(L)=hh(L)+lam0(L)-xmuold
         enddo
      endif

      print*,' '
      
C***********************************************************************
C*******calculation of physical quantities******************************
      print*,'___________________________________'
      totn=0.d0
      do L=1,M
         print*,'Z(',L,')=',Z(L,L),'  n(',L,')=',nm(L)
      totn=totn+nm(L)*2
      enddo
      print*,'ntot=',totn
      print*,'___________________________________'
      print*,' '

cccccc here we overwrite the seed with the converged one
      call writeseed(Z,lam,lam0,xmu,dndmu)
     
      do L=1,M
c         write(20+L,*)xmu,U,Z(L,L),hm(L,L),lam(L)
c     &        ,nm(L),totn,dndmu,SzSz(L,L,1),SzSz(L,L,2)
c         call flush (20+L)
c      write(50+L) U,J,xmu,totn,nm(L),Z(L,L),lam(L),lam0(L),hm(L,L),
c     &        (ffh(L,LL),LL=1,M)
c      call flush (50+L)
      enddo

      write(41,*)U,J,xmu,totn,(Z(L,L),L=1,M),(nm(L),L=1,M)
      call flush(41)
      write(42,*)U,J,xmu,totn,(Sz(L),L=1,M),((SzSz(L1,L2,1),
     &     SzSz(L1,L2,2),L2=L1,M),L1=1,M)
      call flush(42)

 111  format(a1,$)
 333  format(a3,I1,a1,$)
      
      stop
      end

C***********************************************************************
C***********************************************************************

      Subroutine Diagonalization(hm,lam,cm)

      implicit double precision (a-h, o-z)
      logical want_to_see_the_matrix
      integer M,totst
      logical Hund_Ising,zeroTemp, usedos, no_ins_seed, DFTplusSS,fixpop
      include 'param.dat'
      integer ket(2*M),bra(2*M)
      real*8 J
      real*8 epsm(M),tm(M)
      real*8 lagr,lam(M),Sorb(M),hm(M,M),cm(M)
      real*8 H(totst,totst),EVEC(totst,totst)
      real*8 E(totst)
      real*8 WORK(3*totst-1)
      common /diag/ E,EVEC
      common /init/ tm,epsm,U,J,Beta,xmu,Vhyb,Vnl

c*****************writing matrix elements********************************
c      print*,"in subroutine diag, U,hm,lam,J"
c      print*, U,(hm(L),L=1,M),(lam(L),L=1,M),J
      do II=1,totst
         do JJ=1,totst
            H(II,JJ)=0.d0
            do L=1,2*M
               ket(L)=mod(int((II-1)/2**(L-1)),2)
               bra(L)=mod(int((JJ-1)/2**(L-1)),2)
            enddo
            if (II.eq.JJ) then ! --------- diagonal part
               sum=0.d0
               do L=1,2*M
                  sum=sum+ket(L)-0.5
               enddo
               sum_up=0.d0
               sum_down=0.d0
               do L=1,M
                  sum_up=sum_up+ket(2*L-1)-0.5
                  sum_down=sum_down+ket(2*L)-0.5
               enddo
               H(II,JJ)=0.5*(U-2*J)*sum*sum
               H(II,JJ)=H(II,JJ)-0.5*J*(sum_up*sum_up+sum_down*sum_down) 
               sum=0.d0
               lagr=0.d0
               do L=1,M
                  Sorb(L)=ket(2*L-1)+ket(2*L)-1.d0
                  sum=sum+Sorb(L)**2
                  lagr=lagr+lam(L)*(Sorb(L)+1)
               enddo
               H(II,JJ)=H(II,JJ)+lagr+J*sum
c               print*,sum, lagr,Sorb(1),Sorb(2) 
            else  ! --------------- non-diagonal part
               diff=0.d0
               diff2=0.d0
               diff_up=0.d0
               diff_do=0.d0
               icount=0
               do L=1,2*M
                  diff=diff+abs(ket(L)-bra(L))
                  diff2=diff2+ket(L)-bra(L)
                  if (int((L+1)/2).eq.(L+1)/2.d0) then
                     diff_up=diff_up+abs(ket(L)-bra(L))
                  else
                     diff_do=diff_do+abs(ket(L)-bra(L))
                  endif
                  if (abs(ket(L)-bra(L)).eq.1) then 
                     icount=icount+1
                     if (icount.eq.1) index=int((L+1)/2)
                     if (icount.eq.2) index2=int((L+1)/2)
                  endif
               enddo
c               print*,II,JJ,diff,diff_up,diff_do
               if (diff.eq.1) then 
                  if (diff2.eq.1) then 
                     H(II,JJ)=hm(index,index)*(cm(index)+1.d0)
                     else if (diff2.eq.-1) then 
                     H(II,JJ)=hm(index,index)*(cm(index)+1.d0)
                  endif           
               endif           

cccccccccccccccccc  local hybridization terms (SxSx) cccccccccccccccccc
                if (diff.eq.2)then
                  if ((diff_up.eq.2).or.(diff_do.eq.2)) then
                     H(II,JJ)=H(II,JJ)+2*hm(index,index2)
                  endif
               endif

               if (Hund_Ising.eqv.(.false.)) then 
cccccccccccccccccc spin flip term cccccccccccccccccccccccccccccccccccccc
               do L1=1,M-1
                  do L2=L1+1,M
                     If (((ket(2*L1-1).eq.1).and.(ket(2*L1).eq.0).and.
     &                    (ket(2*L2-1).eq.0).and.(ket(2*L2).eq.1).and.
     &                    (bra(2*L1-1).eq.0).and.(bra(2*L1).eq.1).and.
     &                    (bra(2*L2-1).eq.1).and.(bra(2*L2).eq.0)).or.
     &                    ((ket(2*L1-1).eq.0).and.(ket(2*L1).eq.1).and.
     &                    (ket(2*L2-1).eq.1).and.(ket(2*L2).eq.0).and.
     &                    (bra(2*L1-1).eq.1).and.(bra(2*L1).eq.0).and.
     &                  (bra(2*L2-1).eq.0).and.(bra(2*L2).eq.1))) then
                        iflag=1
                        do L3=1,M
                           if ((L3.ne.L1).and.(L3.ne.L2)) then
                              if ((ket(2*L3-1).ne.bra(2*L3-1))
     &                             .or.(ket(2*L3).ne.bra(2*L3))) then
                                 iflag=0
                                 goto 101
                              endif
                           endif
                        enddo
 101                    if (iflag.eq.1)     H(II,JJ)=H(II,JJ)-J
                     endif
                  enddo
               enddo
ccccccccccccccccc pair hopping term cccccccccccccccccccccccccccccccccccc
               do L1=1,M-1
                  do L2=L1+1,M
                     if (((ket(2*L1-1).eq.1).and.(ket(2*L1).eq.1).and.
     &                    (ket(2*L2-1).eq.0).and.(ket(2*L2).eq.0).and.
     &                    (bra(2*L1-1).eq.0).and.(bra(2*L1).eq.0).and.
     &                    (bra(2*L2-1).eq.1).and.(bra(2*L2).eq.1)).or.
     &                    ((ket(2*L1-1).eq.0).and.(ket(2*L1).eq.0).and.
     &                    (ket(2*L2-1).eq.1).and.(ket(2*L2).eq.1).and.
     &                    (bra(2*L1-1).eq.1).and.(bra(2*L1).eq.1).and.
     &                  (bra(2*L2-1).eq.0).and.(bra(2*L2).eq.0))) then
                        iflag=1
                        do L3=1,M
                           if ((L3.ne.L1).and.(L3.ne.L2)) then
                              if ((ket(2*L3-1).ne.bra(2*L3-1))
     &                             .or.(ket(2*L3).ne.bra(2*L3))) then
                                 iflag=0
                                 goto 102
                              endif
                           endif
                        enddo
 102                    if (iflag.eq.1)     H(II,JJ)=H(II,JJ)-J
                     endif
                  enddo
               enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               endif !Hund_Ising (Hund_Ising = .false. e' il normale vertice SU(2))

            endif            
         enddo
       enddo

C********************check of the hamiltonian matrix********************
      want_to_see_the_matrix=.false.
      if (want_to_see_the_matrix) then
         print*,''
         do II=1,totst
            do JJ=1,totst
               print333,H(II,JJ)
               if (mod(jj,2).eq.0)  then
                  print334,'|'
               endif
               print334,' '
            enddo
            print*,' '
            if (mod(II,2).eq.0)  then
               print*,'-------------------------
     &              ------------------------------'
            endif
         enddo
         
 333  format ((f5.2,f5.2),$)
 334  format (a2,$)

      endif

C****************diagonalization*****************************************

      call dsyev('V','U',totst,H,totst,E,WORK,3*totst-1,INFO)

      if (INFO.ne.0) then 
         print*,'INFO=',INFO
         stop 'error in the diagonalization routine'
      endif

      do II=1,totst
         do JJ=1,totst
            EVEC(II,JJ)=H(II,JJ)
            if (want_to_see_the_matrix) then
               print*,'EVEC(',II,',',JJ,')=',evec(II,JJ)
            endif
         enddo
      enddo

      return
      end

C************************************************************************
C************************************************************************

      subroutine AverageValues(Sz,Sx,SxSx,S,SdS,SzSz,cm)

      implicit double precision (a-h, o-z)
      integer M,totst
      logical Hund_Ising,zeroTemp, usedos, no_ins_seed, DFTplusSS,fixpop
      include 'param.dat'
      real*8 J,epsm(M),tm(M)
      real*8 Sz(M),Sx(M),SxSx(M,M),SdS(M,M),SzSz(M,M,2)
      real*8 cm(M),S(M),Xs,Xx,Xz,Xss
      real*8 EVEC(totst,totst)      
      real*8 E(totst)

      common /diag/ E,EVEC

      common /init/ tm,epsm,U,J,Beta,xmu,Vhyb,Vnl

      PartF=0.d0
      do II=1,totst
         if (zeroTemp) then
            PartF=PartF+1
            ndeg=II
            if (abs(E(II+1)-E(1)).gt.tiny) goto 131
         else
            PartF=PartF+dexp(-Beta*(E(II)-E(1)))
         endif
      enddo
 131  continue

                 ! <Sx(L)>, <Sz(L)>
      if (zeroTemp) then
         nfin=ndeg
      else
         nfin=totst
      endif
       do L=1,M
         Sz(L)=0.d0
         Sx(L)=0.d0
         S(L)=0.d0
         do II=1,nfin
            Xz=0.d0
            Xx=0.d0
            Xs=0.d0
            do JJ=1,totst
               occ=mod(int((JJ-1)/2**(2*L-2)),2) !!! 2*L --> paramagnetic case
               if (occ.eq.1) then   
                  JJflip=JJ-2**(2*L-2)
               else
                  JJflip=JJ+2**(2*L-2)
               endif
               Xx=Xx+Evec(JJ,II)*Evec(JJflip,II)*0.5d0
               Xz=Xz+(occ-0.5d0)*Evec(JJ,II)*Evec(JJ,II)
               if (occ.eq.1) then   
               Xs=Xs+Evec(JJ,II)*Evec(JJflip,II)
               else
               Xs=Xs+Evec(JJ,II)*Evec(JJflip,II)*cm(L)
               endif
            enddo
            if (zeroTemp) then
               Sz(L)=Sz(L)+Xz/PartF
               Sx(L)=Sx(L)+Xx/PartF
               S(L)=S(L)+Xs/PartF
            else
               Sz(L)=Sz(L)+Xz*dexp(-Beta*(E(II)-E(1)))/PartF
               Sx(L)=Sx(L)+Xx*dexp(-Beta*(E(II)-E(1)))/PartF
               S(L)=S(L)+Xs*dexp(-Beta*(E(II)-E(1)))/PartF
            endif
         enddo   
      enddo

                ! mean values <Sx(L)Sx(L')> and <S^dag(L)S(L')>
      do L1=1,M
      do L2=1,M
         SxSx(L1,L2)=0.d0
         SdS(L1,L2)=0.d0
         SzSz(L1,L2,1)=0.d0 ! 1 is for up-up correlations
         SzSz(L1,L2,2)=0.d0 ! 2 is for up-down
         do II=1,nfin
            Xx=0.d0
            Xss=0.d0
            Zz=0.d0
            Zzud=0.d0            
            do JJ=1,totst
               occ1=mod(int((JJ-1)/2**(2*L1-2)),2)!!! 2*L --> paramagnetic case
               occ2=mod(int((JJ-1)/2**(2*L2-2)),2)!!! 2*L --> paramagnetic case
               occ2_do=mod(int((JJ-1)/2**(2*L2-1)),2)!!! spin down, for spincorr!
               if (occ1.eq.1) then   
                  JJflip1=JJ-2**(2*L1-2)
                  fac=1
               else
                  JJflip1=JJ+2**(2*L1-2)
                  fac=cm(L1)
               endif
               if (occ2.eq.1) then   
                  JJflip2=JJ-2**(2*L2-2)
                  fac=fac*cm(L2)
               else
                  JJflip2=JJ+2**(2*L2-2)
                  fac=fac
               endif
               Xx=Xx+Evec(JJflip1,II)*Evec(JJflip2,II)*0.25d0 !SxSx!not 2Sx2Sx
               Xss=Xss+Evec(JJflip1,II)*Evec(JJflip2,II)*fac  ! S S !
               Zz=Zz+Evec(JJ,II)*Evec(JJ,II)*
     &              (occ1-0.5d0)*(occ2-0.5d0)
               Zzud=Zzud+Evec(JJ,II)*
     &              Evec(JJ,II)*(occ1-0.5d0)*(occ2_do-0.5d0)
            enddo
            if (zeroTemp) then
               SxSx(L1,L2)=SxSx(L1,L2)+Xx/PartF
               SdS(L1,L2)=SdS(L1,L2)+Xss/PartF
           SzSz(L1,L2,1)=SzSz(L1,L2,1)+Zz/PartF
            SzSz(L1,L2,2)=SzSz(L1,L2,2)+Zzud/PartF           
            else
               SxSx(L1,L2)=SxSx(L1,L2)+Xx*dexp(-Beta*(E(II)-E(1)))/PartF
               SdS(L1,L2)=SdS(L1,L2)+Xss*dexp(-Beta*(E(II)-E(1)))/PartF
           SzSz(L1,L2,1)=SzSz(L1,L2,1)+Zz*dexp(-Beta*(E(II)-E(1)))/PartF
            SzSz(L1,L2,2)=SzSz(L1,L2,2)
     &           +Zzud*dexp(-Beta*(E(II)-E(1)))/PartF           
            endif
        enddo
      enddo
      enddo

ccccccccccccccccc check!

c      if ((abs(SxSx(1,1)-0.25d0).gt.tiny).or.
c     &     (abs(SzSz(1,1,1)-0.25d0).gt.1.d-10).or.
c     & (abs(SxSx(1,2)-SxSx(2,1)).gt.tiny)) then 
c         print*,SxSx(1,1)-0.25d0,SzSz(1,1,1)-0.25d0,SxSx(1,2)-SxSx(2,1)
c         stop "error in the average values"
c      endif
c      print*,"SxSx",SxSx(1,2)

      return
      end

C************************************************************************
C************************************************************************
C     In this routine fermionic average values are calculated, for a 
C     given electronic structure.
C     nf(m,n) contains the local average values <f^+_ims f_ins> 
C     on the diagonal, i.e. the local density for orbital m, and
C     Hamloc(m,n) <f^+_ims f_ins> on the off-diagonal elements.
C     ffH(m,n) is the average value sum_i≠j H^mn_ij <f^+_ims f_jns>
C     involved in the calculation of hm(m)=sum_n <S_n> ffH(m,n)
C************************************************************************

      subroutine Fermionic(Z,nf,ffH,hh,emin,emax)

      implicit double precision (a-h, o-z)
      logical Hund_Ising,zeroTemp, usedos, no_ins_seed, DFTplusSS,fixpop
      include 'param.dat'
      real*8 J, tm(M),epsm(M)
      real*8 Z(M,M),hh(M)
      complex*16 Hamnl(M*Mb,M*Mb)
      real*8 ffH(M,M),nf(M,M)
      real*8 E(M*Mb),fermi(M*Mb)
      real*8 wt(-nkx/2+1:nkx/2,-nky/2+1:nky/2,-nkz/2+1:nkz/2)
      complex*16 Hamk(-nkx/2+1:nkx/2,-nky/2+1:nky/2,
     &     -nkz/2+1:nkz/2,M*Mb,M*Mb)
c      real*8 wt(0:nkx-1,0:nky-1,0:nkz-1)
c      complex*16 Hamk(0:nkx-1,0:nky-1,0:nkz-1,M*Mb,M*Mb)
      complex*16 Hamloc(M*Mb,M*Mb),ev(M*Mb,M*Mb)
      complex*16 WORK(2*M*Mb-1)
      real*8 WORK2(3*M*Mb-2)
      common /init/ tm,epsm,U,J,Beta,xmu,Vhyb,Vnl
      common/ hamilt/ hamloc, hamk, wt

      pi=acos(-1.d0)

      do L1=1,M
         do L2=1,M
            nf(L1,L2)=0.d0
            ffH(L1,L2)=0.d0
         enddo
      enddo

      emin=1000
      emax=-emin

      if ((.not.DFTplusSS).and.(usedos)) then ! using density of states
         de=2.d0/dfloat(ne) ! Dos goes from -1 to 1
      do ie=1,ne
         en=-1.d0+de*ie !energy is then multiplied by tm
         DOS=2/Pi*dsqrt(1-en**2) !Bethe lattice (semicircular DOS of a Caley tree)
         do L1=1,M
            do L2=1,L1
               if (L1.eq.L2) then ! bare non-local hamiltonian
                  Hamnl(L1,L1)=2*tm(L1)*en 
                  Hamloc(L1,L1)=epsm(L1)
               else
                  Hamnl(L1,L2)=2*Vnl*en
                  Hamnl(L2,L1)=Hamnl(L1,L2)
                  Hamloc(L1,L2)=Vhyb
                  Hamloc(L2,L1)=Hamloc(L1,L2)
               endif

               if (L1.eq.L2) then ! full renorm fermionic hamiltonian
                  ev(L1,L1)=Z(L1,L1)*Hamnl(L1,L1)+Hamloc(L1,L1)-hh(L1)
               else
                  ev(L1,L2)=sqrt(Z(L1,L1)*Z(L2,L2))*Hamnl(L1,L2)
     &                 +Z(L1,L2)*Hamloc(L1,L2)
                  ev(L2,L1)=ev(L1,L2)
               endif
            enddo
         enddo
         call zheev('V','U',M,ev,M,E,WORK,2*M-1,WORK2,INFO)
         if (INFO.ne.0) then
            print*,"in Fermionic INFO=",INFO
            stop
         endif
         do L=1,M
            if (E(L).lt.emin) emin=E(L)
            if (E(L).gt.emax) emax=E(L)
            if (zeroTemp) then
               if (E(L).le.0.d0) then
                  fermi(L)=1.d0
               else if ((E(L).gt.0.d0).and.(E(L).le.de)) then
                  fermi(L)=1.d0-E(L)/de ! this ensures a non discretized dependency
               else ! of the integral on hh(L), which can give a singular
                  fermi(L)=0.d0 ! Jacobian matrix in the conj grad routine
               endif
            else
               fermi(L)=1.d0/(dexp(Beta*(E(L)))+1)
            endif
         enddo
         do L=1,M
            do L1=1,M
               do L2=1,M
                  nf(L,L1)=nf(L,L1)+real(ev(L,L2)*conjg(ev(L1,L2))*
     &                 fermi(L2)*DOS*de)
                  ffH(L,L1)=ffH(L,L1)+real(DOS*de*ev(L1,L2)*
     &                 conjg(ev(L,L2))*fermi(L2)*Hamnl(L,L1))
               enddo
            enddo
         enddo
      enddo                     !en

      else                      ! using k-summation
         
         xnorm=0.d0
c sum over the whole brillouin zone
         do ikx=-nkx/2+1,nkx/2
         xk=ikx*(Pi/dfloat(nkx/2))
         do iky=-nky/2+1,nky/2
         yk=iky*(Pi/dfloat(nky/2))
         do ikz=-nkz/2+1,nkz/2
         zk=ikz*(Pi/dfloat(nkz/2))
         weight=1.d0
c sum over a quarter of the Brillouin zone (with weights adjusted for borders)
c         do ikx=0,nkx-1
c         xk=ikx*(Pi/dfloat(nkx))
c         do iky=0,nky-1
c         yk=iky*(Pi/dfloat(nky))
c         do ikz=0,nkz-1
c         zk=ikz*(Pi/dfloat(nkz))
c         weight=wt(ikx,iky,ikz)

         do L1=1,M*Mb
            L1m=mod(L1-1,M)+1
            do L2=1,L1
               L2m=mod(L2-1,M)+1 
                if (L1.eq.L2) then
                  ev(L1,L1)=Hamk(ikx,iky,ikz,L1,L1)*Z(L1m,L1m)
     &                 +Hamloc(L1,L1)-hh(L1m)
               else
                  ev(L1,L2)=sqrt(Z(L1m,L1m)*Z(L2m,L2m))
     &                 *Hamk(ikx,iky,ikz,L1,L2)
     &                 + Z(L1m,L2m)*Hamloc(L1,L2)
                  ev(L2,L1)=conjg(ev(L1,L2))
               endif
            enddo
         enddo
         
         call zheev('V','U',M*Mb,ev,M*Mb,E,WORK,2*M*Mb-1,WORK2,INFO)
         if (INFO.ne.0) then
            print*,"in Fermionic INFO=",INFO
            stop
         endif
         do L=1,M*Mb
            if (E(L).lt.emin) emin=E(L)
            if (E(L).gt.emax) emax=E(L)
            if (zeroTemp) then
               if (E(L).le.0.d0) then
                  fermi(L)=1.d0
               else
                  fermi(L)=0.d0
               endif
            else
               fermi(L)=1.d0/(dexp(Beta*(E(L)))+1)
            endif
         enddo

         do L=1,M
            do L1=1,M*Mb
               L1m=mod(L1-1,M)+1 
               do L2=1,M*Mb
                  if (L1.le.M)  nf(L,L1)=nf(L,L1)
     &                 +real(conjg(ev(L,L2))*ev(L1,L2)*fermi(L2)*weight)
                  ffH(L,L1m)=ffH(L,L1m)+real(conjg(ev(L,L2))*ev(L1,L2)*
     &                 Hamk(ikx,iky,ikz,L,L1)*fermi(L2)*weight)
                  if (L1.gt.M) ffH(L,L1m)=ffH(L,L1m)+real(hamloc(L,L1)*
     &                 conjg(ev(L,L2))*ev(L1,L2)*fermi(L2)*weight)
               enddo
c               print*,ikx,iky,ikz,Hamk(ikx,iky,ikz,L,L1), 
c     &              hamloc(L,L1),weight
            enddo
         enddo
 
         xnorm=xnorm+weight
         enddo       !zk
         enddo      !yk
         enddo     !xk

         do L1=1,M
            do L2=1,M
               nf(L1,L2)=nf(L1,L2)/xnorm
               ffH(L1,L2)= ffH(L1,L2)/xnorm
            enddo
         enddo
         
      endif   ! DOS or k-sum

      do L1=1,M
         do L2=1,M
            if(L1.ne.L2) nf(L1,L2)=nf(L1,L2)*real(Hamloc(L1,L2))
c            print*,nf(L1,L2),ffh(L1,L2)
         enddo
      enddo
c      print*,"nf(1,1),fh(1,1)",nf(1,1),ffh(1,1)
c         print*,""

      return
      end

  
C************************************************************************
C************************************************************************
      subroutine funcv(Mdum,x,f)

      implicit double precision (a-h, o-z)
      logical Hund_Ising,zeroTemp, usedos, no_ins_seed, DFTplusSS,fixpop
      include 'param.dat'
      real*8 J,tm(M),epsm(M)
      real*8 x(M)
      real*8 Z(M,M),hh(M),nm(M),f(M)
      real*8 nf(M,M),ffH(M,M)
      common /newton/ nm,Z
      common /init/ tm,epsm,U,J,Beta,xmu,Vhyb,Vnl
      common / whichorb /Lorb


      pi=acos(-1.d0)

      if ((.not.DFTplusSS).and.((abs(Vhyb)+abs(Vnl)).lt.tiny)) then 
cc         ! separate conj search in absence of hybridization
         hh(Lorb)=x(1)
         call Fermionic(Z,nf,ffH,hh,emin,emax)
cc         print*,Lorb,hh(Lorb),nf(Lorb,Lorb),nm(Lorb)
         f(1)=nf(Lorb,Lorb)-nm(Lorb)
cc            print*,"nf in fermionic", Lorb,nf(Lorb,Lorb),f(1)
cc            print*,"hm", Lorb,hm(Lorb,Lorb)
ccc         write(50+Lorb,*) hh(Lorb),nf(Lorb,Lorb),nm(Lorb)
      else ! common conjugated gradient search with hybridization
         do L=1,M
            hh(L)=x(L)
         enddo
         call Fermionic(Z,nf,ffH,hh,emin,emax)
         
         do L=1,M
            f(L)=nf(L,L)-nm(L)
         enddo
      endif


      return
      end
  
C************************************************************************
C************************************************************************
C This routine saves the k-space Hamiltonian (if usedos=.false.)
C for a user-provided analytical form
C***********************************************************************
      subroutine analytic_kHam()

      implicit double precision (a-h, o-z)
      logical Hund_Ising,zeroTemp, usedos, no_ins_seed, DFTplusSS,fixpop
      include 'param.dat'
      complex*16 Xi,zdum
      real*8 wt(-nkx/2+1:nkx/2,-nky/2+1:nky/2,-nkz/2+1:nkz/2)
      complex*16 Hamk(-nkx/2+1:nkx/2,-nky/2+1:nky/2,
     &     -nkz/2+1:nkz/2,M*Mb,M*Mb)
c      real*8 wt(0:nkx,0:nky,0:nkz)
c      complex*16 Hamk(0:nkx,0:nky,0:nkz,M*Mb,M*Mb)
      complex*16 Hamloc(M*Mb,M*Mb)
      real*8 J,tm(M),epsm(M)

      common /init/ tm,epsm,U,J,Beta,xmu,Vhyb,Vnl
      common/ hamilt/ hamloc, hamk, wt

      Xi=cmplx(0.d0,1.d0)
      pi=dacos(-1.d0)

      do L1=1,M
         do L2=L1,M
            if (L1.eq.L2) then
               hamloc(L1,L2)=epsm(L1) ! here it reads the local part from input (diagonal)
            else
               hamloc(L1,L2)=Vhyb  !here the off-diagonal local part
               hamloc(L2,L1)=hamloc(L1,L2)
            endif
         enddo
      enddo


c sum over the whole Brillouin zone
      do ikx=-nkx/2+1,nkx/2
      xk=ikx*(Pi/dfloat(nkx/2))
      do iky=-nky/2+1,nky/2
      yk=iky*(Pi/dfloat(nky/2))
      do ikz=-nkz/2+1,nkz/2
      zk=ikz*(Pi/dfloat(nkz/2))
c sum over a quarter of the Brillouin zone (incorrect weights at the borders!)
c         do ikx=0,nkx
c         xk=ikx*(Pi/dfloat(nkx))
c         do iky=0,nky
c         yk=iky*(Pi/dfloat(nky))
c         do ikz=0,nkz
c         zk=ikz*(Pi/dfloat(nkz))

         do L1=1,M
            do L2=1,M
               if (L1.eq.L2) then
                  epsk=2*tm(L1)*(cos(xk)+cos(yk)+cos(zk)) ! here we input the dispersion
                  Hamk(ikx,iky,ikz,L1,L1)=epsk
               else
                  Hamk(ikx,iky,ikz,L1,L2)=2*sqrt(3.d0)*Vnl* 
     &                 (cos(xk)-cos(yk))*cos(zk) ! here the off-diagonal part
                  Hamk(ikx,iky,ikz,L2,L1)=
     &                 conjg(Hamk(ikx,iky,ikz,L1,L2))
               endif
            enddo
         enddo

ccc   weights for the integration of a quarter of the BZ (dummy for full integral)
         weight=1.d0
         if (ikx.eq.0) weight=weight*0.5
         if (iky.eq.0) weight=weight*0.5
         if (ikz.eq.0) weight=weight*0.5
         wt(ikx,iky,ikz)=weight
         
      enddo                     !kz
      enddo                     !ky
      enddo                     !kx

      return
      end
C************************************************************************
C************************************************************************
C This routine constructs the k-space Hamiltonian (if DFTplusSS=.true.)
C for a DFT real-space hamiltonian provided in the file '[prefix].rHam',
C and writes it in a file with name [prefix].kHam
C***********************************************************************
      subroutine calc_kHAM(prefix)

      implicit double precision (a-h, o-z)
      logical Hund_Ising,zeroTemp, usedos, no_ins_seed, DFTplusSS,fixpop
      include 'param.dat'
      character(*) prefix
      integer ix(neltotmax),iy(neltotmax),iz(neltotmax)
      integer ior(neltotmax),jor(neltotmax)
      complex*16 ham(neltotmax)
      complex*16 Xi,zdum
      real*8 wt(-nkx/2+1:nkx/2,-nky/2+1:nky/2,-nkz/2+1:nkz/2)
      complex*16 Hamk(-nkx/2+1:nkx/2,-nky/2+1:nky/2,
     &     -nkz/2+1:nkz/2,M*Mb,M*Mb)
c      real*8 wt(0:nkx-1,0:nky-1,0:nkz-1)
c      complex*16 Hamk(0:nkx-1,0:nky-1,0:nkz-1,M*Mb,M*Mb)
      complex*16 Hamloc(M*Mb,M*Mb)
      real*8 J,tm(M),epsm(M)

      common /init/ tm,epsm,U,J,Beta,xmu,Vhyb,Vnl
      common/ hamilt/ hamloc, hamk, wt

      Xi=cmplx(0.d0,1.d0)
      pi=dacos(-1.d0)

      open(3,file=prefix//".kHam",form='unformatted',status="unknown")


      print*,' '
      print*,'---------------------------------------',
     &'----------------------------'
      print*,"calculating the k-space DFT hamiltonian..."
          ! input
      open(12,file=prefix//'.rHam',form='formatted', status='old')
          ! hamiltonian is the matrix of hopping between 0,0,0 and i,j,k
          ! for each orbital (M x M matrix)
         
      icount=0
 100  read(12,*,err=101,end=102)ixc,iyc,izc,iorb,jorb,dum1,dum2
c      print*,ixc,iyc,izc,iorb,jorb,dum1,dum2
      icount=icount+1
      if (icount.gt.neltotmax) stop "neltotmax too small"
      ix(icount)=ixc
      iy(icount)=iyc
      iz(icount)=izc
      ior(icount)=iorb
      jor(icount)=jorb
      if ((iorb.gt.M*Mb).or.(iorb.gt.M*Mb)) then
         print*,"error, Mb too small?"
         goto 101
      endif
      ham(icount)=cmplx(dum1,dum2)

         ! we here identify the local orbital energies
      if ((ixc.eq.0).and.(iyc.eq.0).and.(izc.eq.0)) then
         if (iorb.eq.jorb) then
            epsm(iorb)=dum1
         endif
         Hamloc(iorb,jorb)=dum1
      endif
      
      goto 100
 101  stop "error in reding the input DFT Hamiltonian"
 102  continue
      neltot=icount
         

c sum over the whole Brillouin zone
      do ikx=-nkx/2+1,nkx/2
      xk=ikx*(Pi/dfloat(nkx/2))
      do iky=-nky/2+1,nky/2
      yk=iky*(Pi/dfloat(nky/2))
      do ikz=-nkz/2+1,nkz/2
      zk=ikz*(Pi/dfloat(nkz/2))
c sum over a quarter of the Brillouin zone (incorrect weights at the borders!)
c         do ikx=0,nkx-1
c         xk=ikx*(Pi/dfloat(nkx))
c         do iky=0,nky-1
c         yk=iky*(Pi/dfloat(nky))
c         do ikz=0,nkz-1
c         zk=ikz*(Pi/dfloat(nkz))

      do L1=1,M*Mb
         do L2=1,M*Mb
            hamk(ikx,iky,ikz,L1,L2)=cmplx(0.d0,0.d0)
         enddo
      enddo
      do ic=1,neltot                     
         zdum=exp(Xi*(xk*ix(ic)+yk*iy(ic)+zk*iz(ic)))*ham(ic)
         if (.not.((ix(ic).eq.0).and.(iy(ic).eq.0)
     &        .and.(iz(ic).eq.0))) then
            hamk(ikx,iky,ikz,ior(ic),jor(ic))=
     &           hamk(ikx,iky,ikz,ior(ic),jor(ic))+zdum
         endif
      enddo
      
ccc weights for the integration of a quarter of the BZ (dummy for full integral)
      weight=1.d0
c      if (ikx.eq.0) weight=weight*0.5
c      if (iky.eq.0) weight=weight*0.5
c      if (ikz.eq.0) weight=weight*0.5
      wt(ikx,iky,ikz)=weight
      do L1=1,M*Mb
         do L2=1,M*Mb   
            write(3)ikx,iky,ikz,L1,L2,Hamk(ikx,iky,ikz,L1,L2)
     &           ,Hamloc(L1,L2),wt(ikx,iky,ikz)
         enddo
      enddo

      enddo                     !kz
      enddo                     !ky
      enddo                     !kx

      print*,"--------- end of kHam calculation -------------"
      stop

      return
      end

C************************************************************************
C************************************************************************
C This routine reads the k-space Hamiltonian (if usedos=.false.)
C from file with name [prefix].kHam
C***********************************************************************
      subroutine read_kHAM(prefix)

      implicit double precision (a-h, o-z)
      logical Hund_Ising,zeroTemp, usedos, no_ins_seed, DFTplusSS,fixpop
      include 'param.dat'
      character(*) prefix
      integer ix(neltotmax),iy(neltotmax),iz(neltotmax)
      integer ior(neltotmax),jor(neltotmax)
      complex*16 ham(neltotmax)
      complex*16 Xi,zdum
      real*8 wt(-nkx/2+1:nkx/2,-nky/2+1:nky/2,-nkz/2+1:nkz/2)
      complex*16 Hamk(-nkx/2+1:nkx/2,-nky/2+1:nky/2,
     &     -nkz/2+1:nkz/2,M*Mb,M*Mb)
c      real*8 wt(0:nkx-1,0:nky-1,0:nkz-1)
c      complex*16 Hamk(0:nkx-1,0:nky-1,0:nkz-1,M*Mb,M*Mb)
      complex*16 Hamloc(M*Mb,M*Mb) 
      real*8 J,tm(M),epsm(M)

      common /init/ tm,epsm,U,J,Beta,xmu,Vhyb,Vnl
      common/ hamilt/ hamloc, hamk, wt

      Xi=cmplx(0.d0,1.d0)
      pi=dacos(-1.d0)

      open(3,file=prefix//".kHam",form='unformatted',status="old")

      print*,' '
      print*,'---------------------------------------',
     &     '----------------------------'
      print*,"reading the DFT hamiltonian..."
! input
      icount=0
 100  read(3,end=102, err=101)ikx,iky,ikz,L1,L2,
     &     Hamk(ikx,iky,ikz,L1,L2),Hamloc(L1,L2),wt(ikx,iky,ikz)
      goto 100

 101  stop "problem in reading kham"
 102  return
      end

C***********************************************************************
      subroutine writeseed(Z,lam,lam0,xmu,dndmu)

      implicit double precision (a-h, o-z)
      logical Hund_Ising,zeroTemp, usedos, no_ins_seed, DFTplusSS,fixpop
      include 'param.dat'
      real*8 Z(M,M)
      real*8 lam(M),lam0(M)


      rewind (11)
      write(11,*) "********** Z(L1,L2) *********** "
      do L1=1,M
         do L2=1,L1-1
            write(11,112) '         ' ! print L1-1 blank spaces
         enddo
         write (11,*)(Z(L1,L2),L2=L1,M)
      enddo
      write(11,*) "" 
      write(11,*) "*********** lambda(L),  lambda0(L) ***********"
      do L1=1,M
         write(11,*)lam(L1),lam0(L1)
      enddo

      write(11,*) "" 
      write(11,*) " *********** chemical potential, 
     &    compressibility  ***********"
      write(11,*) xmu, dndmu
      
 112  format(a26,$)

      return
      end


C***********************************************************************
      subroutine calclambdazero(densfin,lam,lam0)

      implicit double precision (a-h, o-z)
      logical Hund_Ising,zeroTemp, usedos, no_ins_seed, DFTplusSS,fixpop
      logical check
      include 'param.dat'
      real*8 J,tm(M),epsm(M)
      real*8 Z(M,M)
      real*8 lam(M),lam0(M),hm(M,M),cm(M),ffh(M,M),nf(M,M)
      real*8 hh(M),x(M)

      common /newton2/ nf,hm,cm,ffh
      common /init/ tm,epsm,U,J,Beta,xmu,Vhyb,Vnl

      U=0.d0
      J=0.d0

      do L1=1,M  
         do L2=1,M
            Z(L1,L2)=1.d0
         enddo
         hh(L1)=0.d0
      enddo

c finds minimum and maximum energy for the bisection procedure
      call Fermionic(Z,nf,ffH,hh,emin,emax)
      xmumax=emax+0.1d0
      xmumin=emin-0.1d0
      print*,"emin,emax = ",emin,emax
      
c      start bisection for xmu

 10   xmu=(xmumin+xmumax)*0.5d0
      print*,"mu=",xmu
      do L=1,M  
         hh(L)=xmu
      enddo
      call Fermionic(Z,nf,ffH,hh,emin,emax)
      totn=0.d0
      do L=1,M  
         totn=totn+2*nf(L,L)
      enddo

      if (abs(totn-densfin).lt.1.d-05) goto 20
      if (totn.lt.densfin) xmumin=xmu
      if (totn.ge.densfin) xmumax=xmu
      goto 10
 20   continue

      print*
      print*,"orb    population"
      do L=1,M  
        print*, L, nf(L,L)
      enddo
      
      !!diagonal effective fields
      do L=1,M
         hm(L,L)=0.d0
         do LL=1,M
            hm(L,L)=hm(L,L)+ffH(L,LL) !! at U=0 <S>=1 by constr. (hm<0)
         enddo
         print*,"hm",L,hm(L,L),"ffh",(ffH(L,LL),LL=1,M)
      enddo
      !!off-diagonal effective fields
      do L1=1,M
          do L2=1,M
             if (L1.ne.L2) hm(L1,L2)=nf(L1,L2) 
         enddo
      enddo

      do L=1,M
         cm(L)=1/sqrt(nf(L,L)-nf(L,L)**2)-1
      enddo
               
ccccccccccccccccc internal iterations only on lambdas
         do L=1,M
            x(L)=lam(L)
         enddo
         call broydn2(x,M,check)
         do L=1,M ! lambda0's are the converged lambda's 
            lam0(L)=x(L)
         enddo

         return
         end

C************************************************************************
C For the conjugated gradient search of the lambda0's
C************************************************************************
      subroutine funcv2(Mdum,x,f)

      implicit double precision (a-h, o-z)
      logical Hund_Ising,zeroTemp, usedos, no_ins_seed, DFTplusSS,fixpop
      include 'param.dat'
      real*8 J,tm(M),epsm(M)
      real*8 x(M)
      real*8 nm(M),f(M)
      real*8 Sx(M),Sz(M),SxSx(M,M),S(M),SdS(M,M),SzSz(M,M,2)
      real*8 lam(M),lam0(M),hm(M,M),cm(M),nf(M,M),ffh(M,M)
      common /newton2/ nf,hm,cm,ffh
      common /init/ tm,epsm,U,J,Beta,xmu,Vhyb,Vnl
      common / whichorb /Lorb

      do L=1,M
         lam(L)=x(L)
      enddo

C*****Diagonalization of the Spin hamiltonian****************************
      call Diagonalization(hm,lam,cm)
C***********************************************************************
C*****Calculation of the average values of Sx(L) and Sz(L)  ************

      call AverageValues(Sz,Sx,SxSx,S,SdS,SzSz,cm)

      do L=1,M
c      print*,"S",L,S(L)
         nm(L)=0.5d0+Sz(L)
         f(L)=nf(L,L)-nm(L)
c         print*,nf(L,L),nm(L),nf(L,L)-nm(L)
c         hm(L,L)=0.d0
cc         do LL=1,M
c            hm(L,L)=hm(L,L)+S(LL)*ffH(L,LL)
c         enddo
c         print*,"hm",L,hm(L,L),"ffh",(ffH(L,LL),LL=1,M)
      enddo

      return
      end
  
C************************************************************************
