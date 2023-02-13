      SUBROUTINE broydn(x,n,check)
      INTEGER n,nn,NP,MAXITS
      REAL*8 x(n),fvec,EPS,TOLF,TOLMIN,TOLX,STPMX
      LOGICAL check
      PARAMETER (NP=40,MAXITS=200,EPS=1.e-7,TOLF=1.e-4,TOLMIN=1.e-6,
     *TOLX=EPS,STPMX=100.)
      COMMON /newtv/ fvec(NP),nn
CU    USES fdjac,fmin,lnsrch,qrdcmp,qrupdt,rsolv
      INTEGER i,its,j,k
      REAL*8 den,f,fold,stpmax,sum,temp,test,c(NP),d(NP),fvcold(NP),
     *g(NP),p(NP),qt(NP,NP),r(NP,NP),s(NP),t(NP),w(NP),xold(NP),fmin
      LOGICAL restrt,sing,skip
      EXTERNAL fmin
      nn=n
      f=fmin(x)
      test=0.
      do 11 i=1,n
        if(abs(fvec(i)).gt.test)test=abs(fvec(i))
11    continue
      if(test.lt..01*TOLF)return
      sum=0.
      do 12 i=1,n
        sum=sum+x(i)**2
12    continue
      stpmax=STPMX*max(sqrt(sum),float(n))
      restrt=.true.
      do 44 its=1,MAXITS
        if(restrt)then
          call fdjac(n,x,fvec,NP,r)
          call qrdcmp(r,n,NP,c,d,sing)
          if(sing) stop 'singular Jacobian in broydn'
          do 14 i=1,n
            do 13 j=1,n
              qt(i,j)=0.
13          continue
            qt(i,i)=1.
14        continue
          do 18 k=1,n-1
            if(c(k).ne.0.)then
              do 17 j=1,n
                sum=0.
                do 15 i=k,n
                  sum=sum+r(i,k)*qt(i,j)
15              continue
                sum=sum/c(k)
                do 16 i=k,n
                  qt(i,j)=qt(i,j)-sum*r(i,k)
16              continue
17            continue
            endif
18        continue
          do 21 i=1,n
            r(i,i)=d(i)
            do 19 j=1,i-1
              r(i,j)=0.
19          continue
21        continue
        else
          do 22 i=1,n
            s(i)=x(i)-xold(i)
22        continue
          do 24 i=1,n
            sum=0.
            do 23 j=i,n
              sum=sum+r(i,j)*s(j)
23          continue
            t(i)=sum
24        continue
          skip=.true.
          do 26 i=1,n
            sum=0.
            do 25 j=1,n
              sum=sum+qt(j,i)*t(j)
25          continue
            w(i)=fvec(i)-fvcold(i)-sum
            if(abs(w(i)).ge.EPS*(abs(fvec(i))+abs(fvcold(i))))then
              skip=.false.
            else
              w(i)=0.
            endif
26        continue
          if(.not.skip)then
            do 28 i=1,n
              sum=0.
              do 27 j=1,n
                sum=sum+qt(i,j)*w(j)
27            continue
              t(i)=sum
28          continue
            den=0.
            do 29 i=1,n
              den=den+s(i)**2
29          continue
            do 31 i=1,n
              s(i)=s(i)/den
31          continue
            call qrupdt(r,qt,n,NP,t,s)
            do 32 i=1,n
              if(r(i,i).eq.0.) stop 'r singular in broydn'
              d(i)=r(i,i)
32          continue
          endif
        endif
        do 34 i=1,n
          sum=0.
          do 33 j=1,n
            sum=sum+qt(i,j)*fvec(j)
33        continue
          g(i)=sum
34      continue
        do 36 i=n,1,-1
          sum=0.
          do 35 j=1,i
            sum=sum+r(j,i)*g(j)
35        continue
          g(i)=sum
36      continue
        do 37 i=1,n
          xold(i)=x(i)
          fvcold(i)=fvec(i)
37      continue
        fold=f
        do 39 i=1,n
          sum=0.
          do 38 j=1,n
            sum=sum+qt(i,j)*fvec(j)
38        continue
          p(i)=-sum
39      continue
        call rsolv(r,n,NP,d,p)
        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
        test=0.
        do 41 i=1,n
          if(abs(fvec(i)).gt.test)test=abs(fvec(i))
41      continue
        if(test.lt.TOLF)then
          check=.false.
          return
        endif
        if(check)then
          if(restrt)then
            return
          else
            test=0.
            den=max(f,.5*n)
            do 42 i=1,n
              temp=abs(g(i))*max(abs(x(i)),1.)/den
              if(temp.gt.test)test=temp
42          continue
            if(test.lt.TOLMIN)then
              return
            else
              restrt=.true.
            endif
          endif
        else
          restrt=.false.
          test=0.
          do 43 i=1,n
            temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
            if(temp.gt.test)test=temp
43        continue
          if(test.lt.TOLX)return
        endif
44    continue
      stop 'MAXITS exceeded in broydn'
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE fdjac(n,x,fvec,np,df)
      INTEGER n,np,NMAX
      REAL*8 df(np,np),fvec(n),x(n),EPS
      PARAMETER (NMAX=40,EPS=1.e-4)
CU    USES funcv
      INTEGER i,j
      REAL*8 h,temp,f(NMAX)
      do 12 j=1,n
        temp=x(j)
        h=EPS*abs(temp)
        if(h.eq.0.)h=EPS
        x(j)=temp+h
        h=x(j)-temp
        call funcv(n,x,f)
        x(j)=temp
        do 11 i=1,n
          df(i,j)=(f(i)-fvec(i))/h
11      continue
12    continue
      return
      END
      FUNCTION fmin(x)
      INTEGER n,NP
      REAL*8 fmin,x(*),fvec
      PARAMETER (NP=40)
      COMMON /newtv/ fvec(NP),n
      SAVE /newtv/
CU    USES funcv
      INTEGER i
      REAL*8 sum
      call funcv(n,x,fvec)
      sum=0.
      do 11 i=1,n
        sum=sum+fvec(i)**2
11    continue
      fmin=0.5*sum
      return
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
      INTEGER n
      LOGICAL check
      REAL*8 f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX
      PARAMETER (ALF=1.e-4,TOLX=1.e-7)
      EXTERNAL func
CU    USES func
      INTEGER i
      REAL*8 a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
     *test,tmplam
      check=.false.
      sum=0.
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=sqrt(sum)
      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif
      slope=0.
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      test=0.
      do 14 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
        f=func(x)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          return
        else if(f.le.fold+ALF*alam*slope)then
          return
        else
          if(alam.eq.1.)then
            tmplam=-slope/(2.*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.)then
              tmplam=-slope/(2.*b)
            else
              disc=b*b-3.*a*slope
              tmplam=(-b+sqrt(disc))/(3.*a)
            endif
            if(tmplam.gt..5*alam)tmplam=.5*alam
          endif
        endif
        alam2=alam
        f2=f
        fold2=fold
        alam=max(tmplam,.1*alam)
      goto 1
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE qrdcmp(a,n,np,c,d,sing)
      INTEGER n,np
      REAL*8 a(np,np),c(n),d(n)
      LOGICAL sing
      INTEGER i,j,k
      REAL*8 scale,sigma,sum,tau
      sing=.false.
      scale=0.
      do 17 k=1,n-1
        do 11 i=k,n
          scale=max(scale,abs(a(i,k)))
11      continue
        if(scale.eq.0.)then
          sing=.true.
          c(k)=0.
          d(k)=0.
        else
          do 12 i=k,n
            a(i,k)=a(i,k)/scale
12        continue
          sum=0.
          do 13 i=k,n
            sum=sum+a(i,k)**2
13        continue
          sigma=sign(sqrt(sum),a(k,k))
          a(k,k)=a(k,k)+sigma
          c(k)=sigma*a(k,k)
          d(k)=-scale*sigma
          do 16 j=k+1,n
            sum=0.
            do 14 i=k,n
              sum=sum+a(i,k)*a(i,j)
14          continue
            tau=sum/c(k)
            do 15 i=k,n
              a(i,j)=a(i,j)-tau*a(i,k)
15          continue
16        continue
        endif
17    continue
      d(n)=a(n,n)
      if(d(n).eq.0.)sing=.true.
      return
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE qrupdt(r,qt,n,np,u,v)
      INTEGER n,np
      REAL*8 r(np,np),qt(np,np),u(np),v(np)
CU    USES rotate
      INTEGER i,j,k
      do 11 k=n,1,-1
        if(u(k).ne.0.)goto 1
11    continue
      k=1
1     do 12 i=k-1,1,-1
        call rotate(r,qt,n,np,i,u(i),-u(i+1))
        if(u(i).eq.0.)then
          u(i)=abs(u(i+1))
        else if(abs(u(i)).gt.abs(u(i+1)))then
          u(i)=abs(u(i))*sqrt(1.+(u(i+1)/u(i))**2)
        else
          u(i)=abs(u(i+1))*sqrt(1.+(u(i)/u(i+1))**2)
        endif
12    continue
      do 13 j=1,n
        r(1,j)=r(1,j)+u(1)*v(j)
13    continue
      do 14 i=1,k-1
        call rotate(r,qt,n,np,i,r(i,i),-r(i+1,i))
14    continue
      return
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE rsolv(a,n,np,d,b)
      INTEGER n,np
      REAL*8 a(np,np),b(n),d(n)
      INTEGER i,j
      REAL*8 sum
      b(n)=b(n)/d(n)
      do 12 i=n-1,1,-1
        sum=0.
        do 11 j=i+1,n
          sum=sum+a(i,j)*b(j)
11      continue
        b(i)=(b(i)-sum)/d(i)
12    continue
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE rotate(r,qt,n,np,i,a,b)
      INTEGER n,np,i
      REAL*8 a,b,r(np,np),qt(np,np)
      INTEGER j
      REAL*8 c,fact,s,w,y
      if(a.eq.0.)then
        c=0.
        s=sign(1.d0,b)
      else if(abs(a).gt.abs(b))then
        fact=b/a
        c=sign(1./sqrt(1.+fact**2),a)
        s=fact*c
      else
        fact=a/b
        s=sign(1./sqrt(1.+fact**2),b)
        c=fact*s
      endif
      do 11 j=i,n
        y=r(i,j)
        w=r(i+1,j)
        r(i,j)=c*y-s*w
        r(i+1,j)=s*y+c*w
11    continue
      do 12 j=1,n
        y=qt(i,j)
        w=qt(i+1,j)
        qt(i,j)=c*y-s*w
        qt(i+1,j)=s*y+c*w
12    continue
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
