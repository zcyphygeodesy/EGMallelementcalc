      subroutine Ellipsoidpara(GRS)
      implicit none
      integer::i,n
	real*8::GRS(6),gm,a,j2,w,f,pwr,qsv,esq,fact,ep2,ds,twoqp
	real*8::b,e2,m,q0,E,U0,ep,fx,U,DU,dE
!---------------------------------------------------------------------
	gm=GRS(1);a=GRS(2);j2=GRS(3);w=GRS(4);f=GRS(5);U0=GRS(6)
	pwr=(w*a)**2*a/gm
      if(j2>0)then
        qsv=3.d0*j2+pwr;esq=qsv
        do i=1,10
          ep2=esq/(1.d0-esq)
          fact=1.d0/(1.d0-esq)**1.5d0
          ds=1.d0;twoqp=0.d0
!c        this loop computes 2*q/e**3
          do n=1,12
            twoqp=twoqp+ds*4.d0*n/((2.d0*n+1.d0)*(2.d0*n+3.d0))*fact
            fact=fact*ep2;ds=-ds
          enddo
          esq=qsv+pwr*(4.d0/15.d0/twoqp-1.d0)
          f=1.d0-dsqrt(1.d0-esq)
        enddo
	  b=(1.d0-f)*a;m = w*w*b*a*a/gm ;E=dsqrt(a**2-b**2)
        U0=gm/E*datan(E/b)+(w*a)**2/3.d0
        GRS(6)=U0;GRS(5)=f
        return
      endif
      if(f>0)then
	  b=(1.d0-f)*a;m = w*w*b*a*a/gm
        ep2=(a/b)**2-1.d0;e2=1.d0-(b/a)**2
        !calculate 2q0=twoqp
        twoqp=0.d0
        do n=1,12
           twoqp=twoqp+(-1.d0)**(n+1)*4.d0*n/((2.d0*n+1.d0)*(2.d0*n+3.d0))*ep2**n*dsqrt(ep2)
        enddo
        q0=twoqp/2.d0
        j2=e2/3.d0*(1.d0-2.d0/15.d0*m*dsqrt(ep2)/q0)
        E=dsqrt(a**2-b**2)
        U0=GM/E*datan(dsqrt(ep2))+(w*a)**2/3.d0
	  GRS(3)=j2;GRS(6)=U0
        return
      endif
      if(U0>0)then
	  !U0=gm/a*sqrt(1+ep2)/ep*arctg(ep)+w**2*a**2/3
        !calculate ep
        f=1.d0/298.2564115287d0;b=a*(1.d0-f)
        E=dsqrt(a**2-b**2);dE=1.d0
        do while(dabs(dE)>1.d-5)
          U=GM/E*datan(E/b)+(w*a)**2/3.d0
          dU=U0-U
          fx=GM/E*(b/a**2-datan(E/b))
          dE=dU/fx
          E=E+dE;b=dsqrt(a**2-E**2)
        enddo
        f=(a-b)/a;e2=1.d0-(b/a)**2
	  m = w*w*b*a*a/gm;ep2=(a/b)**2-1.d0
        !calculate 2q0=twoqp
        twoqp=0.d0
        do n=1,12
           twoqp=twoqp+(-1.d0)**(n+1)*4.d0*n/((2.d0*n+1.d0)*(2.d0*n+3.d0))*ep2**n*dsqrt(ep2)
        enddo
        q0=twoqp/2.d0
        j2=e2/3.d0*(1.d0-2.d0/15.d0*m*dsqrt(ep2)/q0)
	  GRS(3)=j2;GRS(5)=f
        return
      endif
      end
