      subroutine RntGravFdpnm(nmin,maxn,rln,cnm,snm,gvm,GRS,pnm,dpt1,dpt2,gr)
      !输入：nmin，maxn-最小、最大计算阶数
      !输入cnm、snm－模型重力场位系数
      !Input: nmin,maxn-the minimum and maximum calculation degree
      !Input: cnm, snm－geopotential coefficients
      !输出gvm(8)-高程异常(m)、空间异常(mGal)、扰动重力(mGal)、南向、西向垂线偏差(s)、重力梯度(E)、北向、东向水平重力梯度(E)
      !Return: gvm(8)-height anomaly (m), gravity anomaly (mGal), gravity disturbance (mGal), vertical deflection vector
      !     (ʺ, south, west), disturbing gravity gradient (E, radial), tangential gravity gradient vector (E, north, west) 
      implicit none
	integer::nmin,maxn
	real*8::cnm((maxn+2)**2),snm((maxn+2)**2),pnm((maxn+2)**2),dpt1((maxn+2)**2),dpt2((maxn+2)**2)
	real*8::cosml,sinml,tn(5),gr,pi,RAD,rln(3),BLH(3),NFD(5)
	integer n,m,kk,i,astat(6)
	real*8 rr,rlat,rlon,gm,ae,gvm(9),GRS(6),u,t,sp
!---------------------------------------------------------------------------
   	gm=GRS(1);ae=GRS(2);pi=datan(1.d0)*4.d0; RAD=pi/180.d0
	rr=rln(1);rlat=rln(2);rlon=rln(3);gvm=0.d0
      t=dsin(rlat*RAD);u=dcos(rlat*RAD)
	do n=nmin,maxn
        tn(1:5)=0.d0
	  do m=0,n
          kk=n*(n+1)/2+m
          cosml=dcos(dble(m)*rlon*RAD);sinml=dsin(dble(m)*rlon*RAD)
	    tn(1)=tn(1)+(cnm(kk)*cosml+snm(kk)*sinml)*pnm(kk+1)
          tn(2)=tn(2)+(cnm(kk)*cosml+snm(kk)*sinml)*dpt1(kk+1)
	    tn(4)=tn(4)+(cnm(kk)*cosml+snm(kk)*sinml)*dpt2(kk+1)
	    tn(3)=tn(3)+dble(m)*(cnm(kk)*sinml-snm(kk)*cosml)*pnm(kk+1)
	    tn(5)=tn(5)+dble(m**2)*(cnm(kk)*cosml+snm(kk)*sinml)*pnm(kk+1)
        enddo
        sp=dexp(dble(n)*dlog(ae/rr))
	  gvm(1)=gvm(1)+tn(1)*sp                        
	  gvm(2)=gvm(2)+tn(1)*sp*dble(n-1)              
	  gvm(3)=gvm(3)+tn(1)*sp*dble(n+1)              
	  gvm(4)=gvm(4)+tn(2)*sp                        
	  gvm(5)=gvm(5)+tn(3)*sp                        
	  gvm(6)=gvm(6)+tn(1)*sp*dble(n+1)*dble(n+2)    
	  gvm(7)=gvm(7)+tn(4)*sp                        
	  gvm(8)=gvm(8)+tn(5)*sp                        
	enddo
	gvm(9)=gm/rr*gvm(1) !m2/s2                       
	gvm(1)=(gm/rr/gr)*gvm(1)  !m
	gvm(2)=(gm/rr**2)*gvm(2)*1.d5 !mGal
	gvm(3)=(gm/rr**2)*gvm(3)*1.d5 !mGal
	sp=gm/rr**2*gvm(4);gvm(4)=sp/gr/RAD*36.d2 !s
	gvm(5)=gm/rr**2/gr*gvm(5)/u/RAD*36.d2 !s
	gvm(6)=gm/rr**3*gvm(6)*1.d9 !E
	gvm(7)=gm/rr**3*gvm(7)*1.d9-gvm(3)/rr*1.d4 !E
	gvm(8)=-gm/rr**3*gvm(8)/u**2*1.d9-gvm(3)/rr*1.d4+sp/rr*t/u*1.d9 !E
904	return
      end
