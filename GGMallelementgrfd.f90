!  GGMallelementgrfd.f90 
!
!  FUNCTIONS:
!  GGMallelementgrfd - Entry point of console application.
!
!****************************************************************************

      program GGMallelementgrfd
      implicit none
	character*800::line
 	real*8::BLH(3),rec(800),gvm(9),tmp(2),NFD(5),gr,rlat0,ep
      integer::n,m,i,j,nmin,maxn,sn,len,nm,astat(5),kk
	real*8,allocatable::cnm(:),snm(:)
	real*8,allocatable::pnm(:),dpt1(:),dpt2(:)
	integer::status=0
	real*8::GRS(6),djn(80),rgm,ra,t,rln(3),pi,RAD
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378136.3d0; GRS(3)=1.0826359d-3;GRS(4) = 7.292115d-5
      call ELLIPSOIDPARA(GRS)!Calculate GRS(5)=1/f
	pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      !nmin，maxn-最小、最大计算阶数
      !nmin,maxn-the minimum and maximum calculation degree
      nmin=2;maxn=360 
 	allocate(cnm((maxn+2)**2), stat=astat(1))
 	allocate(snm((maxn+2)**2), stat=astat(2))
	if (sum(astat(1:2)) /= 0) goto 908
      cnm=0.d0;snm=0.d0
      !打开地球重力位系数模型
      !Open global geopotential coefficient model.
      open(unit=8,file="EGM2008.gfc",status="old",iostat=status)
      if(status/=0)goto 904
      read(8,'(a)') line
      call PickRecord(line,len,rec,sn)
      rgm=rec(1)*1.d14/GRS(1);ra=rec(2)/GRS(2)
      do while(.not.eof(8))  
        read(8,'(a)') line
        if(len_trim(line)<2)goto 701
        call PickRecord(line,len,rec,sn)
        if(sn<4)goto 701
        n=nint(rec(1));m=nint(rec(2));tmp(1:2)=rec(3:4)
        if(n>maxn)goto 701
        if(m>n)goto 701
        if(n<maxn+1.and.n>0)then
           cnm(n*(n+1)/2+m)=tmp(1);snm(n*(n+1)/2+m)=tmp(2)
        endif
701     continue
      enddo
901   close(8)
      if(n<maxn)maxn=n  
      call normdjn(GRS,djn); call cstonorm(80,djn,GRS,cnm,snm)
 	allocate(pnm((maxn+2)**2),dpt1((maxn+2)**2),dpt2((maxn+2)**2))
      open(unit=8,file="calcpnt.txt",status="old",iostat=status)
      if(status/=0)goto 908
      open(unit=10,file="reslt.txt",status="replace")
      read(8,'(a)') line    !读取头文件
      write(10,'(a)')trim(line)
      rlat0=999.d0;kk=0
      do while(.not.eof(8))  
         read(8,'(a)') line
         call PickRecord(line,len,rec,sn)
         if(sn<4)goto 906   
         BLH(2)=rec(2);BLH(1)=rec(3);BLH(3)=rec(4)
         call BLH_RLAT(GRS,BLH,rln);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
         !一次性计算勒让德函数，提速
         !One-time calculation of Legendre functions for speed-up
         if(dabs(rln(2)-rlat0)>3.d-5)then!!!!!!!!!!!!!!!!!!!!!!!!
           t=dsin(rln(2)*RAD); call BelPnmdt(pnm,dpt1,dpt2,maxn,t);rlat0=rln(2)
         endif!!!!!!!!!!!!!!!!!!!!!
         call RntGravFdpnm(nmin,maxn,rln,cnm,snm,gvm,GRS,pnm,dpt1,dpt2,gr)
         gvm(9)=gvm(6)+gvm(7)+gvm(8);kk=kk+1
         write(10,'(a,10F12.4)')trim(line),(gvm(i),i=1,9)
         if(kk/200*200==kk)write(*, '(a,i9)'), '    Calculated number: ',kk
906	   continue
      enddo
903   close(8)
      close(10)
      deallocate(pnm,dpt1,dpt2)
904   deallocate(cnm,snm)
      write (*,*)'  Complete the computation! The results are saved in the file reslt.txt.'
      pause
908   continue
      end

