c
c advect.H: monthly change in ice thickness due to convergence (meter of ice/s)
c iceprod.H: monthly net ice production (m of ice/s)
c oflux.H: monthly ocean heat flux (m of ice/s) **do not use**
c tice0.H: temperature at the surface of ice, snow, or open ocean (K)
c gice.H: monthly 12-category ice thickness distribution (area fraction)
c heff.H:    monthly ice thickness (m)
c area.H:    month ice concentration (fraction)
c icevel.H:  monthly ice velocity (m/s)
c ssh.H:     monthly SSH (cm)
c uo1_10.H:     upper 10 levels of monthly ocean velocity (cm/s)
c otemp1_10.H:  upper 10 levels of monthly ocean temperature (C)
c osali1_10.H:  upper 10 levels of monthly ocean salinity (multiplied by 1000 to get psu)

      parameter(nx1=360,ny1=120,nx=nx1,ny=ny1,imt=nx1,jmt=ny1,km=30)
      parameter(m1=12)
      parameter(nyear1=2009,nyear2=2009,nday=365,nmon=12)

      dimension advect(imt,jmt),prod(imt,jmt),oflux(imt,jmt)
      dimension tice0(imt,jmt)
      dimension heff(imt,jmt),area(imt,jmt),gice(imt,jmt,m1)
      dimension uice(imt,jmt),vice(imt,jmt),ssh(imt,jmt)
      dimension uo(imt,jmt,km),vo(imt,jmt,km)
      dimension uo1(imt,jmt,km),vo1(imt,jmt,km)
      dimension to(imt,jmt,km),so(imt,jmt,km)
      dimension clon(imt,jmt),clat(imt,jmt),kmt(imt,jmt)
      dimension kmu(imt,jmt),iwork(imt,jmt)
      dimension ulat(imt,jmt),ulon(imt,jmt),HTN(imt,jmt),HTE(imt,jmt)
     &,HUS(imt,jmt),HUW(imt,jmt),angle(imt,jmt),dxt(imt,jmt)
     &,dyt(imt,jmt)
      dimension dz(km),ZDZ (KM),ZDZZ(KM+1)
      character *80 fopen(5),f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12
      character *4 cyear(1900:2100),cyear1(1900:2100)
      character *12 char
      integer SLEN

c advect.H: monthly change in ice thickness due to convergence (meter of ice/s)
c iceprod.H: monthly net ice production (m of ice/s)
c oflux.H: monthly ocean heat flux (m of ice/s)
c tice0.H: temperature at the surface of ice, snow, or open ocean (K)
c gice.H: monthly 12-category ice thickness distribution (area fraction)
c heff.H:    monthly ice thickness (m)
c area.H:    month ice concentration (fraction)
c icevel.H:  monthly ice velocity (m/s)
c ssh.H:     monthly SSH (cm)
c uo1_10.H:     first 10 levels of monthly ocean velocity (cm/s)
c otemp1_10.H:  first 10 levels of monthly ocean temperature (C)
c osali1_10.H:  first 10 levels of monthly ocean salinity (multiplied by 1000 to get psu)

      f1='advect.H'
      f2='iceprod.H'
      f3='oflux.H'
      f4='tice0.H'
      f5='gice.H'
      f6='heff.H'
      f7='area.H'
      f8='icevel.H'
      f9='ssh.H'
      f10='uo1_10.H'
      f11='otemp1_10.H'
      f12='osali1_10.H'

c read lon and lat for scalar fields (like sea ice thickness and concentration)
      open(20,file='grid.dat.rot')
      read(20,'(10f8.2)') ((clon(i,j),i=1,nx1),j=1,ny1)
      read(20,'(10f8.2)') ((clat(i,j),i=1,nx1),j=1,ny1)
      close(20)

c read lon and lat for vector fields (like sea ice and ocean veclocities)
      open(24,file='grid.dat.pop')
        read(24,'(10f8.2)') ((ulat(i,j),i=1,nx),j=1,ny)
        read(24,'(10f8.2)') ((ulon(i,j),i=1,nx),j=1,ny)
c HTN, HTE are lengths of the northern and eastern sides of a scaler grid cell in km, HTN*HTE is the area of a scaler grid cell in km**2
c HUS, HUW are lengths of the southern and western sides of a vector grid cell in km
        read(24,'(10f8.2)') ((HTN  (i,j),i=1,nx),j=1,ny) 
        read(24,'(10f8.2)') ((HTE  (i,j),i=1,nx),j=1,ny)
        read(24,'(10f8.2)') ((HUS  (i,j),i=1,nx),j=1,ny)
        read(24,'(10f8.2)') ((HUW  (i,j),i=1,nx),j=1,ny)
c angle is the angle between latitude line and  grid cell x-coordinate line, needed for vector rotation for plotting the vectors in spherical lat-lon coordinate system
c do rotation like this for vector uo and vo
c        uo1(:,:,k)=uo(:,:,k)*cosd(-angle)+vo(:,:,k)*sind(-angle)
c        vo1(:,:,k)=vo(:,:,k)*cosd(-angle)-uo(:,:,k)*sind(-angle)
        read(24,'(10f8.2)') ((angle(i,j),i=1,nx),j=1,ny)
      close(24)

c read model scaler grid mask; ocean levels > 0, land = 0
      open(20,file='io.dat_360_120.output')
      read(20,'(360i2)') kmt
      close(20)
c get model vector grid mask; ocean levels > 0, land = 0
      IWORK=eoshift(KMT,dim=1,shift=1)
      KMU = min(KMT,IWORK)
      IWORK=eoshift(KMU,dim=2,shift=1)
      KMU = min(KMU,IWORK)  ! KMU = minimum of surrounding KMTs

c read ocean level thickness dz(k) in cm
c zdz(kmt(i,j)) is the ocean depth at scaler grid cell (i,j)
c zdz(kmu(i,j)) is the ocean depth at vector cell (i,j)
      open(20,file='dz.dta30')
      do k=1,km
        read(20,*)dz(k)
        dz(k)=dz(k)*0.01
      end do
      ZDZ(1)=DZ(1)
      do k=2,km
        ZDZ(K)=ZDZ(K-1)+DZ(K)
      end do
      do k=1,km
        zdzz(k)=zdz(k)-0.5*dz(k)
      end do
      close(20)

      do 999 iyear=nyear1,nyear2

      write(unit=cyear(iyear),fmt='(i4)') iyear

      i=slen(f1)
      open(1,file=f1(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')

      i=slen(f2)
      open(2,file=f2(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')

      i=slen(f3)
      open(3,file=f3(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')

      i=slen(f4)
      open(4,file=f4(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f5)
      open(5,file=f5(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f6)
      open(16,file=f6(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f7)
      open(7,file=f7(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f8)
      open(8,file=f8(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f9)
      open(9,file=f9(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f10)
      open(10,file=f10(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f11)
      open(11,file=f11(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
      i=slen(f12)
      open(12,file=f12(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='old')
 
c Note: some machines need a record number to read data like the following
c      irec = 0
c      do imon=1,nmon
c        do k=1,m1 
c          irec = irec+1
c          read(5,rec=irec)((gice(i,j,k),i=1,nx1),j=1,ny1)
c        end do
c      end do

c The following is for machines that do not need a record number to read data
      do imon=1,nmon
        read(1)((advect(i,j),i=1,nx1),j=1,ny1)
        read(2)((prod(i,j),i=1,nx1),j=1,ny1)
        read(3)((oflux(i,j),i=1,nx1),j=1,ny1)
        read(4)((tice0(i,j),i=1,nx1),j=1,ny1)
        do k=1,m1 
          read(5)((gice(i,j,k),i=1,nx1),j=1,ny1)
        end do
        read(16)((heff(i,j),i=1,nx1),j=1,ny1)
        read(7)((area(i,j),i=1,nx1),j=1,ny1)
        read(8)((uice(i,j),i=1,nx1),j=1,ny1)
        read(8)((vice(i,j),i=1,nx1),j=1,ny1)
        read(9)((ssh(i,j),i=1,nx1),j=1,ny1)
c        do k=1,km
        do k=1,10
          read(10)((uo(i,j,k),i=1,nx1),j=1,ny1)
        end do
c        do k=1,km
        do k=1,10
          read(10)((vo(i,j,k),i=1,nx1),j=1,ny1)
        end do
c        do k=1,km
        do k=1,10
          read(11)((to(i,j,k),i=1,nx1),j=1,ny1)
        end do
c        do k=1,km
        do k=1,10
          read(12)((so(i,j,k),i=1,nx1),j=1,ny1)
        end do
      end do

      i=nx/2
      j=ny/2
      write(*,'(2i6, 2f12.3, 2f8.2)') iyear, iday
     &, clon(i,j), clat(i,j), uo(i,j,1), heff(i,j)

      close(2)
      close(3)
      close(4)
      close(5)
      close(16)
999   continue
      
      stop  
      end

      INTEGER FUNCTION slen (string)
C ---
C --- this function computes the length of a character string less
C --- trailing blanks
C --- slen > 0, length of string less trailing blanks
C ---      = 0, character string is blank
C ---
      CHARACTER*(*) string
      CHARACTER*1 cblank
      INTEGER i
      DATA cblank/' '/
C ---
      DO 50 i = LEN(string), 1, -1
         IF (string(i:i) .NE. ' ')  GO TO 100
50    CONTINUE
      i = 0
100   CONTINUE
      slen = i
      RETURN
      END
