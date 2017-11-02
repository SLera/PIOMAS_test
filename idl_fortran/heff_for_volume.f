C Calculate PIOMAS monthly mean sea ice volume for the Arctic

      parameter(nx1=360,ny1=120,nx=nx1,ny=ny1,imt=nx1,jmt=ny1,nt=1)
      parameter(m1=12,idelt=4,nregion=4,cutoff=0.15)
      parameter(nyear1=1979,nyear2=2014,nmon1=1,nmon2=12)

      dimension heff(imt,jmt),h(m1),xn(nx1,ny1),yn(nx1,ny1)
      dimension heff1(imt,jmt), therm(imt,jmt), tair(imt,jmt)
      dimension heffm(imt,jmt,12),work(imt,jmt)
      dimension tair1(imt,jmt),tairm(imt,jmt,12)
      dimension clon(imt,jmt),clat(imt,jmt),kmt(imt,jmt)
      dimension idomain(imt,jmt)
      dimension mask(nx1,ny1,nregion),count(nregion)
     &,count_overice(nregion), hmon(nregion), hyr(nregion)
     &,vmon(nregion), vyr(nregion)
     &,month(12),month_end(13)

      dimension ulat(imt,jmt),ulon(imt,jmt),HTN(imt,jmt),HTE(imt,jmt)
     &,HUS(imt,jmt),HUW(imt,jmt),angle(imt,jmt),dxt(imt,jmt)
     &,dyt(imt,jmt)

      character *80 fopen(5), f1,f2,f3
      character *4 cyear(1900:2100),cyear1(1900:2100)
      character *12 char
      integer SLEN

      data month/31,28,31,30,31,30,31,31,30,31,30,31/
      data month_end/0,31,59,90,120,151,181,212,243,273,304,334,365/

c heff is mean sea ice thickness in a greid cell
      f1='heff.H'

      open(3,file='heff_mon_Data')
      open(4,file='heff_yr_Data')

c clon/clat are the lon/lat at the center of heff grid cell 
      open(20,file='/home/zhang/POP_1.2/360_120/grid.dat.rot')
      read(20,'(10f8.2)') ((clon(i,j),i=1,nx1),j=1,ny1)
      read(20,'(10f8.2)') ((clat(i,j),i=1,nx1),j=1,ny1)
      close(20)

      open(24,file='/home/zhang/POP_1.2/360_120/grid.dat.pop')
        read(24,'(10f8.2)') ((ulat(i,j),i=1,nx),j=1,ny)
        read(24,'(10f8.2)') ((ulon(i,j),i=1,nx),j=1,ny)
        read(24,'(10f8.2)') ((HTN  (i,j),i=1,nx),j=1,ny)
        read(24,'(10f8.2)') ((HTE  (i,j),i=1,nx),j=1,ny)
        read(24,'(10f8.2)') ((HUS  (i,j),i=1,nx),j=1,ny)
        read(24,'(10f8.2)') ((HUW  (i,j),i=1,nx),j=1,ny)
        read(24,'(10f8.2)') ((angle(i,j),i=1,nx),j=1,ny)
      close(24)

C DXT and DYT are in km
      c0=0.0
      c1=1.0
      p5=0.5
      WORK=eoshift(HTN,dim=2,shift=-1)
      DXT = p5*(HTN + WORK)
      WORK=eoshift(HTE,dim=1,shift=-1)

      DYT = p5*(HTE + WORK)
      where (DXT == c0) DXT=c1
      where (DYT == c0) DYT=c1

c kmt is grid mask (kmt=0, land), (kmt>0, ocean)
      open(20,file='io.dat_360_120.output',form='formatted')
      read(20,'(360i2)') kmt
      close(20)

      do i=1,nx1
      do j=1,ny1
      do k=1,nregion
      mask(i,j,k)=0
      end do
      end do
      end do

c mask(3) is arctic basin + barents sea
      do i=99+1,265+1
      do j=33+1,ny1-3+1
      if(kmt(i,j).gt.0) mask(i,j,3)=1
      end do
      end do
      do i=96+1,98+1
      do j=ny1-23,2,-1
      if(kmt(i,j).gt.0) mask(i,j,3)=1
      end do
      end do

      k=0
      jk=0
      do j=ny1-29,ny1-29-25,-1
        k=k+1
        if (mod(k,2) .eq. 1) jk=jk+1
        do i=92-jk+1, 92-jk+50+1
          if(kmt(i,j).gt.0) mask(i,j,3)=1
        end do
      end do

      do i=129-50+1,130+1
      do j=1+1,64+1
        if(kmt(i,j).gt.0) mask(i,j,3)=1
      end do
      end do

c mask(1) is GIN Sea + South of GIN Sea
      do i=1,nx1/2
      do j=2,ny1-1
      if(kmt(i,j).gt.0.and.mask(i,j,3).eq.0) then
      mask(i,j,1)=1
      end if
      end do
      end do

c mask(2) is Bering Sea
      do i=nx1/2-30,nx1-80
      do j=2,33+1
      if(kmt(i,j).gt.0.and.mask(i,j,3).eq.0) then
      mask(i,j,2)=1
      end if
      end do
      end do

c mask(4) for whole PIOMAS domain
      do i=1,nx1
      do j=2,ny1-1
      if(kmt(i,j).gt.0) then
      mask(i,j,4)=1
      end if
      end do
      end do


      do 999 iyear=nyear1,nyear2

      write(*,*) iyear

      write(unit=cyear(iyear),fmt='(i4)') iyear
      i=slen(f1)
      open(1,file=f1(1:i)//cyear(iyear)
     &,access='direct',form='unformatted',recl=nx1*ny1*4
     &,status='unknown')

      do i=1,nregion
      hyr(i)=0.0
      vyr(i)=0.0
      end do

c Note: some machines need a record number to read data as below
c      irec = 0
c      do imon=nmon1,nmon2
c        irec = irec+1
c        read(5,rec=irec)((heff(i,j),i=1,nx1),j=1,ny1)
c      end do

      do imon=nmon1,nmon2

      read(1)((heff(i,j),i=1,nx1),j=1,ny1)

      heffm(:,:,imon)=heff

      do i=1,nregion
      hmon(i)=0.0
      vmon(i)=0.0
      count_overice(i)=0.0
      end do

c hmon is area mean ice thickness; vmon is ice volume
      do i=1,imt
      do j=2,jmt-1
      do ia=1,nregion
      if(mask(i,j,ia).eq.1) then
      vmon(ia)=vmon(ia)+heffm(i,j,imon)*dxt(i,j)*dyt(i,j)
     &*1000.0*1000.0
      end if
      if(mask(i,j,ia).eq.1) then
      count_overice(ia)=count_overice(ia)+dxt(i,j)*dyt(i,j)
      hmon(ia)=hmon(ia)+heffm(i,j,imon)*dxt(i,j)*dyt(i,j)
      end if
      end do
      end do
      end do

      do ia=1,nregion
      count_overice(ia)=max(count_overice(ia),1.0)
      end do

      do ia=1,nregion
      hmon(ia)=hmon(ia)/count_overice(ia)
      vmon(ia)=vmon(ia)*1.0E-12  ! 10**12 m**3
      hyr(ia)=hyr(ia)+hmon(ia)
      vyr(ia)=vyr(ia)+vmon(ia)
      end do

      write(3,'(2i5,4f9.3,4f9.1)')iyear,imon
     &,hmon(1),hmon(2),hmon(3),hmon(4)
     &,vmon(1),vmon(2),vmon(3),vmon(4)

      end do

      do ia=1,nregion
      hyr(ia)=hyr(ia)/(nmon2-nmon1+1.0)
      vyr(ia)=vyr(ia)/(nmon2-nmon1+1.0)
      end do

      write(4,'(i5,4f9.3,4f9.1)')iyear,hyr(1),hyr(2),hyr(3),hyr(4)
     &,vyr(1),vyr(2),vyr(3),vyr(4)
      write(*,'(i5,4f9.3,4f9.1)')iyear,hyr(1),hyr(2),hyr(3),hyr(4)
     &,vyr(1),vyr(2),vyr(3),vyr(4)

      close(1)
      close(2)
      close(12)
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
