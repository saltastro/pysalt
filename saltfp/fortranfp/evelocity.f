      program evelocity

      real data(6500,9,50), fits(6500,8), wgt(9,9)
      real wave0(50), xrc(50), yrc(50), sn(50), norm(50), dn(50)
      real wave(50), flux(50), sigma(50), fit(5), error(5)
      integer isize(2), img(50), jmg(8)
      logical oldf, flag(5), anyf
      character*80 fn
      character*9 tname(8)
      character*1 bell,answ

      data tname /'cont.fit ','inty.fit ','vel.fit  ','disp.fit ',
     &            'econt.fit','einty.fit','evel.fit ','edisp.fit'/

      bell = char(7)

      flag(1) = .true.
      flag(2) = .true.
      flag(3) = .true.
      flag(4) = .true.
      flag(5) = .false.

c	Read calibration file and open images

 10   write (*,1000)
 1000 format ("Enter the image cube descriptor name: "$)
      read (*,1001) fn
 1001 format (a)
      open (10,file=fn,status="old",iostat=ierr)
      if (ierr .ne. 0) then
         write (*,1002) bell
 1002    format (a1,"Error opening calibration file!")
         go to 10
      end if
      read (10,1001) dummy
      read (10,1001) dummy
      read (10,*) rest, sg, sl, f, gain, nfile
      if (nfile .gt. 50) then
         write (*,1003) bell
 1003    format (a1,"Too many images !")
         go to 10
      end if
      fsq = f * f

      read (10,1001) dummy
      do k=1,nfile
         read (10,*) xrc(k),yrc(k),wave0(k),sn(k),norm(k),dn(k),fn
	   sn(k) = sn(k) * gain * norm(k)
	   img(k) = 19 + k
         call ftopen (img(k),fn,0,iblk,ierr)
         if (ierr .ne. 0) then
            write (*,1004) bell, fn
 1004       format (a1,"Error opening image file ",a)
            do i=1,k
               call ftclos (img(i), ierr)
            end do
            go to 10
         end if
c        call imgsiz (img(k),size,idim,ityp,ier)
c        if (size(1) .gt. 300) then
c           write (*,1005) bell
c 1005       format (a1,"Images too large!")
c           do i=1,k
c              call imclos (img(i), ierr)
c           end do
c           go to 10
c        end if
      end do

      close (10)
	call ftgipr (img(1),2,ibitpix,naxis,isize,ier)

c  calculate smoothing kernal for fwhm=3.2pix gaussian

      sum = 0.0
      do i=1,9
         do j=1,9
c            rsq = (i - 4.0)**2 + (j - 4.0)**2
c            if (rsq .le. 9.0) then
c               wgt(i,j) = exp(-0.271 * rsq)
c            else
c               wgt(i,j) = 0.0
c            end if
	      wgt(i,j) = 1.0
            sum = sum + wgt(i,j)
         end do
      end do
c      do i=1,7
c         do j=1,7
c            wgt(i,j) = wgt(i,j) / sum
c         end do
c      end do

c	get area parameters


      
C	read old fits and errors

	do i=1,isize(1)
	   fits(i,1) = -666.0
	end do
	iblk = 1

	write (*,1234)
1234	format ('Use existing fits (y/n)? '$)
	read (*,1001) answ
	if ((answ .eq. "Y") .or. (answ .eq. "y")) then
	   oldf = .true.
	else
	   oldf = .false.
	end if

      do k=1,8
	   jmg(k) = 79 + k
         if (.not. oldf) then
	      call ftinit(jmg(k),tname(k),iblk,ierr)
	      call ftphpr(jmg(k),.true.,-32,2,isize,0,1,.true.,ierr)
            if (ierr .ne. 0) then
               write (*,1008) bell
 1008          format (a1,"Error creating fit map!")
               stop
            end if
	      loc = 1
	      do j=1,isize(2)
	         call ftppre (jmg(k),0,loc,isize(1),fits(1,1),ierr)
	         loc = loc + isize(1)
	      end do
	   else
            call ftopen (jmg(k),tname(k),1,iblk,ierr)
            if (ierr .ne. 0) then
               write (*,1009) bell
 1009          format (a1,"Error opening fit map!")
               stop
            end if
         end if
      end do

c      if (oldf) then
c         write (*,1010)
c 1010    format ("Initialize the old fits (y/n)? "$)
c         read (*,1001) answ
c         if (answ .eq. "y") then
c            write (*,1011)
c 1011       format ("Are you sure (y/n)? "$)
c            read (*,1001) answ
c            if (answ .eq. "y") then
c               do k=1,8
c                  do j=1,size(2)
c                     call impl2r (jmg(k),fits(1,1),j,ierr)
c                  end do
c               end do
c            end if
c         end if
c      end if

      write (*,1006)
 1006 format ("Enter fit region center coordinates, and radius: "$)
      read (*,*) xc,yc,rmax
      rmaxsq = rmax * rmax

      iylo = max(5,ifix(yc - rmax - 1))
      iyhi = min(isize(2)-4,ifix(yc + rmax + 1))
      do 80 iy=iylo,iyhi

         yr = iy - yc
         dx = rmaxsq - yr * yr
         if (dx .lt. 0.0) go to 80
         dx = sqrt(dx)

         do k=1,nfile
	   loc = (iy-5) * isize(1) + 1
            do j=1,9
	        call ftgpve(img(k),0,loc,isize(1),0,data(1,j,k),anyf,ierr)
              loc = loc + isize(1)
            end do
         end do
         loc = (iy-1) * isize(1) + 1
         do k=1,8
            call ftgpve(jmg(k),0,loc,isize(1),0,fits(1,k),anyf,ierr)
         end do

         do 70 ix=int(xc-dx),int(xc+dx)

c            if (fits(ix,5) .gt. -500.0)  go to 70
	      if (data(ix,5,1) .lt. 0.0) go to 70

C	compute fluxes and velocities

            do k=1,nfile
               radsq = (ix - xrc(k))**2 + (iy - yrc(k))**2
               wave(k) = wave0(k) / sqrt(1.0 + (radsq / fsq))
               sum = 0.0
               ssum = 0.0
               do ii=1,9
                  do jj = 1,9
                     ipix = ix + jj - 5
                     flx = data(ipix,ii,k) * norm(k) * gain
                     ssig = abs(flx) + sn(k)**2 + (flx*dn(k))**2
                     sum = sum + flx * wgt(ii,jj)
                     ssum = ssum + wgt(ii,jj)**2 * ssig
                  end do
               end do
               flux(k) = sum
               sigma(k) = sqrt(ssum)
            end do
            call evinit (wave,flux,nfile,fit)
            fit(4) = sg
            fit(5) = sl
            call evfit(wave,flux,sigma,nfile,fit,flag,error,chi)
            if (chi .gt. 0.0) then
               fits(ix,1) = fit(1)
               fits(ix,2) = fit(2)
               fits(ix,3) = (fit(3) - rest) * 2.997925e+05 / rest
               fits(ix,4) = abs(fit(4)) * 2.997925e+05 / fit(3)

               fits(ix,5) = error(1)
               fits(ix,6) = error(2)
               fits(ix,7) = error(3) * 2.997925e+05 / rest
               fits(ix,8) = error(4) * 2.997925e+05 / fit(3)
            end if
 70      continue
         loc = (iy-1) * isize(1) + 1
         do k=1,8
            call ftppre (jmg(k),0,loc,isize(1),fits(1,k),ierr)
         end do

 80   continue


      do i=1,nfile
         call ftclos (img(i),ierr)
      end do

      do i=1,8
         call ftclos (jmg(i),ierr)
      end do

      stop
      end
