      subroutine eprofile(plottype, icdname)

c 16.08.10. nsl updated to be compatible with pyraf pipeline

c  RSS version  


      real dataarray(6500)
      real xrc(50), yrc(50), wave0(50), sky(50), norm(50), dn(50)
      real wave(50), flux(50), sig(50), ftop(50), fbot(50)
	integer img(50), isize(2)

      real xpl(100),ypl(100)
      real error(5), dy(5), fit(5), chisq, wgt(9,9)
      
      character*80 fname, icdname
      character*4 plottype
      logical flag(5), anyf

c      call pgbegin (0,"?",1,1)
      call pgbegin (0,plottype,1,1)
      call pgask (.false.)

c calculate a fixed smoothing kernal (9x9 for now)

      sum = 0.0
      do i=1,9
         do j=1,9
c            rsq = (i - 4.0)**2 + (j - 4.0)**2
c            if (rsq .le. 9.0) then
c               wgt(i,j) = exp(-0.271 * rsq)
c            else
c               wgt(i,j) = 0.0
c            end if
	      wgt (i,j) = 1.0 / 81.0
            sum = sum + wgt(i,j)
         end do
      end do
c      do i=1,7
c         do j=1,7
c            wgt(i,j) = wgt(i,j) / sum
c         end do
c      end do

c  read calibration data on images
 
      write(*,1000)
 1000 format ("Line Profile Fitting and Display")
c 10   write (*,1001)
c 1001 format ("Enter the name of the image cube descriptor file: "$)
c      read (*,1002) fname
 1002 format (a)
c      open (10,file=fname,status="old",err=10)
      ier = 0
      open (10,file=icdname,status="old",iostat=ier)
      if (ier .ne. 0) then
         write (*, '("Error opening image descriptor file")')
         stop
      end if      

      read (10,1002) dummy
      read (10,1002) dummy
      read (10,*) rwave,sg,sl,f,gain,nimage
      fsq = f * f

      read (10,1002) dummy
      do i=1,nimage
         read (10,*) xrc(i),yrc(i),wave0(i),sky(i),norm(i),dn(i),fname
         sky(i) = sky(i) * gain * norm(i)
c	   call getlu (img(i))
	   img(i) = 19 + i
         call ftopen (img(i),fname,0,iblk,ier)
         if (ier .ne. 0) then
            write (*,1003) fname
 1003       format ("Error opening file ",a)
            if (i .gt. 1) then
               do j=1,i
                  call ftclos (img(j),ier)
               end do
            end if
c            go to 10
c nsl commented this out.. dont see why have to keep specifying the descriptor file
c as this has a list of filenames in it itself. Need to check usage with Ted.

         end if
      end do

      close (10)
      ier = 0
      call ftgipr (img(1),2,ibitpix,naxis,isize,ier)

c	get coordinates of next pixel to display

   20 write (*,1010)
 1010 format (" Enter pixel coordinates (0,0 to stop): "$)
      read (*,*) ixp,iyp

      if (ixp .eq. 0) then
         do i=1,nimage
            call ftclos (img(i),ier)
         end do
         call pgend ()
         return
c         stop, nsl added return
      end if

      ixp = max (5,ixp)
      ixp = min (isize(1)-4,ixp)
      iyp = max (5,iyp)
      iyp = min (isize(2)-4,iyp)
      
c      write (*,1011)
c 1011 format ("Enter smoothing size and fit order: "$)
c      read (*,*) ibox,iord

c	read data for that pixel
        
      print*, isize(1), isize(2), 'array size'

      do i=1,nimage
         radsq = (ixp - xrc(i))**2 + (iyp - yrc(i))**2
         wave(i) = wave0(i) / sqrt(1.0 + (radsq / fsq))
         sum = 0.0
         ssum = 0.0
         do ii=1,9
            line = iyp - 5 + ii
            loc = (line - 1) * isize(1) + 1
            call ftgpve (img(i),0,loc,isize(1),0.0,dataarray,anyf,ier)
            
            do jj = 1,9
               ipix = ixp - 5 + jj
               flx = (dataarray(ipix) * norm(i) * gain) - sky(i)
               ssig = abs(flx) + abs(sky(i)) + (flx*dn(i))**2
               sum = sum + flx * wgt(ii,jj)
               ssum = ssum + wgt(ii,jj)**2 * ssig
            end do
         end do
         flux(i) = sum
         sig(i) = sqrt(ssum)
         ftop(i) = flux(i) + sig(i)
         fbot(i) = flux(i) - sig(i)
      end do
      
c     
c     fit profile
c
      
      do i = 1, 4
         flag( i ) = .true.
      end do
      flag(5) = .false.
      
      call evinit (wave, flux, nimage, fit)
      fit (4) = sg
      fit(5) = sl
      call evfit (wave, flux, sig, nimage, fit, flag, error, chisq)
      write(*, '(a)') 'Fit coeff and error' >' (coefficients 1,2,3,4)' 
      write (*,1021) (fit(i),error(i),i=1,4)
 1021 format (4(5x,f7.1,1x,f7.1))
      
c     plot fit
      
      wmin = wave(1)
      wmax = wave(1)
      fmin = fbot(1)
      fmax = ftop(1)
      do i=2,nimage
         wmin = min (wave(i),wmin)
         wmax = max (wave(i),wmax)
         fmin = min (fbot(i),fmin)
         fmax = max (ftop(i),fmax)
      end do
      
      wav = wmin
      winc = (wmax - wmin) / 99.0
      do i=1,100
         xpl(i) = wav
         call evoigt (wav, fit, 5, ypl(i), dy)
         fmin = min (ypl(i),fmin)
         fmax = max (ypl(i),fmax)
         wav = wav + winc
      end do
      
      wmin = wmin - 1.0
      wmax = wmax + 1.0
      temp = 0.05 * (fmax - fmin)
      fmin = fmin - temp
      fmax = fmax + temp
      
      call pgenv (wmin,wmax,fmin,fmax,0,1)
      call pglabel ("wavelength","intensity","")
      call pgpoint (nimage,wave,flux,0)
      call pgerry (nimage,wave,ftop,fbot,1.0)
      call pgline (100,xpl,ypl)
      
      go to 20      
        
      end subroutine eprofile
      
