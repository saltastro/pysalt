c---------------------------------------------------------------------------------
      subroutine skyring(ax, ay, ar, rx, ry, filename, pmin, pmax,
     >     plottype, outfile)
c---------------------------------------------------------------------------------

c nsl 18.08.10: 
c This program subtracts sky rings from images.
c Removed all IRAF dependencies
c nsl 19.08.10 made compatible with pysalt package

      real sum(64),sumsq(64),ave(64),rad(64),xxx(64)
      real error(64),rms(64),fit(2),spl(64)
      real smooth(130),yetop(64),yebot(64)
c      real  indata(1024)
      real indata(7324), pmin, pmax, outdata(7324)
      integer num(64),imsize(2), imdim
      character*40 filename,outfile
      character*24 comment
      character*1 answer,bell
      character*3 plottype
      integer      blocksize, readwrite, inunit, outunit
      integer      group, nullvalj, npixels, firstpix, bitpix, ier
      integer      keyval, ipix
      logical      anynull, simple, extend
      real         nullvale, den
      integer      nelements, naxis, naxes(2), line
      real ax, ay, ar,rx,ry
      real fact,rmaxsq, rysq, dx,arsq,dy,ymax,xmax 
      real arad,rrad, rrsq, smlen, ymin, const, fac, sky
      integer ir

      print*, 'in skyring'

      bell = char(7)

c  setup initial parameters

      call pgbegin (0, plottype,1,1)
      call pgask (.false.)

c      write (*,1000)
c 1000 format ('Enter aperture xc, yc, and radius: '$)
c      read (*,*) ax,ay,ar
      arsq = ar * ar
      print*, ar, arsq, ' ar and arsq'
c      write (*,1001)
c 1001 format ('Enter ring xc, yc: '$)
c      read (*,*) rx,ry
      dx = abs(rx-ax)
      dy = abs(ry-ay)
      if (dy .gt. 0.0) then
         ymax = ar / sqrt(1.0 + (dx*dx)/(dy*dy))
      else
         ymax = ar
      endif

      xmax = sqrt(arsq - ymax*ymax)
      ymax = ymax + dy
      xmax = xmax + dx
      rmaxsq = xmax*xmax + ymax*ymax
      fact = 64.0 / rmaxsq
      print*, ax, ay, ar, rx, ry, ' ax ay ar rx ry'
      print*, dx, dy, xmax, ymax, rmaxsq, fact 

      print*, ' set up limits'
c  read sky image

 1011 format (a)

C  Get an unused Logical Unit Number to use to open the FITS file.

      ier = 0
      call getlu(inunit)
      print*, ier, 'ier after ftgiou'
      ier = 0
      readwrite = 0
      blocksize = 1
      print*, 'about to read'

      call ftopen(inunit,filename,readwrite,blocksize,ier)

      print*, inunit, ier
      if (ier .ne. 0) then
         write (0,1012) bell,filename
 1012    format (a1,'*** Error opening image: ',a,' ***')
c         go to 1
         return
      endif

      print*, 'opened image'

      call ftgknj(inunit,'NAXIS',1,2,imsize,imdim,ier)
c      call imgsiz (img, imsize, imdim, imtype, ier)

      print*, imdim, ' imdim'

      call FTGKYJ(inunit,'BITPIX', keyval,comment,ier)
      bitpix=keyval
      print*, 'keyval', keyval

      if (ier .ne. 0) then
         write (0,1013) bell
 1013    format (a1,'*** Error reading size of image ***')
c         call imclos (img, ier)
         call ftclos (inunit, ier)
c         go to 1
         return
      end if
      if (imdim .ne. 2) then
         write (0,1014) bell,imdim
 1014    format (a1,'*** Wrong image dimension: ',i2,' ***')
c         call imclos (img, ier)
         call ftclos (inunit, ier)
c         go to 1
         return
      end if
c      if (.not.((imtype .eq. 3) .or. (imtype .eq. 6))) then
      if ((keyval .ne. 16) .and. (keyval .ne. -32) 
     >     .and. (keyval .ne. 32)) then
         ier = 1
         write (0,1015) bell,keyval
 1015    format (a1,'*** Wrong image type: ',i2,' ***')
c         call imclos (img, ier)
         call ftclos (inunit, ier)
c         go to 1
         return
      end if

c      write (*,1016)
c 1016 format ("Enter min and max values to accept: "$)
c      read (*,*) pmin,pmax

c  initialize arrays

      do i=1,64
         sum(i) = 0.0
         sumsq(i) = 0.0
         num(i) = 0
      end do

c  compute mean sky and error

      npixels=imsize(1) ! no in x direction
      firstpix=1
      nullvalj=-999
      nullvale=-999.0
      group=1

      print*, imsize(1), imsize(2), ' image size'

      do line=1,imsize(2)
c         call imgl2r (img, data, line, ier)
         ier = 0
         firstpix= 1 + (line-1)*(npixels)
         call ftgpve(inunit,group,firstpix,npixels,nullvale,
     &        indata,anynull,ier)

         arad = (line - ay)**2
         rrad = (line - ry)**2
         do ipix=1,imsize(1)

            if (((ipix - ax)**2 + arad) .le. arsq) then
c               print*, ' will do', indata(ipix)
               if ((indata(ipix) .ge. pmin) .and. 
     $             (indata(ipix) .le. pmax)) then
                  rrsq = rrad + (ipix - rx)**2
                  ir = 1.0 + rrsq * fact
c                  print*, rrsq, rrad, ipix, rx, fact, ir
                  if (ir .le. 64) then
                     sum(ir) = sum(ir) + indata(ipix)
                     sumsq(ir) = sumsq(ir) + indata(ipix) * indata(ipix)
                     num(ir) = num(ir) + 1
                  end if
               end if
            end if
         end do
      end do

      do i=1,64
         xxx(i) = i
         ave(i) = sum(i) / num(i)
         rad(i) = sqrt((i - 0.5) / fact)
         rms(i) = sqrt((sumsq(i)/num(i)) - ave(i)*ave(i))
         den = num(i) - 1.0
         if (den .lt. 1.0) den = 1.0
         error(i) = rms(i) / sqrt(den)
         yetop(i) = ave(i) + error(i)
         yebot(i) = ave(i) - error(i)
c         print*, num(i), ave(i), rad(i), rms(i), error(i), 
c     >        yetop(i), yebot(i)
      end do


c  fit sky spectrum

      fit(2) = (ave(64) - ave(1)) / 63.0 ! gradient
      fit(1) = ave(1) - fit(2) ! y intercept




c  smooth spectrum

 20   write (*,1020)
 1020 format ("Enter smoothing length: "$)
      read (*,*) smlen
      ymin = yebot(1)
      ymax = yetop(1)
      do i=1,64 
         ymin = min(ymin,yebot(i))
         ymax = max(ymax,yetop(i))
         smooth(i) = ave(i) - fit(1) - fit(2) * xxx(i)
         smooth(i+64) = 0.0
      end do

      do j=1,64
         print*, fit(1), fit(2), xxx(j), ave(j), smooth(j), ' smoothed'
      enddo

c      call realft (smooth,64,1)
      call realft (smooth,128,1)
      smooth(1) = smooth(1) / 64.0
      const = (smlen/128)**2
      fac = 1.0
      do j=1,63
         k = 2 * j + 1
         if (fac .ne. 0.0) then
            fac = max(0.0,(1.0-const*j*j)/64.0)
c            print*, j,k,fac
            smooth(k) = fac * smooth(k)
            smooth(k+1) = fac * smooth(k+1)
         else
c            print*, j, k, fac, 'zero'
            smooth(k) = 0.0
            smooth(k+1) = 0.0
         end if
      end do
c      fac = amax1(0.0,(1.0-0.25*smlen*smlen)/64.0)
c      print*, fac, 'for smooth 2'
c nic changed this to make it compatible with stuff in do loop above
      fac = max(0.0,(1.0-0.25*const)/64.0)

      smooth(2) = fac * smooth(2)
c      call realft (smooth,64,-1)
      call realft (smooth,128,-1)
      do i=1,64
         smooth(i) = smooth(i) + fit(1) + fit(2) * xxx(i)
         ymin = min(ymin,smooth(i))
         ymax = max(ymax,smooth(i))
      end do


c  plot results

      call pgenv (0.0,65.0,ymin,ymax,0,1)
      call pglabel (" "," ",filename)
      call pgerry (64,xxx,yetop,yebot,1.0)
      call pgpoint (64,xxx,ave,0)
      call pgline (64,xxx,smooth)

      write (*,1021)
 1021 format ('Change smoothing length (y/n)? '$)
      read (*,1011) answer
      if (answer .eq. 'y') go to 20
      

c 21   write (*,1022)
c 1022 format ("Enter filename for subtracted image: "$)
c      read (*,1011) outfile
c      call imopnc (outfile, img, jmg, ier)

      ier = 0
      call deletefile(outfile, ier)
      call ftgiou(outunit,ier)
c     create a new empty file
      call ftinit(outunit, outfile, blocksize, ier)

      ier=0
c     copy the header from the old file to the new file
      call ftcphd(inunit, outunit, ier)

      if (ier .ne. 0) then
         write (0,1012) bell,outfile
         return
c         go to 21
      end if

      call spline (xxx,smooth,64,spl)

c     for output image:
      simple =.true.
      bitpix = -32 ! this is what the primary reduced data appear to be.
      naxis=2
      naxes(1)=imsize(1)
      naxes(2)=imsize(2)
      extend =.true.

      do line=1,imsize(2)
c         call imgl2r (img, indata, line, ier)
            firstpix= 1 + (line-1)*(npixels)
            call ftgpve(inunit,group,firstpix,npixels,nullvale,
     &           indata,anynull,ier)

         arad = (line - ay)**2
         rrad = (line - ry)**2
         do ipix=1,imsize(1)
            if (((ipix - ax)**2 + arad) .le. arsq) then
               if (indata(ipix) .gt. -32000.) then
                  rrsq = rrad + (ipix - rx)**2
                  ar = 1.0 + rrsq * fact
                  if (ar .gt. 64.0) ar = 64.0
c                  sky = smooth(ir) + (smooth(ir+1) - smooth(ir))
c     $                 * (ar - ir)
                  call splint (xxx,smooth,spl,64,ar,sky)
                  indata(ipix) = indata(ipix) - sky
               end if
            end if
         end do
         call ftppre(outunit, group, firstpix, npixels, indata,ier)
c         call impl2r (jmg,indata,line,ier)
      end do

c      call imclos (img,ierr)
c      call imclos (jmg,ierr)

      call ftclos (outunit, ier)


c      write (*,1030)
c 1030 format ('Another file (y/n)? '$)
c      read (*,1011) answer
c      if (answer .eq. 'y') go to 1

      call pgend ()
      return
c
      end subroutine skyring

