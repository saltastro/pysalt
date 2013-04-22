      
      real    data(1000,512),biw(512),sb(512)
      real    yetop(512),yebot(512),smooth(1024),fit(2),wave(512)
      real    dat(1000),xxx(512),spl(512)
      integer iaxlen(7),num(512)
      character*80 filename,outfile
      character*1  answer,bell

      data a,b,c,d,f /6739.832,119.9739,-0.57941,0.08697,5631.551/

      bell = char(7)
      call pgbegin (20,'?',1,1)
      call pgask (.false.)
 
C	get centers
 
      write (*,1000)
 1000 format ("Enter aperture xc, yc, and rmax: "$)
      read (*,*) axc,ayc,armax
      armxsq = armax * armax
      write (*,1001)
 1001 format ("Enter ring xc and yc: "$)
      read (*,*) rxc,ryc
      rmax = armax + sqrt((axc-rxc)**2 + (ayc-ryc)**2)
      rrmxsq = rmax * rmax
      write (*,1100)
 1100 format ("Enter calibration a, b, c, d, f: "$)
      read (*,*) a,b,c,d,f
      write (*,1101)
 1101 format ("Enter number of radial bins: "$)
      read (*,*) nbin

c	file loop

 10   write (*,1002)
 1002 format ("Enter image filename: "$)
      read (*,1003) filename
 1003 format (a)
      call imopen (filename, 1, im, ier)
      if (ier .ne. 0) then
         write (*,1004) bell
 1004    format (a1,"Problem opening file")
         call imclos (im, ier)
         go to 10
      end if
      call imgsiz (im, iaxlen, naxis, itype, ier)
      if (ier .ne. 0) then
         write (*,1005) bell
 1005    format (a1,"Problem getting image size")
         call imclos (im, ier)
         go to 10
      end if
      if (naxis .ne. 2) then
         write (*,1006) bell
 1006    format (a1,"Wrong image dimension")
         call imclos (im, ier)
         go to 10
      end if
      nx = iaxlen(1)
      ny = iaxlen(2)
      if ((itype .ne. 3) .and. (itype .ne. 6)) then
         write (*,1007) bell
 1007    format (a1,"Wrong image type")
         call imclos (im, ier)
         go to 10
      end if

      do i=1,nbin
         num(i) = 0
      end do

      do iy=1,ny
         call imgl2r (im, dat, iy, ier)
         if (ier .ne. 0) then
            write (*,1008) bell
 1008       format (a1,"Error reading image")
            call imclos (im, ier)
            go to 10
         end if
         ars = (ayc - iy)**2
         if (ars .le. armxsq) then
            rrs = (ryc - iy)**2
            do ix=1,nx
               arsq = ars + (axc - ix)**2
               if (arsq .le. armxsq) then
                  rrsq = rrs + (rxc - ix)**2
                  ir = nbin * rrsq / rrmxsq + 0.5
                  if (ir .lt. 1)  ir = 1
                  if (ir .gt. nbin) ir = nbin
                  if (num(ir) .lt. 1000) then
                     num(ir) = num(ir) + 1
                     data(num(ir),ir) = dat(ix)
                  end if
               end if
            end do
         end if
      end do

      write (*,1012)
 1012 format ("Enter z for file: "$)
      read (*,*) z
      wave0 = a + b*z + c*z*z + d*z*z*z
      loc = 0
      do ir=1,nbin
         if (num(ir) .gt. 10) then
            loc = loc + 1
            xxx(loc) = ir
            rad = sqrt(rrmxsq * ir / nbin)
            wave(loc) = wave0 * cos(atan2(rad,f))
            call biwgt (data(1,ir),num(ir),biw(loc),sb(loc))
            sb(loc) = sb(loc) / sqrt(num(ir) - 1.0)
            num(loc) = num(ir)
            yetop(loc) = biw(loc) + sb(loc)
            yebot(loc) = biw(loc) - sb(loc)
         end if
      end do

c  remove linear trend for FFT

      fit(2) = (biw(loc) - biw(1)) / (loc - 1)
      fit(1) = biw(1) - fit(2)


c  smooth spectrum

 220  write (*,1220)
 1220 format ("Enter smoothing length: "$)
      read (*,*) smlen
      ymin = yebot(1)
      ymax = yetop(1)
      do i=1,1024
         smooth(i) = 0.0
      end do
      do i=1,loc
         ymin = amin1(ymin,yebot(i))
         ymax = amax1(ymax,yetop(i))
         smooth(i) = biw(i) - fit(1) - fit(2) * i
      end do
      call realft (smooth,nbin,1)
      smooth(1) = smooth(1) / nbin
      const = (smlen/(2.0*nbin))**2
      fac = 1.0
      do j=1,nbin-1
         k = 2 * j + 1
         if (fac .ne. 0.0) then
            fac = amax1(0.0,(1.0-const*j*j)/nbin)
            smooth(k) = fac * smooth(k)
            smooth(k+1) = fac * smooth(k+1)
         else
            smooth(k) = 0.0
            smooth(k+1) = 0.0
         end if
      end do
      fac = amax1(0.0,(1.0-0.25*smlen*smlen)/nbin)
      smooth(2) = fac * smooth(2)
      call realft (smooth,nbin,-1)
      do i=1,loc
         smooth(i) = smooth(i) + fit(1) + fit(2) * i
         ymin = amin1(ymin,smooth(i))
         ymax = amax1(ymax,smooth(i))
      end do


c  plot results

      call pgenv (0.0,loc+1.0,ymin,ymax,0,1)
      call pglabel (" "," ",filename)
      call pgerry (loc,xxx,yetop,yebot,1.0)
      call pgpoint (loc,xxx,biw,0)
      call pgline (loc,xxx,smooth)

      write (*,1221)
 1221 format ('Change smoothing length (y/n)? '$)
      read (*,1003) answer
      if (answer .eq. 'y') go to 220
      

 221  write (*,1222)
 1222 format ("Enter filename for subtracted image: "$)
      read (*,1003) outfile
      call imopnc (outfile, im, jm, ier)
      if (ier .ne. 0) then
         write (0,1004) bell
         go to 221
      end if

      call spline (xxx,smooth,loc,spl)

      do line=1,ny
         call imgl2r (im, dat, line, ier)
         arad = (line - ayc)**2
         rrad = (line - ryc)**2
         do ipix=1,nx
            if (((ipix - axc)**2 + arad) .le. armxsq) then
               if (dat(ipix) .gt. -32000.) then
                  rrsq = rrad + (ipix - rxc)**2
                  ar = nbin * rrsq / rrmxsq + 0.5
                  if (ar .gt. loc) ar = loc
                  call splint (xxx,smooth,spl,loc,ar,sky)
                  dat(ipix) = dat(ipix) - sky
               end if
            end if
         end do
         call impl2r (jm,dat,line,ier)
      end do

      call imclos (im, ier)
      if (ier .ne. 0) then
         write (*,1009) bell
 1009    format (a1,"Error closing image")
      end if

      call imclos (jm, ier)
      if (ier .ne. 0) then
         write (*,1009) bell
      end if

      write (*,1010)
 1010 format ("Enter filename for output: "$)
      read (*,1003) filename
      open (10,file=filename,form="formatted")
      do ir=1,loc
         rad = sqrt(rrmxsq * ir / nbin)
         write (10,1011) num(ir),rad,wave(ir),biw(ir),
     $                   sb(ir),smooth(ir)
 1011    format (i5,f7.2,f10.3,3f10.2)
      end do
      close (10)

C	done - check if more rings

      write (*,"('More rings (y/n)? '$)")
      read (*,*) answer
      if ((answer.eq.'y').or.(answer.eq.'Y')) go to 10

      call pgend ()

      stop
      end
 
