c------------------------------------------------------------------------------
      subroutine ringfilter(dir, axc, ayc, arad,rxcen,rycen,
     >     icut, iwide, itmax, converg, wscale, filename)
c------------------------------------------------------------------------------

c  This program determines the center coordinates of a ring, bins the
c  ring radially and computes its power spectrum, and allows the user 
c  to select a smoothing filter for the ring.

c nsl:
c 12.05.10 replaced all references to imfort routines with fitsio routines
c replaced isize with naxes
c replaced calls to getrfp.f with getpfp.f
c 06.07.10 updated to make compatible with SALT pyraf package

      implicit none

      integer   NRING, NSEG, MAXPTS, MAXIN, MAXDIV
c      parameter (NRING=256, NSEG=8, MAXPTS=512, MAXIN=2048, MAXDIV=3)
c updated 12.05.10 by nic 
      parameter (NRING=256, NSEG=8, MAXPTS=512, MAXIN=3172, MAXDIV=3)
      integer   ISHORT, IREAL, IRW, IRO
      parameter (ISHORT=3, IREAL=6, IRW=3, IRO=1)
      integer   POINTCOLOR, LINECOLOR
      parameter (POINTCOLOR=7, LINECOLOR=5)

      real dat(MAXIN), newdata(MAXPTS,NRING,NSEG)
      real x(NRING), y(NRING+2), dy(NRING), pl(NRING+2), fr(NRING)
      real ryp(NRING), ryh(NRING), ryl(NRING), drp(NRING)
      real pl2(2*NRING + 4), ddy
      real a(5), da(5), dvda(5), rad(NSEG)
      real xmin, xmax, frmin, frmax, axc, ayc, arad, rxcen, rycen
      real converg, wscale, arsq, rxc, ryc, rmax, rmaxsq, rx, ry, rysq
      real dx, rsq, chisq, xcen, ycen, fwhm, cont, width, oldstep, rave
      real rerr, dxc, dyc, ang, step, rymin, rymax, vgt, pmin, pmax
      real pamp, arg, yscale
      integer numb(NRING,NSEG), isam(NRING,NSEG), isize(7)
      integer i, j, ier, itmax, isample, jsample, lnblnk, ibeg, iend
      integer img, naxis, itype, iylo, iyhi, iy, ixlo, ixhi, ix, ir
      integer iter, ndiv, iseg, nhalf, loc, icut, iwide, ilo, ihi, count
      logical imgopen, flag(5), query, conver
      character*80 filename, param, dir
      character*1 bell
      integer blocksize, inunit, readwrite, idim, keyval, bitpix
      character*24 comment
      integer npixels, naxes(2), firstpix, nullvalj, group
      real nullvale
      logical anynull

c  initialization
      
      count = 0

      bell = char(7)
      imgopen = .false.
      do i=1,5
         flag(i) = .true.
      end do
      do i=1,NRING
         x(i) = (i - 0.5) / NRING
         fr(i) = i - 1.0
      end do
      xmin = 0.0
      xmax = 1.0
      frmin = 0.0
      frmax = NRING ! / 2.0
      write (*, '(/"Ring Power Spectrum Filtering"//)')


c  get the options for the ring fitting operation

c      call getpfp ("axc", param)
c      read (param, *, iostat=ier) axc
c      do while (ier .ne. 0)
c         write (*, '("Enter aperture x center: "$)')
c         read (*, *, iostat=ier) axc
c      end do
c
c      call getpfp ("ayc", param)
c      read (param, *, iostat=ier) ayc
c      do while (ier .ne. 0)
c         write (*, '("Enter aperture y center: "$)')
c         read (*, *, iostat=ier) ayc
c      end do
c
c      call getpfp ("arad", param)
c      read (param, *, iostat=ier) arad
c      do while (ier .ne. 0)
c         write (*, '("Enter aperture radius: "$)')
c         read (*, *, iostat=ier) arad
c      end do
c
c      call getpfp ("rxc", param)
c      read (param, *, iostat=ier) rxcen
c      do while (ier .ne. 0)
c         write (*, '("Enter ring x center: "$)')
c         read (*, *, iostat=ier) rxcen
c      end do
c
c      call getpfp ("ryc", param)
c      read (param, *, iostat=ier) rycen
c      do while (ier .ne. 0)
c         write (*, '("Enter ring y center: "$)')
c         read (*, *, iostat=ier) rycen
c      end do
c
c      call getpfp ("calring_itmax", param)
c      read (param, *, iostat=ier) itmax
c      if (ier .ne. 0) itmax = 10
c
c      call getpfp ("calring_conv", param)
c      read (param, *, iostat=ier) converg
c      if (ier .ne. 0) converg = 0.01
c
c      call getpfp ("calring_fitwidth", param)
c      read (param, *, iostat=ier) wscale
c      if (ier .ne. 0) wscale = 0.50
c
      arsq = arad * arad
      isample = 1.10 + 3.14159 * arsq  / (NSEG * MAXPTS * NRING)
      jsample = 1.10 + 3.14159 * arsq  / (MAXPTS * NRING)
c forcing plot to xw (nsl)
      call pgbegin (0,"/xw",2,2)
      call pgask (.false.)


c  get the next ring image

      if (lnblnk(dir) .gt. 0) then
         ier = chdir(dir)
         if (ier .ne. 0) write (*, '("Error changing directory")')
      else
         ier = 0
      end if

c 100  filename = " "
c      do while (lnblnk(filename) .le. 0)
c         write (*, '("Enter ring filename:  "$)')
c         read (*, '(a)') filename
c      end do
      call strlim (filename, ibeg, iend)

C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(inunit,ier)
      blocksize=1 ! this is ignored
      ier=0
      readwrite=0 ! read only
      call ftopen(inunit,filename,readwrite,blocksize,ier)
c      call imopen (filename, IRO, img, ier)
      if (ier .ne. 0) then
         write (*, '(a,"Error opening ring image ")') bell
         go to 300
      end if
      imgopen = .true.

      call ftgknj(inunit,'NAXIS',1,2,naxes,idim,ier)
c      call imgsiz (img, isize, naxis, itype, ier)
      if (ier .ne. 0) then
         write (*, '(a,"Error determining image size")') bell
         go to 300
      end if

c      if ((itype .ne. ISHORT) .and. (itype .ne. IREAL)) then
      call FTGKYJ(inunit,'BITPIX', keyval,comment,ier)
      bitpix=keyval

      call FTGKYJ(inunit,'NAXIS', keyval,comment,ier)
      naxis=keyval

      if ((bitpix .ne. 16) .and. (bitpix .ne. -32)) then
         write (*, '(a,"Wrong image type ",a)') 
     >     bell, filename(ibeg:iend)
         go to 300
      end if

      if (naxis .ne. 2) then
         write (*, '(a,"Wrong image dimension: ",i2)')
     $        bell, naxis
         go to 300
      end if

      if (naxes(1) .gt. MAXIN) then
         write (*, '(a,"Ring image too large: ",i2)') 
     >        bell, naxes(1)
         go to 300
      end if

      print*, naxes

c  determine initial ring parameters

      rxc = rxcen
      ryc = rycen

      do i=1,NRING
         isam(i,1) = 0
         numb(i,1) = 0
      end do

      rmax = arad - sqrt((axc - rxc)**2 + (ayc - ryc)**2)
      rmaxsq = rmax * rmax
      iylo = max(nint(ryc - rmax), 1)
      iyhi = min(nint(ryc + rmax), naxes(2))

C  Initialize variables for image reading

      npixels=naxes(1) ! no in x direction
      firstpix=1
      nullvalj=-999
      nullvale=-999.0
      group=1


      do iy=iylo,iyhi
c         call imgl2r (img, dat, iy, ier)
         firstpix= 1 + (iy-1)*(npixels)
         call ftgpve(inunit,group,firstpix,npixels,nullvale,
     &        dat,anynull,ier)

         if (ier .ne. 0) then
            write (*, '(a,"Error while reading ring file")') bell
            go to 300
         end if

         ry = iy - ryc
         rysq = ry * ry
         dx = sqrt(rmaxsq - rysq)
         ixlo = max (nint(rxc - dx), 1)
         ixhi = min (nint(rxc + dx), naxes(1))

         do ix=ixlo,ixhi
            rx = ix - rxc
            rsq = rysq + rx * rx
            ir = 1.0 + NRING * rsq / rmaxsq
            if (ir .le. NRING) then
               isam(ir,1) = isam(ir,1) + 1
               if (mod(isam(ir,1),jsample) .eq. 0) then
                  numb(ir,1) = numb(ir,1) + 1
                  newdata(numb(ir,1),ir,1) = dat(ix)
               end if
            end if
         end do
      end do
 
      do i=1,NRING
         call biwgt (newdata(1,i,1),numb(i,1),y(i),dy(i))
      end do

      call evinit (x, y, NRING, a)
      call evfit (x, y, dy, NRING, a, flag, da, chisq)
      call evstat (a, ycen, fwhm)
      xcen = a(3) * NRING
      cont = a(1)
      width = fwhm * NRING * wscale


c  loop to measure centers and radius

      iter = 0
      conver = .false.
      oldstep = 1000.0
      ndiv = 0
      write (*, '("Measuring ring center "$)')
      call flush (6)

      do while ((iter .lt. itmax) .and. (.not. conver) .and. 
     $          (ndiv .lt. MAXDIV))

         do j=1,NSEG
            do i=1,NRING
               isam(i,j) = 0
               numb(i,j) = 0
            end do
         end do

         rmax = arad - sqrt ((axc - rxc)**2 + (ayc - ryc)**2)
         rmaxsq = rmax * rmax
         iylo = max (nint(ryc - rmax), 1)
         iyhi = min (nint(ryc + rmax), naxes(2))

         do iy=iylo,iyhi

c            call imgl2r (img, dat, iy, ier)
            firstpix= 1 + (iy-1)*(npixels)
            call ftgpve(inunit,group,firstpix,npixels,nullvale,
     &           dat,anynull,ier)
            if (ier .ne. 0) then
               write (*, '(a1,"Error while reading ring file")') bell
               go to 300
            end if

            ry = iy - ryc
            rysq = ry * ry
            dx = sqrt (rmaxsq - rysq)
            ixlo = max (nint(rxc - dx), 1)
            ixhi = min (nint(rxc + dx), naxes(1))

            do ix=ixlo,ixhi

               rx = ix - rxc
               rsq = rysq + rx * rx
               ir = 1.0 + NRING * rsq / rmaxsq
               if (ir .le. NRING) then

                  iseg = (atan2 (ry, rx) + 3.14159) * NSEG / 6.28319 + 1
                  iseg = max (iseg, 1)
                  iseg = min (iseg, NSEG)
               
                  isam(ir,iseg) = isam(ir,iseg) + 1
                  if (mod(isam(ir,iseg),isample) .eq. 0) then
                     numb(ir,iseg) = numb(ir,iseg) + 1
                     newdata(numb(ir,iseg),ir,iseg) = dat(ix)
                  end if

               end if
            end do
         end do
 
         rave = 0.0
         rerr = 0.0
         dxc  = 0.0
         dyc  = 0.0
         
         do j=1,NSEG
            
            do i=1,NRING
               call biwgt (newdata(1,i,j), numb(i,j), y(i), dy(i))
            end do
            
            call center (y, NRING, xcen, width, cont)
            if (xcen .lt.0.0) then
               write (*, '(a, "Error fitting ring center")') bell
               go to 300
            end if

            rad(j) = rmax * sqrt ((xcen - 0.5) / NRING)
            rave = rave + rad(j)
            rerr = rerr + rad(j) * rad(j)
            ang = 6.28319 * (j - 0.5) / NSEG
            dxc = dxc + rad(j) * cos (ang)
            dyc = dyc + rad(j) * sin (ang)

         end do

         dxc = 2.0 * dxc / NSEG
         dyc = 2.0 * dyc / NSEG
         rxc = rxc - dxc
         ryc = ryc - dyc
         rave = rave / NSEG
         rerr = rerr / NSEG - rave * rave
         if (rerr .gt. 0.0) then
            rerr = sqrt (rerr)
         else
            rerr = 0.005
         end if
         step = sqrt (dxc**2 + dyc**2)

         conver = step .lt. (converg * max (1.0, rerr))

         if (step .lt. oldstep) then
            oldstep = step
            ndiv = 0
         else
            ndiv = ndiv + 1
         end if

         iter = iter + 1
         write (*, '("."$)')
         call flush (6)

      end do
      write (*, '(" done")')
      write (*, '("Ring X center: ", f7.2, 2x, "Y center: ", f7.2,
     $     2x, "error: ", f4.2)') rxc, ryc, rerr


c  get final ring profile

      do i=1,NRING
         isam(i,1) = 0
         numb(i,1) = 0
      end do

      rmax = arad - sqrt((axc - rxc)**2 + (ayc - ryc)**2)
      rmaxsq = rmax * rmax
      iylo = max(nint(ryc - rmax), 1)
      iyhi = min(nint(ryc + rmax), naxes(2))

      do iy=iylo,iyhi
            firstpix= 1 + (iy-1)*(npixels)
            call ftgpve(inunit,group,firstpix,npixels,nullvale,
     &           dat,anynull,ier)

c         call imgl2r (img, dat, iy, ier)
         if (ier .ne. 0) then
            write (*, '(a,"Error while reading ring file")') bell
            go to 300
         end if
         
         ry = iy - ryc
         rysq = ry * ry
         dx = sqrt(rmaxsq - rysq)
         ixlo = max (nint(rxc - dx), 1)
         ixhi = min (nint(rxc + dx), naxes(1))

         do ix=ixlo,ixhi

            rx = ix - rxc
            rsq = rysq + rx * rx
            ir = 1.0 + NRING * rsq / rmaxsq
            if (ir .le. NRING) then
               isam(ir,1) = isam(ir,1) + 1
               if (mod(isam(ir,1),jsample) .eq. 0) then
                  numb(ir,1) = numb(ir,1) + 1
                  newdata(numb(ir,1),ir,1) = dat(ix)
               end if
            end if
         end do
      end do
 
      do i=1,NRING
         call biwgt (newdata(1,i,1),numb(i,1),ryp(i),drp(i))
         if (numb(i,1) .gt. 1) drp(i) = drp(i) / sqrt(numb(i,1) - 1.0)
         ryh(i) = ryp(i) + drp(i)
         ryl(i) = ryp(i) - drp(i)
      end do

      call evinit (x, ryp, NRING, a)
      call evfit (x, ryp, drp, NRING, a, flag, da, chisq)

      rymin = ryl(1)
      rymax = ryh(1)
      do i=2,NRING
         rymin = min(rymin, ryl(i))
         rymax = max(rymax, ryh(i))
      end do
      yscale = rymax - rymin
      rymin = rymin - 0.05 * yscale
      rymax = rymax + 0.05 * yscale


c determine filter params

 200  call pgsch (1.25)
      call pgenv (xmin, xmax, rymin, rymax, 0, 1)
      call pglabel ("Radius\u2\d", "Intensity", "Raw Profile")
      call pgsch (1.00)
      call pgsci (POINTCOLOR)
      call pgpoint (NRING, x, ryp, 17)
      call pgerry (NRING, x, ryh, ryl, 1.0)
      call pgsci (1)

      do i=1,NRING
         pl(i) = ryp(i)
      end do

      nhalf = NRING / 2


      do i=1,2*NRING + 4
         pl2(i)=0.0
      enddo
      do i=1,NRING 
         pl2(i)=pl(i)
      enddo      

c      call realft (pl, nhalf, 1)
      call realft (pl2, NRING, 1)

      do i=1,NRING
         pl(i)=pl2(i)
      enddo

      pl(1) = alog10(pl(1))
      pmin = pl(1)
      pmax = pl(1)
      loc = 3
      do i=2,nhalf
         pl(i) = alog10(sqrt(pl(loc)**2 + pl(loc+1)**2))
         pmin = min(pmin, pl(i))
         pmax = max(pmax, pl(i))
         loc = loc + 2
      end do
      pmax = pmax + 0.1
      pmin = pmin - 0.1

      call pgsch (1.25)
      call pgenv (frmin, frmax/2, pmin, pmax, 0, 20)
      call pglabel ("Frequency", "Power", "Power Spectrum")
      call pgsch (1.00)
      call pgsci (POINTCOLOR)
      call pgpoint (nhalf, fr, pl, 17)
      call pgsci (1)

      if (count .gt. 0) then
         write (*, '("Enter filter cut-off frequency and width: "$)')
         read (*,*) icut, iwide
      endif

c      loc = min (max(icut, 2), nhalf)
c      ilo = max (loc - iwide, 2)
c      ihi = min (loc + iwide, nhalf)
c      step = 3.14159 / (ihi - ilo + 2.0)

      loc = min (max(icut, 2), NRING)
      ilo = max (loc - iwide, 2)
      ihi = min (loc + iwide, NRING)
      step = 3.14159 / (ihi - ilo + 2.0)

      do i=1,ilo-1
         pl(i) = pmax - 0.05
      end do
      pamp = 0.5 * (pmax - pmin - 0.10)
      arg = step
      do i=ilo,ihi
         pl(i) = pmin + 0.05 + pamp * ( 1.0 + cos(arg))
         arg = arg + step
      end do
c      if (ihi .lt. nhalf) then
c         do i=ihi+1,nhalf
      if (ihi .lt. NRING) then
         do i=ihi+1,NRING
            pl(i) = pmin + 0.05
         end do
      end if
      call pgsci (LINECOLOR)
      call pgslw (2)
      call pgline (NRING, fr, pl)
      call pgsci (1)
      call pgslw (1)

c this is where the lowpass filtering takes place - lowpass has been edited to make it work properly

      do i=1,NRING
         call evoigt (x(i), a, 5, vgt, dvda)
         y(i) = ryp(i) - vgt
      end do
      call lowpass (y, NRING, icut, iwide)
      do i=1,NRING
         call evoigt (x(i), a, 5, vgt, dvda)
         y(i) = y(i) + vgt
      end do

c plot original and smoothed data

      print*, xmin, xmax, rymin, rymax
      call pgenv (xmin, xmax, rymin, rymax, 0, 1)
      call pglabel ("Radius\u2\d", "Intensity", "Raw Profile")
      call pgsch (1.00)
      call pgsci (POINTCOLOR)
      call pgpoint (NRING, x, ryp, 17)
      call pgerry (NRING, x, ryh, ryl, 1.0)
      call pgsci (1)

      call pgsch (1.25)
      call pgenv (xmin, xmax, rymin, rymax, 0, 1)
      call pglabel ("Radius\u2\d", "Intensity", "Filtered Profile")
      call pgsch (1.00)
      call pgsci (POINTCOLOR)
      call pgpoint (NRING, x, y, 17)
      call pgsci (1)

      ier = 1
      do while (ier .ne. 0)
         write (*, '("Change the filter (y/n)? "$)')
         read (*, '(a)') param
         call yesno (param, query, ier)
      end do

      if (query) then
         count = count + 1
         go to 200
      endif

      write (*, '("Filter frequency: ", i3, 3x, "width: " i3)') 
     $     icut, iwide

c  ring completed - loop for more

 300  if (imgopen) then
         call ftclos (inunit, ier)
         imgopen = .false.
         call ftfiou(inunit, ier)
      endif

c      ier = 1
c      do while (ier .ne. 0)
c         write (*, '("Do another ring (y/n)? "$)')
c         read (*, '(a)') param
c         call yesno (param, query, ier)
c      end do
c      if (query) go to 100


c  done - close up

      write (*, '("Processing completed")')
      call pgend ()

      
      end subroutine ringfilter



