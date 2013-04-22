c------------------------------------------------------------------------------
      subroutine nightring(dir,nightringlogfile, outfilename, label, 
     >     axc,ayc,arad,rxcen,rycen,filter,icut, iwide,
     >     plot, plottype,itmax,converg, wscale,
     >     filename,
     >     cala, calb, calc, cald, calf,verb)
c------------------------------------------------------------------------------

c  This program measures the center coordinates and radius of one or more
c  night-time calibration rings.  It calculates the current wavelength
c  zero point (A) for use in tracking the etalon drifts through the night.

c nsl:
c 12.05.10 replaced all references to imfort routines with fitsio routines
c replaced isize with naxes
c replaced calls to getrfp to getpfp
c ut read as a string from the headers rather than real
c Updated to read which etalon is in use and take correct Z from fits header.
c 010710, updated to be compatible with SALT pyraf package. 

      implicit none

      integer   NRING, NSEG, MAXPTS, MAXIN, MAXDIV
c      parameter (NRING=256, NSEG=8, MAXPTS=512, MAXIN=2048, MAXDIV=3)
c updated 12.05.10 by nic
      parameter (NRING=256, NSEG=8, MAXPTS=512, MAXIN=3172, MAXDIV=3)

      integer   ISHORT, IREAL, IRW, IRO
      parameter (ISHORT=3, IREAL=6, IRW=3, IRO=1)

      real dat(MAXIN), newdata(MAXPTS,NRING,NSEG)
      real x(NRING), y(NRING+2), rry(NRING), ryl(NRING), ryh(NRING)
      real a(5), da(5), rad(NSEG), xv(2), yv(2), dy(NRING)
      real axc, ayc, arad, arsq, rxcen, rycen, rxc, ryc, converg, wscale
      real rmax, rmaxsq, ry, rysq, dx, rx, rsq, chisq, ycen, fwhm, xcen
      real cont, width, oldstep, rave, rerr, dxc, dyc, ang, step
      real ymin, ymax,z, wave
      real cala, calb, calc, cald, calf, anew
      integer numb(NRING,NSEG), isam(NRING,NSEG)
      integer i,ier, ibeg, iend, icut, iwide, itmax, isample
      integer jsample, naxis, iylo, iyhi, iy, ixlo, ixhi, ix
      integer ir, iter, ndiv, j, iseg, izval, lnblnk
      logical imgopen,flag(5), query, plot, list, conver, verb, filter
      logical anynull
      character*80 filename, label,outfilename,nightringlogfile
      character*80 dir
      character*40 param
      character*24 datetime, comment
      character*1 bell
      character*4 plottype
      integer      blocksize, readwrite, inunit, keyval, naxes(2)
      integer      group, nullvalj, npixels, firstpix, bitpix 
      real         nullvale
      character*4 etalonz
      character*20 etstate
      character*8 etalonl
      real ut
      integer iyr, imo, ida, ihr, imn
      real*8 secs


c      write(*, '(a)')  dir, nightringlogfile, outfilename, label
c      print*, axc, ayc, rxcen, rycen, icut, iwide
c      print*, filter, plot
c      write(*, '(a)') plottype, filename
c      print*, itmax, converg, wscale, verb
c      print*, cala, calb, calc, cald, calf

c  initialization

      bell = char(7)
      imgopen = .false.
      do i=1,5
         flag(i) = .true.
      end do
      do i=1,NRING
         x(i) = (i - 0.5) / NRING
      end do
c      nargs = iargc ()
      write (*, '(/"Calibration Ring Fitting"//)')

c nic added this to change to directory containing the data
      ier = 1
      if (lnblnk(dir) .gt. 0) then
         ier = chdir(dir)
         if (ier .ne. 0) write (*, '("Error changing directory")')
      else
         ier = 0
      end if      

c  start the log file

      open (10, file=nightringlogfile, access="append", iostat=ier)
      if (ier .ne. 0) then
         write (*, '(a,"Error opening log file")') bell
         stop
      end if
      write (10, '(/"Calibration Ring Fitting")')
      call fdate (datetime)
      write (10, '(a24)') datetime


c  open the output file

c      if (nargs .ge. 2) then
c         call getarg (2, outfilename)
c      else
c         outfilename = " "
c         do while (lnblnk(outfilename) .le. 0)
c            write (*, '("Enter a filename for the output: "$)')
c            read (*, '(a)') outfilename
c         end do
c      end if

      call strlim (outfilename, ibeg, iend)
      inquire (file=outfilename, exist=query)
      open (11, file=outfilename, access="append", iostat=ier)
      if (ier .ne. 0) then
         write (*, '(a,"Error opening output file ",a)')
     $        bell, outfilename(ibeg:iend)
         write (10, '("Error opening output file ",a)')
     $        outfilename(ibeg:iend)
         go to 666
      end if
      write (10, '("Output file: ",a)') outfilename(ibeg:iend)

      if (.not. query) then
c         if (nargs .ge. 2) then
c            label = datetime
c         else
c            write (*, '("Enter a comment for the output file:")')
c            read (*, '(a)') label
            if (lnblnk(label) .le. 0) then
               label = datetime
            end if
c         end if

         write (11, '(a)') label(1:lnblnk(label))
         write (11, '(" radius  err   xc     yc      z     ut",
     $        "      wave      a     file")')
      end if


c  get the options for the ring fitting operation

c      call getpfp ("axc", param)
c      read (param, *, iostat=ier) axc
c      do while (ier .ne. 0)
c         write (*, '("Enter aperture x center: "$)')
c         read (*, *, iostat=ier) axc
c      end do

c      call getpfp ("ayc", param)
c      read (param, *, iostat=ier) ayc
c      do while (ier .ne. 0)
c         write (*, '("Enter aperture y center: "$)')
c         read (*, *, iostat=ier) ayc
c      end do

c      call getpfp ("arad", param)
c      read (param, *, iostat=ier) arad
c      do while (ier .ne. 0)
c         write (*, '("Enter aperture radius: "$)')
c         read (*, *, iostat=ier) arad
c      end do

c      call getpfp ("rxc", param)
c      read (param, *, iostat=ier) rxcen
c      do while (ier .ne. 0)
c         write (*, '("Enter ring x center: "$)')
c         read (*, *, iostat=ier) rxcen
c      end do

c      call getpfp ("ryc", param)
c      read (param, *, iostat=ier) rycen
c      do while (ier .ne. 0)
c         write (*, '("Enter ring y center: "$)')
c         read (*, *, iostat=ier) rycen
c      end do

c      call getpfp ("calibration_a", param)
c      read (param, *, iostat=ier) cala
c      do while (ier .ne. 0)
c         write (*, '("Enter calibration A value: "$)')
c         read (*, *, iostat=ier) cala
c      end do

c      call getpfp ("calibration_b", param)
c      read (param, *, iostat=ier) calb
c      do while (ier .ne. 0)
c         write (*, '("Enter calibration B value: "$)')
c         read (*, *, iostat=ier) calb
c      end do

c      call getpfp ("calibration_c", param)
c      read (param, *, iostat=ier) calc
c      do while (ier .ne. 0)
c         write (*, '("Enter calibration C value: "$)')
c         read (*, *, iostat=ier) calc
c      end do

c      call getpfp ("calibration_d", param)
c      read (param, *, iostat=ier) cald
c      do while (ier .ne. 0)
c         write (*, '("Enter calibration D value: "$)')
c         read (*, *, iostat=ier) cald
c      end do

c      call getpfp ("calibration_f", param)
c      read (param, *, iostat=ier) calf
c      do while (ier .ne. 0)
c         write (*, '("Enter calibration F value: "$)')
c         read (*, *, iostat=ier) calf
c      end do

c      call getpfp ("calring_filter", param)
c      call yesno (param, filter, ier)
c      do while (ier .ne. 0)
c         write (*, '("Filter the profiles (y/n)? "$)')
c         read (*, '(a)') param
c         call yesno (param, filter, ier)
c      end do

c      if (filter) then

c         call getpfp ("calring_filter_freq", param)
c         read (param, *, iostat=ier) icut
c         do while (ier .ne. 0)
c            write (*, '("Enter filter frequency: "$)')
c            read (*, *, iostat=ier) icut
c         end do

c         call getpfp ("calring_filter_width", param)
c         read (param, *, iostat=ier) iwide
c         do while (ier .ne. 0)
c            write (*, '("Enter filter cutoff width: "$)')
c            read (*, *, iostat=ier) iwide
c         end do

c      end if

c      call getpfp ("calring_plot", param)
c      call yesno (param, plot, ier)
c      do while (ier .ne. 0)
c         write (*, '("Plot the fits (y/n)? "$)')
c         read (*, '(a)') param
c         call yesno (param, plot, ier)
c      end do

c      call getpfp ("calring_verbose", param)
c      call yesno (param, verb, ier)
c      if (ier .ne. 0) verb = .false.

c      call getpfp ("calring_itmax", param)
c      read (param, *, iostat=ier) itmax
c      if (ier .ne. 0) itmax = 10

c      call getpfp ("calring_conv", param)
c      read (param, *, iostat=ier) converg
c      if (ier .ne. 0) converg = 0.01

c      call getpfp ("calring_fitwidth", param)
c      read (param, *, iostat=ier) wscale
c      if (ier .ne. 0) wscale = 0.50

      arsq = arad * arad
      isample = 1.10 + 3.14159 * arsq  / (NSEG * MAXPTS * NRING)
      jsample = 1.10 + 3.14159 * arsq  / (MAXPTS * NRING)
      if (plot) then
         call pgbegin (0,plottype,1,1)
c         call pgbegin (0,"?",1,1)
         call pgask (.false.)
      end if


c  get the input filename (or list)

c 100  if (nargs .ge. 1) then
c         call getarg (1, filename)
c      else
c         filename = " "
c         do while (lnblnk(filename) .le. 0)
c            write (*, '("Enter ring filename or list: "$)')
c            read (*, '(a)') filename
c         end do
c      end if

      call strlim (filename, ibeg, iend)
      if (filename(ibeg:ibeg) .eq. "@") then
         open (12, file=filename(ibeg+1:iend), status="old", 
     $         iostat=ier)
         if (ier .ne. 0) then
            write (10, '("Error opening image list ",a)')
     $           filename(ibeg+1:iend)
            write (*, '(a,"Error opening image list ",a)')
     $           bell, filename(ibeg+1:iend)
c            if (nargs .ge. 1) then
c               go to 666
c            else
c               go to 100
c            end if
         end if
         list = .true.
      else
         list = .false.
      end if


c  open next ring image

 200  if (list) then
         read (12, '(a)', iostat=ier) filename
         if (ier .ne. 0) then
            if (ier .gt. 0) then
               write (*, '(a, "Error reading ring filename from list")')
     $              bell
               write (10, '("Error reading ring filename from list")')
            end if
            close (12)
            list = .false.
            go to 300
         end if
         call strlim (filename, ibeg, iend)
         write (*, '("Processing ring: ",a)') filename(ibeg:iend)
      end if
      write (10, '("Processing ring: ",a)') filename(ibeg:iend)

      readwrite=0 !read only
      blocksize=1 ! this is ignored

      call ftgiou(inunit,ier) ! get a LUN
      call ftopen(inunit,filename,readwrite,blocksize,ier)
c      call imopen (filename, IRO, img, ier)
      if (ier .ne. 0) then
         write (*, '(a,"Error opening ring image ")') bell
         write (10, '("Error opening ring image ")')
         go to 300
      end if
      imgopen = .true.

      call ftgknj(inunit,'NAXIS',1,2,naxes,naxis,ier)
c      call imgsiz (img, isize, naxis, itype, ier)
      if (ier .ne. 0) then
         write (*, '(a,"Error determining image size")') bell
         write (10, '("Error determining image size")')
         go to 300
      end if

c     if ((itype .ne. ISHORT) .and. (itype .ne. IREAL)) then
      call FTGKYJ(inunit,'BITPIX', keyval,comment,ier)
      bitpix=keyval
      if ((bitpix .ne. 16) .and. (bitpix .ne. -32)) then      
         write (*, '(a,"Error: wrong image type: ",i1)') bell, bitpix
         write (10, '("Error: wrong image type: ",i1)') bitpix
         go to 300
      end if

      if (naxis .ne. 2) then
         write (*, '(a,"Error: wrong image dimension: ",i2)')
     $        bell, naxis
         write (10, '("Error: wrong image dimension: ",i2)') naxis
         go to 300
      end if

      if (naxes(1) .gt. MAXIN) then
         write (*, '(a,"Ring image too large: ",i2)') bell, naxes(1)
         write (10, '("Ring image too large: ",i2)') naxes(1)
         go to 300
      end if


c  determine initial ring parameters

      rxc = rxcen
      ryc = rycen
c      print*, rxc, ryc, 'initial centres'

      do i=1,NRING
         isam(i,1) = 0
         numb(i,1) = 0
      end do

      rmax = arad - sqrt((axc - rxc)**2 + (ayc - ryc)**2)
      rmaxsq = rmax * rmax
      iylo = max(nint(ryc - rmax), 1)
      iyhi = min(nint(ryc + rmax), naxes(2))

c added nsl 05/07/10
      iy=0

C  Initialize variables for image reading

      npixels=naxes(1) ! no in x direction
      firstpix=1
      nullvalj=-999
      nullvale=-999.0
      group=1


      do iy=iylo,iyhi
         ier=0
         firstpix= 1 + (iy-1)*(npixels)
         call ftgpve(inunit,group,firstpix,npixels,nullvale,
     &        dat,anynull,ier)

c         call imgl2r (img, dat, iy, ier)
         if (ier .ne. 0) then
            write (*, '(a,"Error while reading ring file")') bell
            write (10, '("Error while reading ring file")')
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

      if (filter) call lowpass (y, NRING, icut, iwide)

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

c         print*, npixels, 'npixels', maxin
         do iy=iylo,iyhi
            ier=0
            firstpix= 1 + (iy-1)*(npixels)
            call ftgpve(inunit,group,firstpix,npixels,nullvale,
     &           dat,anynull,ier)
c            call imgl2r (img, dat, iy, ier)
            if (ier .ne. 0) then
               write (*, '(a, "Error while reading ring file")') bell
               write (10, '("Error while reading ring file")')
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
            
            if (filter) call lowpass (y, NRING, icut, iwide)
            
            call center (y, NRING, xcen, width, cont)
            if (xcen .lt. 0.0) then
               write (*, '(a, "Error fitting ring center")') bell
               write (10, '("Error fitting ring center")')
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

         if (verb) then
            if (iter .eq. 0) then
               write (*, '("  radius   error   x center   y center",
     $              "    step")')
            end if
            write (*, '(1x, f7.2, 3x, f5.2, 3x, 2(f7.2, 4x), f5.2)')
     $           rave, rerr, rxc, ryc, step
         end if


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


c  write output data

      write (*, '("Radius: ",f7.2, " +/- ", f6.2, 5x, "X center: ",
     $     f7.2, 5x, "Y center: ", f7.2)') rave, rerr, rxc, ryc
      

c read these values from the ring file
      ier=0
      call ftgkys (inunit, "UTC-OBS", param, comment, ier)

      if (ier .eq. 0) then
         call fts2tm (param,iyr,imo,ida,ihr,imn,secs,ier)
      end if
      ut = real(ihr) + real(imn)/60.0 + secs/3600.0

      if (ier .ne. 0) then
         ut = -1.0
      end if


c      call imgkwi (img, "FPZ", izval, ier)

      call FTGKYS(inunit,"ET-STATE", etstate,comment,ier)
      write(*,'(a)') etstate(6:7)
      if ( (etstate(1:2) .eq. 'S2') .or.
     >    (etstate(6:7) .eq. 'S2') ) then                !Et1
         etalonz='ET1Z'
         etalonl='ET1WAVE0'
      else if ( (etstate(1:2) .eq. 'S3') .or. 
     >        (etstate(6:7) .eq. 'S3') ) then              !Et2
         etalonz='ET2Z'
         etalonl='ET2WAVE0'
      else
         print*, 'problem with which etalon to choose'
         etalonz='ETZ' ! this doesnt exist, forces user input
         etalonl='ETWAVE0' ! this doesnt exist, forces user input
      endif

      call ftgkyj (inunit, etalonz, izval,comment,ier)
      if (ier .ne. 0) then
         write (*, '("Enter Z value for this ring: "$)')
         read (*,*) izval
      end if
      z = izval / 1000.0

      ier=0
      call ftgkye (inunit, etalonl, wave, comment, ier)
      do while (ier .ne. 0)
         write (*, '("Enter wavelength (", f8.3, "): "$)') wave
         read (*, '(a)') param
         if (lnblnk(param) .gt. 0) then
           read (param, *, iostat=ier) wave
         else
            ier = 0
         end if
      end do

c  calculate new zero point

      anew = wave * sqrt (1.0 + (rave / calf)**2) - z * (calb + z * 
     $        (calc + z * cald))
      write (*, '("New fpwave0 value: ", f8.3)') anew
      write (11, '(1x, f6.2, 1x, f4.2, 2(1x, f7.2), 1x, f6.3, 1x, 
     $        f10.5, 1x, f8.3, 1x, f8.3, 1x, a)') rave, rerr, rxc, ryc, 
     $        z, ut, wave, anew, filename(ibeg:iend)
 

c  optionally, plot fit

      if (plot) then
         do i=1,NRING
            isam(i,1) = 0
            numb(i,1) = 0
         end do

         rmax = arad - sqrt((axc - rxc)**2 + (ayc - ryc)**2)
         rmaxsq = rmax * rmax
         iylo = max(nint(ryc - rmax), 1)
         iyhi = min(nint(ryc + rmax), naxes(2))

         do iy=iylo,iyhi
            ier=0
            firstpix= 1 + (iy-1)*(npixels)
            call ftgpve(inunit,group,firstpix,npixels,nullvale,
     &           dat,anynull,ier)
c            call imgl2r (img, dat, iy, ier)
            if (ier .ne. 0) then
               write (*, '(a, "Error while reading ring file")') bell
               write (10, '("Error while reading ring file")')
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

            if (numb(i,1) .gt. 1.0) dy(i) = dy(i) / 
     $         sqrt(numb(i,1) - 1.0)
            rry(i) = y(i)
            ryh(i) = y(i) + dy(i)
            ryl(i) = y(i) - dy(i)
         end do

         if (filter) call lowpass (y, NRING, icut, iwide)



         ymin = y(1)
         ymax = y(1)
         do i=1,NRING
            ymin = min(ymin, y(i), ryl(i))
            ymax = max(ymax, y(i), ryh(i))
         end do
         
         xv(1) = ((NRING * (rave / rmax)**2) + 0.5) / NRING
         xv(2) = xv(1)
         yv(1) = ymin
         yv(2) = ymax
         ymin = 0.90 * ymin
         ymax = 1.10 * ymax
         call pgenv (0.0, 1.0, ymin, ymax, 0, 1)
         call pglabel ("radius\u2\d", "intensity", 
     $        filename(ibeg:iend))
         call pgpoint (NRING, x, rry, 17)
         call pgerry (NRING, x, ryh, ryl, 1.0)
         call pgsci (7)
         if (filter) call pgline (NRING, x, y)
         call pgsci (2)
         call pgline (2, xv, yv)
         call pgsci (1)
      end if




c  ring completed - loop for more

 300  if (imgopen) then
         call ftclos (inunit, ier)
         call ftfiou (inunit, ier)
         imgopen = .false.
      endif

      if (list) go to 200

c      if (nargs .lt. 1) then
c         ier = 1
c         do while (ier .ne. 0)
c            write (*, '("Do another ring or list (y/n)? "$)')
c            read (*, '(a)') param
c            call yesno (param, query, ier)
c         end do
c         if (query) go to 100
c      end if


c  done - close up shop

      close (11)
      write (*, '("Processing completed")')
      write (10, '("Processing completed")')
      close (10)
      if (plot) call pgend ()

c      stop


c  fatal error - close files and terminate

 666  inquire (unit=10, opened=query)
      if (query) then
         write (10, '("Execution terminating")')
         close (10)
      end if

      inquire (unit=11, opened=query)
      if (query) close (11)

      inquire (unit=12, opened=query)
      if (query) close (12)

      if (imgopen) then
         call ftclos (inunit, ier)
         call ftfiou (inunit, ier)
      endif

      if (plot) call pgend ()

      
      return
      end subroutine nightring


