c------------------------------------------------------------------------------
      program calring_MR
c------------------------------------------------------------------------------

c  This program measures the center coordinates and radius of one or more
c  PFIS FP calibration rings.


c	use DFPORT
c	use DFLIB

      parameter (NRING=512, MAXIN=17622192)

	real*8 sec
      real ring(MAXIN)
      real x(NRING), y(NRING+2), dy(NRING)
      real xv(2), yv(2)
 
      integer isize(2)

      logical  query, plot, list, filter, fixed
      character*80 filename, label, param, comment
      character*24 datetime

	common /pfppar/ axc, ayc, arad, rxcen, rycen, icut, iwide, itmax,
     $                converg, wscale, rlo, rhi


c  initialization

      iargs = iargc ()
      write (*, '(/"Calibration Ring Fitting"//)')

	if (iargs .eq. 0) then
	   ier = 1
	   do while (ier .ne. 0)
	      write (*,'("Enter directory path name: "$)')
	      read (*,'(a)') filename
	      if (lnblnk(filename) .gt. 0) then
	         ier = chdir(filename)
	         if (ier .ne. 0) write (*, '("Error changing directory")')
	      else
	         ier = 0
	      end if
	   end do
	end if

c  start the log file

	call getlu (lulog)
      open (lulog, file="pfp.log", access="append", iostat=ier)
      if (ier .ne. 0) then
         write (*, '("Error opening log file")')
         stop
      end if
      write (lulog, '(/"Calibration Ring Fitting")')
      call fdate (datetime)
      write (lulog, '(a24)') datetime


c  open the output file

	call getlu (luout)
      if (iargs .ge. 2) then
         call getarg (2, filename)
      else
         filename = " "
         do while (lnblnk(filename) .le. 0)
            write (*, '("Enter a filename for the output: "$)')
            read (*, '(a)') filename
         end do
      end if
      call strlim (filename, ibeg, iend)
      inquire (file=filename, exist=query)
      open (luout, file=filename, access="append", iostat=ier)
      if (ier .ne. 0) then
         write (*, '("Error opening output file ",a)')
     $        filename(ibeg:iend)
         write (lulog, '("Error opening output file ",a)')
     $        filename(ibeg:iend)
         close (lulog)
	   stop
      end if
      write (lulog, '("Output file: ",a)') filename(ibeg:iend)
      if (.not. query) then
         if (iargs .ge. 2) then
            label = datetime
         else
            write (*, '("Enter a comment for the output file:")')
            read (*, '(a)') label
            if (lnblnk(label) .le. 0) then
               label = datetime
            end if
         end if
         write (luout, '(a)') label(1:lnblnk(label))
         write (luout, '("  radius   err     xc      yc      z    ",
     $        "  ut     wave   dn  file")')
      end if


c  get the options for the ring fitting operation

      call getpfp ("axc", param)
      read (param, *, iostat=ier) axc
      do while (ier .ne. 0)
         write (*, '("Enter aperture x center: "$)')
         read (*, *, iostat=ier) axc
      end do

      call getpfp ("ayc", param)
      read (param, *, iostat=ier) ayc
      do while (ier .ne. 0)
         write (*, '("Enter aperture y center: "$)')
         read (*, *, iostat=ier) ayc
      end do

      call getpfp ("arad", param)
      read (param, *, iostat=ier) arad
      do while (ier .ne. 0)
         write (*, '("Enter aperture radius: "$)')
         read (*, *, iostat=ier) arad
      end do

      call getpfp ("rxc", param)
      read (param, *, iostat=ier) rxcen
      do while (ier .ne. 0)
         write (*, '("Enter ring x center: "$)')
         read (*, *, iostat=ier) rxcen
      end do

      call getpfp ("ryc", param)
      read (param, *, iostat=ier) rycen
      do while (ier .ne. 0)
         write (*, '("Enter ring y center: "$)')
         read (*, *, iostat=ier) rycen
      end do

      call getpfp ("calring_filter", param)
      call yesno (param, filter, ier)
      do while (ier .ne. 0)
         write (*, '("Filter the profiles (y/n)? "$)')
         read (*, '(a)') param
         call yesno (param, filter, ier)
      end do

      if (filter) then
         call getpfp ("calring_filter_freq", param)
         read (param, *, iostat=ier) icut
         do while (ier .ne. 0)
            write (*, '("Enter filter frequency: "$)')
            read (*, *, iostat=ier) icut
         end do

         call getpfp ("calring_filter_width", param)
         read (param, *, iostat=ier) iwide
         do while (ier .ne. 0)
            write (*, '("Enter filter cutoff width: "$)')
            read (*, *, iostat=ier) iwide
         end do
      end if

      call getpfp ("calring_plot", param)
      call yesno (param, plot, ier)
      do while (ier .ne. 0)
         write (*, '("Plot the fits (y/n)? "$)')
         read (*, '(a)') param
         call yesno (param, plot, ier)
      end do
      if (plot) then
         call pgbegin (0,"?",1,2)
         call pgask (.false.)
      end if

	call getpfp ("calring_itmax", param)
	read (param, *, iostat=ier) itmax
	if (ier .ne. 0) itmax = 50

	call getpfp ("calring_conv", param)
	read (param, *, iostat=ier) converg
	if (ier .ne. 0) converg = 0.01

	call getpfp ("calring_fitwidth", param)
	read (param, *, iostat=ier) wscale
	if (ier .ne. 0) wscale = 0.5

	call getpfp ("calring_rlo", param)
	read (param, *, iostat=ier) rlo
	if (ier .ne. 0) rlo = -1.0

	call getpfp ("calring_rhi", param)
	read (param, *, iostat=ier) rhi
	if (ier .ne. 0) rhi = -1.0

      call getpfp ("calring_fixed", param)
      call yesno (param, fixed, ier)
      do while (ier .ne. 0)
         write (*, '("Fix the ring center (y/n)? "$)')
         read (*, '(a)') param
         call yesno (param, fixed, ier)
      end do


c  get the input filename (or list)

 100  if (iargs .ge. 1) then
         call getarg (1, filename)
      else
         filename = " "
         do while (lnblnk(filename) .le. 0)
            write (*, '("Enter ring filename or list: "$)')
            read (*, '(a)') filename
         end do
      end if

      call strlim (filename, ibeg, iend)
      if (filename(ibeg:ibeg) .eq. "@") then
	   call getlu (lulist)
         open (lulist, file=filename(ibeg+1:iend), status="old", 
     $         iostat=ier)
         if (ier .ne. 0) then
            if (iargs .ge. 1) then
               write (lulog, '("Error opening image list ",a)')
     $              filename(ibeg+1:iend)
               write (lulog, '("Execution terminating")')
               close (lulog)
	         close (luout)
	         call pgend ()
	         stop
            else
               write (*, '("Error opening image list ",a)')
     $              filename(ibeg+1:iend)
               go to 100
            end if
         end if
         list = .true.
      else
         list = .false.
      end if


c  open next ring image

 200  if (list) then
         read (lulist, '(a)', iostat=ier) filename
         if (ier .ne. 0) then
            if (ier .gt. 0) then
               write (*, '("Error reading ring filename from list")')
               write (lulog,'("Error reading ring filename from list")')
            end if
            close (lulist)
            list = .false.
            go to 300
         end if
         call strlim (filename, ibeg, iend)
         write (*, '("Processing ring: ",a)') filename(ibeg:iend)
      end if
      write (lulog, '("Processing ring: ",a)') filename(ibeg:iend)

	call getlu (luimg)
	ier = 0
      call ftopen (luimg, filename, 1, iblk, ier)
      if (ier .ne. 0) then
         write (*, '("Error opening ring image")')
         write (lulog, '("Error opening ring image")')
         go to 300
      end if

      call ftgipr (luimg, 2, ibitpix, naxis, isize, ier)
      if (ier .ne. 0) then
         write (*, '("Error determining image size")')
         write (lulog, '("Error determining image size")')
	   call ftclos (luimg, ier)
         go to 300
      end if

      if (naxis .ne. 2) then
         write (*, '("Error: wrong image dimension: ",i4)') naxis
         write (lulog, '("Error: wrong image dimension: ",i4)') naxis
	   call ftclos (luimg, ier)
         go to 300
      end if

	npix = isize(1) * isize(2)
        print*, isize(1), isize(2), npix, ' naxis 1, naxis 2, no of pix'

      if (npix .gt. MAXIN) then
         write (*, '("Ring image too large: ",i9)') npix
         write (lulog, '("Ring image too large: ",i9)') npix
	   call ftclos (luimg, ier)
         go to 300
      end if

	call ftgkys (luimg, "UTC-OBS", param, comment, ier)
	if (ier .eq. 0) then
	   call fts2tm (param,iyr,imo,ida,ihr,imn,sec,ier)
	end if
      do while (ier .ne. 0)
	   param = " "
	   do while (lnblnk(param) .le. 0)
	      write (*,'("Enter UT of ring (hh:mm:ss): ")')
	      read (*,'(a)') param
	   end do
	   ier = 0
	   call fts2tm (param,iyr,imo,ida,ihr,imn,sec,ier)
	end do
	ut = ihr + imn/60.0 + sec/3600.0

c NOTE HARDWIRED FOR ETALON 1 AT THE MOMENT

      call ftgkyj (luimg, "ET1Z", izval, comment, ier)
      if (ier .ne. 0) then
         write (*, '("Enter Z value for this ring: ")')
         read (*,*) izval
      end if
      z = izval / 1000.0

      call ftgkye (luimg, "WAVE", wave, comment, ier)
      do while (ier .ne. 0)
         write (*, '("Enter wavelength (", f8.3, "): ")') wave
         read (*, '(a)') param
         if (lnblnk(param) .gt. 0) then
            read (param, *, iostat=ier) wave
         else
            ier = 0
         end if
      end do

	call ftgpve (luimg, 0, 1, npix, 0, ring, ext, ier)
	if (ier .ne. 0) then
	   write (*,'("Error reading ring image")')
	   write (lulog,'("Error reading ring image")')
	   call ftclos (luimg, ier)
	   go to 300
	end if
	call ftclos (luimg, ier)

c  display and query limits

      if (plot) then
         call ringpro (ring, isize(1), isize(2), .false., x, y, dy)
	   xmin = x(1)
	   xmax = x(1)
         ymin = y(1)
         ymax = y(1)
         do i=2,NRING
            xmin = min(xmin, x(i))
            xmax = max(xmax, x(i))
            ymin = min(ymin, y(i))
            ymax = max(ymax, y(i))
         end do
         xv(1) = rave
         xv(2) = xv(1)
         yv(1) = ymin
         yv(2) = ymax
	   xmin = xmin - 0.05 * (xmax - xmin)
	   xmin = max(xmin, 0.0)
	   xmax = xmax + 0.05 * (xmax - xmin)
	   ymin = ymin - 0.05 * (ymax - ymin)
	   ymax = ymax + 0.05 * (ymax - ymin)
         call pgenv (xmin, xmax, ymin, ymax, 0, 1)
         call pglabel ("radius", "intensity", filename(ibeg:iend))
         call pgline (NRING, x, y)
	   write (*,'("Enter min fit r (0): ")')
         read (*, '(a)') param
         if (lnblnk(param) .gt. 0) then
            read (param, *, iostat=ier) rlo
         else
            rlo = -1.0
         end if
	   write (*,'("Enter max fit r (all): ")')
         read (*, '(a)') param
         if (lnblnk(param) .gt. 0) then
            read (param, *, iostat=ier) rhi
         else
            rhi = -1.0
         end if
	   if (rlo .gt. 0.0) xmin = rlo
	   if (rhi .gt. 0.0) xmax = rhi
         call pgenv (xmin, xmax, ymin, ymax, 0, 1)
         call pglabel ("radius", "intensity", filename(ibeg:iend))
         call pgline (NRING, x, y)
	end if


c  fit the ring

	call ringcen (ring,isize(1),isize(2),fixed,rave,rerr,rxc,ryc)
	if (rave .lt. 0.0) then
	   ier = nint(rave)
	   write (*,'("Error fitting ring: ",i3)') ier
	   write (lulog,'("Error fitting ring: ",i3)') ier
	   go to 250
	end if

	rxsav = rxcen
	rysav = rycen
	rxcen = rxc
	rycen = ryc

c  write output data

      write (*, '("Radius: ",f7.2, " +/- ", f6.2, 5x, "X center: ",
     $     f7.2, 5x, "Y center: ", f7.2)') rave, rerr, rxc, ryc
      write (luout, '(1x, f7.2, 1x, f6.2, 2(1x, f7.2), 1x, f6.3, 1x, 
     $     f7.4, 1x, f8.3, 1x, " 0", 1x, a)') rave, rerr, rxc, ryc, 
     $     z, ut, wave, filename(ibeg:iend)


c  optionally, plot fit

250   if (plot) then
         call ringpro (ring, isize(1), isize(2), .false., x, y, dy)
	   xmin = x(1)
	   xmax = x(1)
         ymin = y(1)
         ymax = y(1)
         do i=2,NRING
            xmin = min(xmin, x(i))
            xmax = max(xmax, x(i))
            ymin = min(ymin, y(i))
            ymax = max(ymax, y(i))
         end do
         xv(1) = rave
         xv(2) = xv(1)
         yv(1) = ymin
         yv(2) = ymax
	   xmin = xmin - 0.05 * (xmax - xmin)
	   xmin = max(xmin, 0.0)
	   xmax = xmax + 0.05 * (xmax - xmin)
	   ymin = ymin - 0.05 * (ymax - ymin)
	   ymax = ymax + 0.05 * (ymax - ymin)
         call pgenv (xmin, xmax, ymin, ymax, 0, 1)
         call pglabel ("radius", "intensity", filename(ibeg:iend))
         call pgline (NRING, x, y)
         call pgsci (2)
         if (rave .gt. 0.0) call pgline (2, xv, yv)
         call pgsci (1)

	   scale = xmax * xmax
	   do i=1,NRING
	      x(i) = x(i) * x(i) / scale
	   end do
	   xv(1) = xv(1) * xv(1) / scale
	   xv(2) = xv(1)
	   call pgenv (0.0, 1.0, ymin, ymax, 0, 1)
         call pglabel ("radius\u2\d", "intensity", filename(ibeg:iend))
         call pgline (NRING, x, y)
         call pgsci (2)
         if (rave .gt. 0.0) call pgline (2, xv, yv)
         call pgsci (1)
      end if

	rxcen = rxsav
	rycen = rysav

c  ring completed - loop for more

 300  continue

      if (list) go to 200

      if (iargs .lt. 1) then
         ier = 1
         do while (ier .ne. 0)
            write (*, '("Do another ring or list (y/n)? ")')
            read (*, '(a)') param
            call yesno (param, query, ier)
         end do
         if (query) go to 100
      end if


c  done - close up shop

      close (luout)
      write (*, '("Processing completed")')
      write (lulog, '("Processing completed")')
      close (lulog)
      if (plot) call pgend ()

      stop
      end


