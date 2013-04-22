c------------------------------------------------------------------------------
      subroutine calring(dir,calringlog, output, comment,
     >     axc,ayc,arad,rxc,ryc,filter, icut, iwide,plot,plottype,
     >     itmax,converg,wscale, rlo, rhi, rfixed,fixed,filename)
c------------------------------------------------------------------------------

c  This program measures the center coordinates and radius of one or more
c  PFIS FP calibration rings.
c 13.05.10 Updated to read which etalon is in use and take correct Z and wave value from
c fits header. (nsl).
c 22.06.10 Updated to parse pfp parameters from python, removing common block.(nsl)
c 25.06.10 Updated to allow a fixed rlo and rhi, so that user isnt prompted
c each time. (nsl)
c 05.05.11 Updated due to change in RSS keywords (nsl)

c	use DFPORT
c	use DFLIB

      parameter (NRING=512, MAXIN=17622192)

      real*8 sec
      real ring(MAXIN)
      real x(NRING), y(NRING+2), dy(NRING)
      real xv(2), yv(2)
 
      integer isize(2)

      logical  query, plot, list, filter, fixed, rfixed
      character*80 filename, label, param, comment,filename2,output,dir
      character*24 datetime, calringlog
      character*4 etalonz, plottype
      character*8 etalonl
      character*20 etstate

c      common /pfppar/ axc, ayc, arad, rxcen, rycen, icut, iwide, itmax,
c     $                converg, wscale, rlo, rhi


c  initialization

c      iargs = iargc ()
      write (*, '(/"Calibration Ring Fitting"//)')

c     if (iargs .eq. 0) then
      ier = 1
c     do while (ier .ne. 0)
c     write (*,'("Enter directory path name: "$)')
c     read (*,'(a)') filename2
      if (lnblnk(dir) .gt. 0) then
         ier = chdir(dir)
         if (ier .ne. 0) write (*, '("Error changing directory")')
      else
         ier = 0
      end if
c     end do
c     end if

c  start the log file

      call getlu (lulog)
      open (lulog, file=calringlog, access="append", iostat=ier)
      if (ier .ne. 0) then
         write (*, '("Error opening calring log file")')
         stop
      end if
      write (lulog, '(/"Calibration Ring Fitting")')
      call fdate (datetime)
      write (lulog, '(a24)') datetime


c  open the output file

	call getlu (luout)
c      if (iargs .ge. 2) then
c         call getarg (2, filename)
c      else
c         filename = " "
c         do while (lnblnk(filename) .le. 0)
c            write (*, '("Enter a filename for the output: "$)')
c            read (*, '(a)') filename
c         end do
c      end if
      call strlim (output, ibeg, iend)
      inquire (file=output, exist=query)
      open (luout, file=output, access="append", iostat=ier)
      if (ier .ne. 0) then
         write (*, '("Error opening output file ",a)')
     $        output(ibeg:iend)
         write (lulog, '("Error opening output file ",a)')
     $        output(ibeg:iend)
         close (lulog)
	   stop
      end if
      write (lulog, '("Output file: ",a)') output(ibeg:iend)
      if (.not. query) then
c         if (iargs .ge. 2) then
c            label = datetime
c         else
c            write (*, '("Enter a comment for the output file:")')
c            read (*, '(a)') label
            label = comment
            if (lnblnk(label) .le. 0) then
               label = datetime
            end if
c         end if
         write (luout, '(a)') label(1:lnblnk(label))
         write (luout, '("  radius   err     xc      yc      z    ",
     $        "  ut     wave   dn  file")')
      end if
c comment will only be added to a new file.

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
      if (plot) then
         call pgbegin (0,plottype,1,2)
c         call pgbegin (0,"?",1,2)
         call pgask (.false.)
      end if

c	call getpfp ("calring_itmax", param)
c	read (param, *, iostat=ier) itmax
c	if (ier .ne. 0) itmax = 50

c	call getpfp ("calring_conv", param)
c	read (param, *, iostat=ier) converg
c	if (ier .ne. 0) converg = 0.01

c	call getpfp ("calring_fitwidth", param)
c	read (param, *, iostat=ier) wscale
c	if (ier .ne. 0) wscale = 0.5

c	call getpfp ("calring_rlo", param)
c	read (param, *, iostat=ier) rlo
c	if (ier .ne. 0) rlo = -1.0

c	call getpfp ("calring_rhi", param)
c	read (param, *, iostat=ier) rhi
c	if (ier .ne. 0) rhi = -1.0

c      call getpfp ("calring_fixed", param)
c      call yesno (param, fixed, ier)
c      do while (ier .ne. 0)
c         write (*, '("Fix the ring center (y/n)? "$)')
c         read (*, '(a)') param
c         call yesno (param, fixed, ier)
c      end do


c  get the input filename (or list)

c 100  if (iargs .ge. 1) then
c         call getarg (1, filename)
c      else
c        filename = " "
c         do while (lnblnk(filename) .le. 0)
c            write (*, '("Enter ring filename or list: "$)')
c            read (*, '(a)') filename
c         end do
c      end if

      call strlim (filename, ibeg, iend)
      if (filename(ibeg:ibeg) .eq. "@") then
	   call getlu (lulist)
         open (lulist, file=filename(ibeg+1:iend), status="old", 
     $         iostat=ier)
         if (ier .ne. 0) then
c            if (iargs .ge. 1) then
               write (lulog, '("Error opening image list ",a)')
     $              filename(ibeg+1:iend)
               write (lulog, '("Execution terminating")')
               close (lulog)
	         close (luout)
	         call pgend ()
	         stop
c            else
c               write (*, '("Error opening image list ",a)')
c     $              filename(ibeg+1:iend)
c               go to 100
c            end if
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

	call ftgkys (luimg, "TIME-OBS", param, comment, ier)
	if (ier .eq. 0) then
	   call fts2tm (param,iyr,imo,ida,ihr,imn,sec,ier)
c           write(*, '(a)') param
c           print*, iyr, imo, ida, ihr, imn, sec, ier
	else
c for old data files
           ier=0
           call ftgkys (luimg, "UTC-OBS", param, comment, ier)   
           if (ier .eq. 0) then
              call fts2tm (param,iyr,imo,ida,ihr,imn,sec,ier)
           endif
        endif

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

C     Updated to read which etalon is in use if only one and to read from header.
C 5.5.11 Changed string extraction here as the fits headers have changed in new data. 

      call FTGKYS(luimg,"ET-STATE", etstate,comment,ier)
      write(*,'(a)') etstate(1:2)
      if ( (etstate(6:7) .eq. 'S2') .or. 
     >     (etstate(1:2) .eq. 'S2') ) then !Et1
         etalonz='ET1Z'
         etalonl='ET1WAVE0'
      else if ( (etstate(6:7) .eq. 'S3') .or. 
     >        (etstate(1:2) .eq. 'S3') ) then !Et2
         etalonz='ET2Z'
         etalonl='ET2WAVE0'
      else
         print*, 'problem with which etalon to choose'
         etalonz='ETZ' ! this doesnt exist, forces user input
         etalonl='ETWAVE0'  ! this doesnt exist, forces user input
      endif
      ier=0
      call ftgkyj (luimg, etalonz, izval, comment, ier)
      if (ier .ne. 0) then
         write (*, '("Enter Z value for this ring: "$)')
         read (*, *, iostat=ier) izval
      end if
      z = izval / 1000.0

      ier=1 ! set to false to force the user to enter the wavelength
c      call ftgkye (luimg, etalonl, wave, comment, ier)
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
         call ringpro (ring, isize(1), isize(2), .false., x, y, dy,
     >    axc,ayc,arad,rxc,ryc,icut,iwide,itmax,converg,wscale,rlo,rhi)
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
         
         if ((xmin .eq. xmax) .or. (ymin .eq. ymax)) then
            print*, 'plotting limit errors'
            go to 300
         endif
         


         call pgenv (xmin, xmax, ymin, ymax, 0, 1)
         call pglabel ("radius", "intensity", filename(ibeg:iend))
         call pgline (NRING, x, y)

c 25.06.10 - only do this if values not fixed in pyraf epar 
c allows user to use without inputting each time

         if (.not. rfixed) then
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
         endif


         if (rlo .gt. 0.0) xmin = rlo
         if (rhi .gt. 0.0) xmax = rhi
         call pgenv (xmin, xmax, ymin, ymax, 0, 1)
         call pglabel ("radius", "intensity", filename(ibeg:iend))
         call pgline (NRING, x, y)
	end if


c  fit the ring

	call ringcen (ring,isize(1),isize(2),fixed,rave,rerr,rxc,ryc,
     > axc,ayc,arad,rxc,ryc,icut,iwide,itmax,converg,wscale,rlo,rhi)
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
         call ringpro (ring, isize(1), isize(2), .false., x, y, dy,
     >     axc,ayc,arad,rxc,ryc,icut,iwide,itmax,converg,wscale,rlo,rhi)
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

c      if (iargs .lt. 1) then
c         ier = 1
c         do while (ier .ne. 0)
c            write (*, '("Do another ring or list (y/n)? ")')
c            read (*, '(a)') param
c            call yesno (param, query, ier)
c         end do
c         if (query) go to 100
c      end if


c  done - close up shop

      close (luout)
      write (*, '("Processing completed")')
      write (lulog, '("Processing completed")')
      close (lulog)
      if (plot) call pgend ()

      
      end subroutine calring


