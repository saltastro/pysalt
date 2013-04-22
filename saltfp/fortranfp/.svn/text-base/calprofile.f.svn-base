c------------------------------------------------------------------------------
      subroutine calprofile(dir,calprofilelogfile,outfile,label,
     >     axc,ayc,arad,rxcen,rycen,filter,icut,iwidth, plottype,
     >     itmax,converg,wscale,rlo,rhi,rfixed,fixed,filename, 
     >     cala, calb, calc, cald, calf)
c------------------------------------------------------------------------------

c  This program fits a Voigt profile to one or more calibration rings,
c  displays the fit, and records the fit parameters for later use.
c 13.05.10 Updated to read which etalon is in use and take correct Z from
c fits header. (nsl).
c 23.06.10 Updated to parse pfp parameters from python, removing common block.(nsl)
c 25.06.10 Updated to allow a fixed rlo and rhi, so that user isnt prompted
c each time. (nsl)

c	use DFPORT
c	use DFLIB

 
      parameter (NRING=512, MAXIN=17622192)
 
      real ring(MAXIN)
      real x(NRING), y(NRING+2), dy(NRING)
      real a(5), dvda(5), da(5), calfit(5), xv(2), yv(2)
      real cala, calb, calc, cald, calf
      integer isize(2), naxis
      logical flag(5), query, filter, list, fixed, rfixed
      character*80 filename,label,param, comment,calprofilelogfile,dir,
     >     outfile
      character*24 datetime
      character*4 etalonz, plottype
      character*20 etstate

 
c      common /pfppar/ axc, ayc, arad, rxcen, rycen, icut, iwide, itmax,
c     $                converg, wscale, rlo, rhi


c      write(*, '(a)') dir
c      write(*, '(a)') calprofilelogfile
c      write(*, '(a)') outfile
c      write(*, '(a)') label
c      print*, axc,ayc,arad,rxcen,rycen
c      print*, filter,icut,iwidth
c      write(*, '(a)') plottype
c      print*, itmax,converg,wscale,rlo,rhi,rfixed,fixed
c      write(*, '(a)') filename 
c      print*, cala, calb, calc, cald, calf

c  initialization

      do i=1,5
         flag(i) = .true.
      end do
c      iargs = iargc ()
      write (*, '(/"Calibration Ring Profile "//)')

	ier = 1
c	do while (ier .ne. 0)
c	   write (*,'("Enter directory path name: "$)')
c	   read (*,'(a)') dir
	   if (lnblnk(dir) .gt. 0) then
	      ier = chdir(dir)
	      if (ier .ne. 0) write (*, '("Error changing directory")')
	   else
	      ier = 0
	   end if
c	end do
	
c  start the log file

	call getlu (lulog)
      open (lulog, file=calprofilelogfile, access="append", iostat=ier)
      if (ier .ne. 0) then
         write (*, '("Error opening log file")')
         stop
      end if
      write (lulog, '(/"Calibration Ring Profile")')
      call fdate (datetime)
      write (lulog, '(a24)') datetime


c  open the output file

      call getlu (luout)
c      if (iargs .ge. 2) then
c         call getarg (2, outfile)
c      else
c         filename = " "
c         do while (lnblnk(outfile) .le. 0)
c            write (*, '("Enter a filename for the output: "$)')
c            read (*, '(a)') outfile
c         end do
c      end if

      inquire (file=outfile, exist=query)
      open (luout, file=outfile, access="append", iostat=ier)
      if (ier .ne. 0) then
         write (*, '("Error opening output file ",a)')
     $        outfile(1:lnblnk(outfile))
         write (lulog, '("Error opening output file ",a)')
     $        outfile(1:lnblnk(outfile))
         write (lulog, '("Execution terminating")')
	   close (lulog)
	   stop       
      end if
      if (.not. query) then
c         if (iargs .ge. 2) then
c            label = datetime
c         else
c            write (*, '("Enter a comment for the output file: ")')
c            read (*, '(a)') label
            if (lnblnk(label) .le. 0) label = datetime
c         end if
         write (luout, '(a)') label(1:lnblnk(label))
      end if
c comment will only be added to a new file.

! note to ted - compress this output format?
      write (luout, 50)
 50   format (4x, "Continuum", 10x, "Line", 13x, "Center", 8x, 
     $     "Gaussian",  5x, "Lorentzian")
      write (luout, 51)
 51   format(4x, "Peak",7x,"FWHM",3x,"R. Centre X", 3x,"R. Centre Y")

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
c         write (*, '("Enter ring X center: "$)')
c         read (*, *, iostat=ier) rxcen
c      end do

c      call getpfp ("ryc", param)
c      read (param, *, iostat=ier) rycen
c      do while (ier .ne. 0)
c         write (*, '("Enter ring Y center: "$)')
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
c            write (*,'("Enter filter frequency: "$)')
c            read (*, *, iostat=ier) icut
c         end do

c         call getpfp ("calring_filter_width", param)
c         read (param, *, iostat=ier) iwide
c         do while (ier .ne. 0)
c            write (*, '("Enter filter cutoff width: "$)')
c            read (*, *, iostat=ier) iwide
c         end do
c	else
c	   icut = 0
c	   iwide = 0
c      end if

c	call getpfp ("calring_itmax", param)
c	read (param, *, iostat=ier) itmax
c	if (ier .ne. 0) itmax = 50
c
c	call getpfp ("calring_conv", param)
c	read (param, *, iostat=ier) converg
c	if (ier .ne. 0) converg = 0.01
c
c	call getpfp ("calring_fitwidth", param)
c	read (param, *, iostat=ier) wscale
c	if (ier .ne. 0) wscale = 0.5
c
c	call getpfp ("calring_rlo", param)
c	read (param, *, iostat=ier) rlo
c	if (ier .ne. 0) rlo = -1.0
c
c	call getpfp ("calring_rhi", param)
c	read (param, *, iostat=ier) rhi
c	if (ier .ne. 0) rhi = -1.0

c      call getpfp ("calibration_a", param)
c      read (param, *, iostat=ier) calfit(1)
c      do while (ier .ne. 0)
c         write (*, '("Enter calibration A value: "$)')
c         read (*, *, iostat=ier) calfit(1)
c      end do
      calfit(1) = cala

c      call getpfp ("calibration_b", param)
c      read (param, *, iostat=ier) calfit(2)
c      do while (ier .ne. 0)
c         write (*, '("Enter calibration B value: "$)')
c         read (*, *, iostat=ier) calfit(2)
c      end do
      calfit(2) = calb


c      call getpfp ("calibration_c", param)
c      read (param, *, iostat=ier) calfit(3)
c      do while (ier .ne. 0)
c         write (*, '("Enter calibration C value: "$)')
c         read (*, *, iostat=ier) calfit(3)
c      end do
      calfit(3) = calc


c      call getpfp ("calibration_d", param)
c      read (param, *, iostat=ier) calfit(4)
c      do while (ier .ne. 0)
c         write (*, '("Enter calibration D value: "$)')
c         read (*, *, iostat=ier) calfit(4)
c      end do
      calfit(4) = cald


c      call getpfp ("calibration_f", param)
c      read (param, *, iostat=ier) calfit(5)
c      do while (ier .ne. 0)
c         write (*, '("Enter calibration F value: "$)')
c         read (*, *, iostat=ier) calfit(5)
c      end do
      calfit(5) = calf

c      call getpfp ("calring_fixed", param)
c      call yesno (param, fixed, ier)
c      do while (ier .ne. 0)
c         write (*, '("Fix the ring centers (y/n)? "$)')
c         read (*, '(a)') param
c         call yesno (param, fixed, ier)
c      end do



c      call pgbegin (0, "?", 1, 1)
      call pgbegin(0, plottype, 1, 1)
      call pgask (.false.)


c  get the input filename (or list)

c 100  if (iargs .ge. 1) then
c         call getarg (1, filename)
c      else
c         filename = " "
c         do while (lnblnk(filename) .le. 0)
c            write (*, '("Enter ring filename or list:  "$)')
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
               close (luout)
	         close (lulog)
	         call pgend ()
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


c  ring processing loop

 200  if (list) then
         read (lulist, '(a)', iostat=ier) filename
         if (ier .ne. 0) then
            if (ier .gt. 0) then
               write (*, '(a)') "Error reading ring filename from list"
               write (lulog,'(a)') 
     >              "Error reading ring filename from list"
            end if
            close (lulist)
            list = .false.
            go to 300
         end if
         call strlim (filename, ibeg, iend)
         write (*, 20) "Processing ring: ", filename(ibeg:iend)
      end if
 20   format(a20, a20)
	call getlu (luimg)
	ier = 0
        write(*,'(a)') filename
	call ftopen (luimg, filename, 1, iblk, ier)
	if (ier .ne. 0) then
		write (*,15) "Error opening ring image, status = ", ier
		write (lulog,15) "Error opening ring image, status = ",ier
		go to 300
	end if
 15     format(a45,i1)

        naxis=0
        istat=0
	call ftgipr (luimg, 2, ibitpix, naxis, isize, istat)
c        print*, ibitpix, naxis, isize, istat
      if (istat .ne. 0) then
         write (*, '("Error determining image size")')
         write (lulog, '("Error determining image size")')
	   call ftclos (luimg, ier)
         go to 300
      end if
c      print*, naxis, ' naxis'     

 25   format(a40, i4)
      if (naxis .ne. 2) then
         write (*,25) "Error: wrong image dimension: ", naxis
         write (lulog, 25) "Error: wrong image dimension: ",naxis
	   call ftclos (luimg, ier)
         go to 300
      end if
	npix = isize(1) * isize(2)
	if (npix .gt. MAXIN) then
         write (*,15) "Ring image too large: ", isize(1)
         write (lulog,15) "Ring image too large: ",isize(1)
	   call ftclos (luimg, ier)
         go to 300
      end if

C     Updated to read which etalon is in use if only one and to read from header.

      call FTGKYS(luimg,"ET-STATE", etstate,comment,ier)
      write(*,'(a)') etstate(6:7)
      if (etstate(1:2) .eq. 'S2') then !Et1
         etalonz='ET1Z'
      else if (etstate(1:2) .eq. 'S3') then !Et2
         etalonz='ET2Z'
      else
         print*, 'problem with which etalon to choose'
         etalonz='ETZ' ! this doesnt exist, forces user input
      endif
      ier=0
      call ftgkyj (luimg, etalonz, izval, comment, ier)
      do while (ier .ne. 0)
         write (*, '("Enter Z value: "$)')
         read (*, *, iostat=ier) izval
      end do
      z = izval / 1000.0
      
      call ftgpve (luimg, 0, 1, npix, 0, ring, query, ier)
      call ftclos (luimg, ier)
 
c  plot profile and get limits

      call ringpro (ring,isize(1),isize(2),.false.,x,y,dy,
     >     axc,ayc,arad,rxcen,rycen,icut,iwide,itmax,
     >     converg,wscale,rlo,rhi)
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
      yscale = ymax - ymin
      ymin = ymin - 0.05 * yscale
      ymax = ymax + 0.05 * yscale
      xscale = xmax - xmin
      xmin = max(0.0, xmin - 0.05 * xscale)
      xmax = xmax + 0.05 * xscale

      call pgenv (xmin, xmax, ymin, ymax, 0, 1)
      call pglabel ("radius", "intensity", filename(ibeg:iend))
      call pgpoint (NRING, x, y, 0)

      if (.not. rfixed) then
	write (*,'("Enter lower and upper radial limits: "$)')
	read (*,'(a)') param
	ier = 1
	read (param, *, iostat=ier) rlo, rhi
        if (rhi .le. rlo) then
           write(*, '("Error: rlow is greater than rhigh")')
           ier = 1
        endif
	do while (ier .ne. 0) 
           write(*, '("Please re-enter lower and upper radial limits")')
           read (*,'(a)') param
           ier = 1
           read (param, *, iostat=ier) rlo, rhi
           if (rhi .le. rlo) then
              write(*,'("Error:rlo is greater than or equal to rhi")')
              ier = 1
           endif
c     rlo = -1.0
c     rhi = -1.0
        enddo        
      endif

c  measure ring center 
      call ringcen (ring,isize(1),isize(2),fixed,rave,rerr,rxc,ryc,
     >     axc,ayc,arad,rxcen,rycen,icut,iwide,itmax,
     >     converg,wscale,rlo,rhi)

	if (rave .gt. 0.0) then
	   rxcen = rxc
	   rycen = ryc
	end if
	wcen = wave(rave, z, calfit)

	call ringpro (ring,isize(1),isize(2),.false.,x,y,dy,
     >       axc,ayc,arad,rxcen,rycen,icut,iwide,itmax,
     >       converg,wscale,rlo,rhi)


c  plot the profile

	x(1) = wave(x(1), z, calfit)
	xmin = x(1)
	xmax = x(1)
        ymin = y(1)
        ymax = y(1)
        do i=2,NRING
	   x(i) = wave(x(i), z, calfit)
           xmin = min(xmin, x(i))
           xmax = max(xmax, x(i))	    
           ymin = min(ymin, y(i))
           ymax = max(ymax, y(i))
        end do
        yscale = ymax - ymin
        ymin = ymin - 0.05 * yscale
        ymax = ymax + 0.05 * yscale
        xmin = xmin - 0.5
        xmax = xmax + 0.5
        
      call pgenv (xmin, xmax, ymin, ymax, 0, 1)
      call pglabel ("wavelength", "intensity", filename(ibeg:iend))
      call pgpoint (NRING, x, y, 0)

	xv(1) = wcen
	xv(2) = wcen
	yv(1) = ymin
	yv(2) = ymax
	call pgsci (3)
	call pgline (2, xv, yv)
	call pgsci (1)

c  fit the profile and plot fit

      rmax = arad - sqrt ((axc - rxc)**2 + (ayc - ryc)**2)
      rmaxsq = rmax * rmax
	if (rlo .gt. 0.0) then
	   ilo = max(1,nint(NRING*rlo*rlo/rmaxsq))
	else
	   ilo = 1
	end if
	if (rhi .gt. 0.0) then
	   ihi = min(NRING,nint(NRING*rhi*rhi/rmaxsq))
	else
	   ihi = NRING
	end if
	num = ihi - ilo + 1

      call evinit (x(ilo), y(ilo), num, a)
	write (*,*) a

      call evfit (x(ilo), y(ilo), dy(ilo), num, a, flag, da, chisq)
      call evstat (a, ycen, fwhm)
      do i=ilo,ihi
         call evoigt (x(i), a, 5, y(i), dvda)
      end do
      call pgsci (2)
      call pgline (num, x(ilo), y(ilo))
      call pgsci (1)

	wlo = wcen - (fwhm/2.0)
	whi = wcen + (fwhm/2.0)
	x1 = wfind(x,NRING,wlo)
	x2 = wfind(x,NRING,wcen)
	x3 = wfind(x,NRING,whi)
	call integrate (y,NRING,x3,x2,area1)
	call integrate (y,NRING,x2,x1,area2)
	asym1 = (area1 - area2) / (area1 + area2)
	asym2 = wcen - a(3)

      write (luout, '(a)') filename
      write (luout, 210) a(1), da(1), a(2), da(2), a(3), da(3),
     $     a(4), da(4), a(5), da(5), ycen, fwhm,rxcen,rycen,asym1,asym2
      write (*, 210) a(1), da(1), a(2), da(2), a(3), da(3),
     $     a(4), da(4), a(5), da(5),ycen,fwhm,rxcen,rycen,asym1,asym2
 210  format (1x, f7.1, 1x, f7.1, 2x, f9.1, 1x, f8.1, 2x, f8.3, 1x, 
     $     f6.3, 2(2x, f6.2, 1x, f6.2), 1x, f7.1, 1x, f5.2, 
     $     4x, f7.2, 4x, f7.2,
     $     1x, f6.3, 1x, f6.3)

 300  continue
      if (list) go to 200

c      if (iargs .lt. 1) then
c         ier = 1
c         do while (ier .ne. 0)
c            write (*, '("Do another ring or list (y/n)? "$)')
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
      call pgend ()

      end subroutine calprofile


      function wave (r, z, calfit)

      real wave, r, z, calfit(5)

      wave = (calfit(1) + z *(calfit(2) + z*(calfit(3) + z* calfit(4))))
     $       / sqrt(1.0 + (r / calfit(5))**2)

      return
      end


	function wfind (x, n, w)

	real x(n)

	if ((w .lt. x(n)) .or. (w .gt. x(1))) then
	   wfind = -1.0
	   return
	end if

	loc = n
	do while (x(loc) .ge. w)
	   loc = loc - 1
	end do

	wfind = loc - ((w - x(loc))/(x(loc-1) - x(loc)))

	return
	end 
