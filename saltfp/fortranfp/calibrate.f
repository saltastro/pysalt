c------------------------------------------------------------------------------
      subroutine calibrate(plottype,filename, output, logfile)
c------------------------------------------------------------------------------

c   This program determines the wavelength calibration of the Fabry-Perot    
c   interactively, using the ring measurements from the CALRING program.     
c                                                                           



      parameter (MAXRINGS=100, NPLOT=200)

      real z(MAXRINGS), wave(MAXRINGS), r(MAXRINGS), dr(MAXRINGS)
      real t(MAXRINGS), dn(MAXRINGS), xc(MAXRINGS), yc(MAXRINGS)

      real x(4,MAXRINGS), y(MAXRINGS), dy(MAXRINGS)
      real fit(7), err(7), covar(7,7)
      real ff(7), dwdf(7)
      real xplot(NPLOT), yplot(NPLOT)
      integer idn, iplot, i, j, nargs
	integer pgbegin, pgcurse
      logical flags(7), iflag(MAXRINGS), multiorder

      character*80 filename, param, label, fname(MAXRINGS), output, 
     >     logfile
      character*24 datetime, xlabel, plottype
      character*1 answ, bell

      external lambda


c  initialization

      bell = char(7)
      iplot = 1
      do i=1,MAXRINGS
         iflag(i) = .true.
      end do
!      nargs = iargc ()

      write (*, '(/"Fabry-Perot Wavelength Calibration"//)')

c      if (pgbegin (0,"?",2,2) .ne. 1) then
      if (pgbegin (0,plottype,2,2) .ne. 1) then
		write (*,'("PGPlot error")')
		stop
	end if
      call pgask (.false.)


c  start the log file

      open (10, file=logfile, access="append", iostat=ier)
      if (ier .ne. 0) then
         write (*, '(a1,"Error opening log file")') bell
         stop
      end if
!      write (10, '(/"Wavelength Calibration")')
      call fdate (datetime)
      write (10, '(a24)') datetime

c  get the input data

!      if (nargs .gt. 0) then
!         call getarg (1, filename)
!         len = lnblnk(filename)
!      else
         len = 0
c         do while (len .eq. 0)
c            write (*, '("Enter input filename: "$)')
c            read (*, '(a)') filename
            len = lnblnk(filename)
c         end do
!      end if
      ier = 1
      do while (ier .ne. 0)
         open (20, file=filename, status="old", iostat=ier)
         if (ier .ne. 0) then
            write (*, '(a1,"Error opening file ",a)')
     $           bell, filename(1:len)
            len = 0
            do while (len .eq. 0)
               write (*, '("Enter input filename: "$)')
               read (*, '(a)') filename
               len = lnblnk(filename)
            end do
         end if
      end do
      write (10, '("Input file: ",a)') filename(1:len)

      read (20,'(a)') label
      read (20,'(a)') label

      multiorder = .false.
      nrings = 1
      do while (nrings .le. MAXRINGS)
         read (20, *, end=100)
     $        r(nrings), dr(nrings), xc(nrings), yc(nrings), 
     $        z(nrings), t(nrings), wave(nrings), idn, fname(nrings)
         dn(nrings) = idn
         if (idn .ne. 0) multiorder = .true.
         nrings = nrings + 1
      end do

 100  close (20)
      nrings = nrings - 1

      wmin = wave(1)
      wmax = wave(1)
      zmin = z(1)
      zmax = z(1)
      tmin = t(1)
      tmax = t(1)
      do i=1,nrings
         wmin = min(wmin,wave(i))
         wmax = max(wmax,wave(i))
         zmin = min(zmin,z(i))
         zmax = max(zmax,z(i))
         tmin = min(tmin,t(i))
         tmax = max(tmax,t(i))
      end do

      if ((tmax - tmin) .gt. 12.0) then
         do i=1,nrings
            if (t(i) .gt. 12.0) t(i) = 24.0 - t(i)
         end do
         tmin = t(1)
         tmax = t(1)
         do i=1,nrings
            tmin = min(tmin,t(i))
            tmax = max(tmax,t(i))
         end do
      end if

      tmid = (tmin + tmax) / 2.0
      do i=1,nrings
         t(i) = t(i) - tmid
      end do
      tmin = tmin - tmid
      tmax = tmax - tmid
      scale = tmax - tmin
      tlo = tmin - 0.05 * scale
      thi = tmax + 0.05 * scale

      scale = zmax - zmin
      zlo = zmin - 0.05 * scale
      zhi = zmax + 0.05 * scale

      rmin = 0.0
!      rmax = arad + sqrt((axc - rxc)**2 + (ayc - ryc)**2)
	rmax = 950
      scale = rmax - rmin
      rlo = rmin - 0.05 * scale
      rhi = rmax + 0.05 * scale

c  initial guesses for the fit

      c = 0.0
      d = 0.0
      e = 0.0
!      f = 24.0 * arad
      f = 11000.0
      fsq = f * f

      do i=1,nrings-1
         do j=i+1,nrings
            if ((wave(i).ne.wave(j)).and.(dn(i).eq.dn(j))) go to 200
         end do
      end do
      write (*, '(a,"Error - not enough rings")') bell
      go to 666
 200  wi = wave(i) * sqrt(1.0 + (r(i)**2 / fsq))
      wj = wave(j) * sqrt(1.0 + (r(j)**2 / fsq))
      b = (wj - wi) / (z(j) - z(i))
      a = wi - b * z(i)

      if (multiorder) then
         do j=1,nrings
            if (dn(j) .ne. dn(i)) go to 210
         end do
 210     temp = wave(j) * sqrt(1.0 + (r(j)**2 / fsq)) / (a + b*z(j))
         ord = (temp * dn(j) - dn(i)) / (1.0 - temp)
         a = a * (ord + dn(i))
         b = b * (ord + dn(i))
      else
         ord = 1.0
      end if

      fit(1) = a
      fit(2) = b
      fit(3) = c
      fit(4) = d
      fit(5) = f
      fit(6) = e
      fit(7) = ord

c start with linear fit

      flags(1) = .true.
      flags(2) = .true.
      flags(3) = .false.
      flags(4) = .false.
      flags(5) = .true.
      flags(6) = .false.
      if (multiorder) then
         flags(7) = .true.
      else
         flags(7) = .false.
      end if

c fitting loop

 300  continue

      nfit = 0
      do i=1,nrings
         if (iflag(i)) then
            nfit = nfit + 1
            x(1,nfit) = z(i)
            x(2,nfit) = r(i)
            x(3,nfit) = t(i)
            x(4,nfit) = dn(i)
            y(nfit) = wave(i)
            dy(nfit) = 1.0
         end if
      end do
      ndof = nfit
      do i=1,7
         if (flags(i)) ndof = ndof - 1
      end do


      call mrqfit2d (x, y, dy, 4, nfit, fit, flags, 7, covar, chisq,
     $             lambda)
      rms = sqrt(chisq / ndof)

      do i=1,7
         if (covar(i,i) .gt. 0.0) then
            err(i) = rms * sqrt(covar(i,i))
         else
            err(i) = 0.0
         end if
      end do

      a = fit(1) / fit(7)
      aerr = err(1) / fit(7)
      b = fit(2) / fit(7)
      berr = err(2) / fit(7)
      c = fit(3) / fit(7)
      cerr = err(3) / fit(7)
      d = fit(4) / fit(7)
      derr = err(4) / fit(7)
      e = fit(6) / fit(7)
      eerr = err(6) / fit(7)
      f = fit(5)
      ferr = err(5)

      write (*, '("Wave0 : ", f8.3, " +/- ", f6.3)') a, aerr
      write (*, '("dW/dZ : ", f8.4, " +/- ", f6.4)') b, berr
      write (*, '("dW/dZ2: ", f8.4, " +/- ", f6.4)') c, cerr
      write (*, '("dW/dZ3: ", f8.4, " +/- ", f6.4)') d, derr
      write (*, '("F     : ", f8.2, " +/- ", f7.2)') f, ferr
      write (*, '("dW/dT : ", f8.4, " +/- ", f8.4)') e, eerr
      write (*, '("n     : ", f8.4, " +/- ", f6.4)') fit(7), err(7)
      write (*, '("RMS Residual: ", f5.3)') rms

      do i=1,nrings
         x(1,i) = z(i)
         x(2,i) = r(i)
         x(3,i) = t(i)
         x(4,i) = dn(i)
      end do

c  plot non-linear residual vs z
      
      ff(1) = fit(1)
      ff(2) = fit(2)
      ff(3) = 0.0
      ff(4) = 0.0
      ff(5) = fit(5)
      ff(6) = fit(6)
      ff(7) = fit(7)
      ymin = 0.0
      ymax = 0.0
      do i=1,nrings
         call lambda (x(1,i), ff, 7, wfit, dwdf)
         yplot(i) = wave(i) - wfit
         if (iflag(i)) then
            ymin = min (ymin, yplot(i))
            ymax = max (ymax, yplot(i))
         end if
      end do
      scale = ymax - ymin
      ylo = ymin - 0.05 * scale
      yhi = ymax + 0.05 * scale

	call pgsci (1)
      call pgenv (zlo, zhi, ylo, yhi, 0, 0)
      call pgsch (1.25)
      call pglabel ("Z", "Wavelength Residual", "Non-Linear")
      call pgsch (1.50)
      do i=1,nrings
         if (iflag(i)) then
            call pgsci (3)
         else
            call pgsci (2)
            yplot(i) = max (ymin, yplot(i))
            yplot(i) = min (ymax, yplot(i))
         end if
         call pgpoint (1, z(i), yplot(i), 0)
      end do
      call pgsch (1.00)
      call pgsci (1)

      xinc = (zmax - zmin) / (NPLOT - 1)
      xx = zmin
      do i=1,NPLOT
         xplot(i) = xx
         yplot(i) = (xx**2 *(fit(3) + xx*fit(4))) / fit(7)
         xx = xx + xinc
      end do
      call pgline (NPLOT, xplot, yplot)


c  plot wave / cen wave vs r

      ymin = 1.0
      ymax = 1.0
      do i=1,nrings
         yplot(i) = (fit(1)+z(i)*(fit(2)+z(i)*(fit(3)+z(i)*fit(4)))
     $       + t(i)*fit(6)) / (fit(7)+dn(i))
         yplot(i) = yplot(i) / wave(i)
         if (iflag(i)) then
            ymin = min (ymin, yplot(i))
            ymax = max (ymax, yplot(i))
         end if
      end do
      ymax = max (ymax, sqrt(1.0 + rmax**2/fit(5)**2))
      scale = ymax - ymin
      ylo = ymin - 0.05 * scale
      yhi = ymax + 0.05 * scale

      call pgenv (rlo, rhi, ylo, yhi, 0, 0)
      call pgsch (1.25)
      call pglabel ("R", "Wave(0) / Wave(R)", "Radial")
      call pgsch (1.50)
      do i=1,nrings
         if (iflag(i)) then
            call pgsci (3)
         else
            call pgsci (2)
            yplot(i) = max (ymin, yplot(i))
            yplot(i) = min (ymax, yplot(i))
         end if
         call pgpoint (1, r(i), yplot(i), 0)
      end do
      call pgsch (1.00)
      call pgsci (1)

      xinc = (rmax) / (NPLOT - 1)
      xx = 0.0
      do i=1,NPLOT
         xplot(i) = xx
         yplot(i) = sqrt(1.0 + xx**2 / fit(5)**2)
         xx = xx + xinc
      end do
      call pgline (NPLOT, xplot, yplot)


c  plot a0 vs t

      ymin = 1.0e+06
      ymax = -1.0e+6
      do i=1,nrings
         yplot(i) = (wave(i)*sqrt(1.0+r(i)**2/fit(5)**2)*(fit(7)+dn(i))
     $        -(z(i)*(fit(2)+ z(i)*(fit(3) + z(i)*fit(4)))))
     $        / (fit(7))
         if (iflag(i)) then
            ymin = min (ymin, yplot(i))
            ymax = max (ymax, yplot(i))
         end if
      end do
      scale = ymax - ymin
      ylo = ymin - 0.05 * scale
      yhi = ymax + 0.05 * scale

      call pgenv (tlo, thi, ylo, yhi, 0, 0)
      call pgsch (1.25)
      call pglabel ("Time", "Wavelength", "Temporal Drift")
      call pgsch (1.50)
      do i=1,nrings
         if (iflag(i)) then
            call pgsci (3)
         else
            call pgsci (2)
            yplot(i) = max (ymin, yplot(i))
            yplot(i) = min (ymax, yplot(i))
         end if
         call pgpoint (1, t(i), yplot(i), 0)
      end do
      call pgsch (1.0)
      call pgsci (1)

      xinc = (tmax - tmin) / (NPLOT - 1)
      xx = tmin
      do i=1,NPLOT
         xplot(i) = xx
         yplot(i) = (fit(1) + fit(6) * xx) / fit(7)
         xx = xx + xinc
      end do
      call pgline (NPLOT, xplot, yplot)


c  plot residual vs z

      ymin = 0.0
      ymax = 0.0
      do i=1,nrings
         call lambda (x(1,i), fit, 7, wfit, dwdf)
         yplot(i) = wave(i) - wfit
         if (iflag(i)) then
            ymin = min (ymin, yplot(i))
            ymax = max (ymax, yplot(i))
         end if
      end do
      scale = ymax - ymin
      ylo = ymin - 0.05 * scale
      yhi = ymax + 0.05 * scale
      if (iplot .eq. 1) then
         xlo = zlo
         xhi = zhi
         xlabel = "Z"
      else if (iplot .eq. 2) then
         xlo = rlo
         xhi = rhi
         xlabel = "Radius"
      else
         xlo = tlo
         xhi = thi
         xlabel = "Time"
      end if
      xscale = xhi - xlo

      call pgenv (xlo, xhi, ylo, yhi, 0, 0)
      call pgsch (1.25)
      call pglabel (xlabel, "Wavelength Residual", "Full Fit Residuals")
      call pgsch (1.50)
      do i=1,nrings
         if (iflag(i)) then
            call pgsci (3)
         else
            call pgsci (2)
            yplot(i) = max (ymin, yplot(i))
            yplot(i) = min (ymax, yplot(i))
         end if
         if (iplot .eq. 1) then
            xplot(i) = z(i)
         else if (iplot .eq. 2) then
            xplot(i) = r(i)
         else
            xplot(i) = t(i)
         end if
         call pgpoint (1, xplot(i), yplot(i), 0)
      end do
      call pgsch (1.00)
      call pgsci (1)


 310  write (*, '("Enter option: "$)')
      read (*,'(a1)') answ

      if (answ .eq. "1") then
         flags(3) = .false.
         flags(4) = .false.
         fit(3) = 0.0
         fit(4) = 0.0
         go to 300
      else if (answ .eq. "2") then
         flags(3) = .true.
         flags(4) = .false.
         fit(4) = 0.0
         go to 300
      else if (answ .eq. "3") then
         flags(3) = .true.
         flags(4) = .true.
         go to 300
      else if ((answ .eq. "T") .or. (answ .eq. "t")) then
         if (.not. flags(6)) then
            flags(6) = .true.
         else
            flags(6) = .false.
            fit(6) = 0.0
         end if
         go to 300
      else if ((answ .eq. "Z") .or. (answ .eq. "z")) then
         iplot = 1
         go to 300
      else if ((answ .eq. "R") .or. (answ .eq. "r")) then
         iplot = 2
         go to 300
      else if ((answ .eq. "H") .or. (answ .eq. "h")) then
         iplot = 3
         go to 300
      else if ((answ .eq. "E") .or. (answ .eq. "e")) then
         write (*, '("Mark point with cursor")')
	   xcur = (xlo + xhi) / 2.0
	   ycur = (ylo + yhi) / 2.0
         ier = pgcurse (xcur, ycur, answ)
         dmin = 1.0e+6
         do i=1,nrings
            xd = (xcur - xplot(i))/ xscale
            yd = (ycur - yplot(i))/ scale
            dist = xd * xd + yd * yd
            if (dist .lt. dmin) then
               loc = i
               dmin = dist
            end if
         end do
         iflag(loc) = .NOT. iflag(loc)
         go to 300
      else if ((answ .eq. "Q") .or. (answ .eq. "q")) then
         go to 400
      else
         write (*,320)
 320     format ("Fit:  1 linear    2 quadratic    3 cubic    T time"/
     $        "Edit Point:  E"/
     $        "Residual Plot:  Z z-value    R radius    H time"/
     $        "Done:  Q")
         go to 310
      end if


c  output results

 400  continue
!	if (nargs .gt. 1) then
!         call getarg (2, filename)
!         len = lnblnk(filename)
!      else
         len = 0
         len = lnblnk(output)
         if (len .eq. 0 ) then
            print*, 'Error: Output file not specified'
            stop
         endif
c         do while (len .eq. 0)
c            write (*, '("Enter a name for the output file: "$)')
c            read (*, '(a)') output
c            len = lnblnk(output)
c         end do
!      end if

      ier = 1
      do while (ier .ne. 0)
!         open (11, file=filename, status="old", fileopt="eof", 
!     $        iostat=ier)
         open (20, file=output, status="old", iostat=ier)
         if (ier .ne. 0) then
            open (20, file=output, status="new", iostat=ier)
            if (ier .ne. 0) then
               write (*, '(a1,"Error opening output file ",a)')
     $              bell, output(1:len)
               len = 0
               do while (len .eq. 0)
                  write (*, '("Enter output filename: "$)')
                  read (*, '(a)') output
                  len = lnblnk(output)
               end do
            end if
         end if
      end do
 !     write (10, '("Output file: ",a)') filename(1:len)

!      label = "Wavelength Calibration Fitting    " // datetime
      label = "Wavelength Calibration Fitting    "
      write (20, '(a/)') label(1:lnblnk(label))

      write (20, 410) a, aerr, b, berr, c, cerr, d, derr, e, eerr,
     $     f, ferr
      write (*, 410) a, aerr, b, berr, c, cerr, d, derr, e, eerr,
     $     f, ferr
 410  format (10x, "Wave = (A + B*Z + C*Z**2 + D*Z**3) "
     $     "(1.0 + R**2/F**2)**(-1/2)"/
     $     25x, "A = ", f8.3, " +/- ", f5.3 /
     $     25x, "B = ", f8.4, " +/- ", f6.4 /
     $     25x, "C = ", f8.4, " +/- ", f6.4 /
     $     25x, "D = ", f8.4, " +/- ", f6.4 /
     $     25x, "E = ", f8.4, " +/- ", f8.4 /
     $     25x, "F = ", f8.2, " +/- ", f7.2 )
      if (multiorder) then
         write (20, 420) fit(7), err(7)
         write (*, 420) fit(7), err(7)
 420     format (25x, "Order  ", f7.3, " +/- ", f5.3)
      end if
      write (20, 430) rms
      write (*, 430) rms
 430  format (25x, "RMS fit residual: ", f5.3//
     $     4x, "Z", 7x, "R", 6x, "dR", 5x, "Xc", 6x, "Yc", 7x, "T",
     $     5x, "dN", 4x, "Wave", 7x, "Fit", 5x, "Error")

      do i=1,nrings
         call lambda (x(1,i), fit, 7, wfit, dwdf)
         resid = wave(i) - wfit
         idn = dn(i)
         if (iflag(i)) then
            write (20, 440) z(i), r(i), dr(i), xc(i), yc(i), t(i), 
     $           idn, wave(i), wfit, resid
            write (*, 440) z(i), r(i), dr(i), xc(i), yc(i), t(i), 
     $           idn, wave(i), wfit, resid
 440        format (1x, f6.3, 2x, f6.2, 2x, f5.2, 2(1x, f7.2), 2x, 
     $           f6.3, 2x, i3, 2(2x, f8.3), 2x, f6.3)
         else
            write (20, 450) z(i), r(i), dr(i), xc(i), yc(i), t(i), 
     $           idn, wave(i), wfit, resid
            write (*, 450) z(i), r(i), dr(i), xc(i), yc(i), t(i), 
     $           idn, wave(i), wfit, resid
 450        format ("(", f6.3, 2x, f6.2, 2x, f5.2, 2(1x, f7.2), 2x, 
     $           f6.3, 2x, i3, 2(2x, f8.3), 2x, f6.3, ")")
         end if
      end do

      call pgend ()
      close (20)
!      write (10, '("Successful completion")')
      close (10)
    

c  Error Termination

 666  continue
!	write (10, '("Error Termination")')
!      close (10)
      call pgend ()
      close(20)
      close(10)
      
      end subroutine calibrate




      subroutine lambda (x, p, npar, wave, dwdp)

c       independent variables: z = x(1)  r = x(2)  t = x(3)  dn = x(4)
c       parameters: a = p(1)  b = p(2)  c = p(3)  d = p(4)
c                   f = p(5)  e = p(6)  n = p(7)

      real  x(4), p(7), wave, dwdp(7)
      real  f1, f2, f3

      f1 = 1.0 / (p(7) + x(4))
      f2 = p(1) + x(1) * (p(2) + x(1) * (p(3) + x(1) * p(4)))
     $     + p(6) * x(3)
      f3 = 1.0 / sqrt(1.0 + x(2) * x(2) / (p(5) * p(5)))

      wave = f1 * f2 * f3

      dwdp(1) = f1 * f3
      dwdp(2) = dwdp(1) * x(1)
      dwdp(3) = dwdp(2) * x(1)
      dwdp(4) = dwdp(3) * x(1)
      dwdp(5) = wave * f3 * f3 * x(2) * x(2) / (p(5) * p(5) * p(5))
      dwdp(6) = dwdp(1) * x(3)
      dwdp(7) = - f1 * wave
      
      return
      end subroutine lambda 

