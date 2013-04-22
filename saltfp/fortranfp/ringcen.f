c------------------------------------------------------------------------------
      subroutine ringcen (image,imgx,imgy,fixed,rave,rerr,rxc,ryc,
     > axc,ayc,arad,rxcen,rycen,icut,iwide,itmax,converg,wscale,rlo,rhi)
c------------------------------------------------------------------------------

c  This routine measures the center coordinates and radius of a
c  PFIS FP calibration ring.  


      parameter (NRING=512, NSEG=4, MAXPTS=1024, MAXDIV=3)
      parameter (BADLO=-600.0)

      real image(imgx,imgy)
      real x(NRING), y(NRING+2), dy(NRING)
      real dat(MAXPTS,NRING,NSEG)
      real a(5), da(5), rad(NSEG)
      integer numb(NRING,NSEG), isam(NRING,NSEG)
      logical flag(5), conver, fixed

c      common /pfppar/ axc, ayc, arad, rxcen, rycen, icut, iwide, itmax,
c     $                converg, wscale, rlo, rhi


c  initialization

      do i=1,5
         flag(i) = .true.
      end do

	rxc = rxcen
	ryc = rycen
      arsq = arad * arad
      isample = 1.10 + 3.14159 * arsq  / (NSEG * MAXPTS * NRING)
      rmax = arad - sqrt ((axc - rxc)**2 + (ayc - ryc)**2)
      rmaxsq = rmax * rmax

c  determine initial ring parameters

	call ringpro (image,imgx,imgy,.false.,x,y,dy,
     > axc,ayc,arad,rxcen,rycen,icut,iwide,itmax,converg,wscale,rlo,rhi)
      do i=1,NRING
         x(i) = (i - 0.5) / NRING
      end do

	if (icut .gt. 0) call lowpass (y, NRING, icut, iwide)

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
      call evfit (x(ilo), y(ilo), dy(ilo), num, a, flag, da, chisq)
	if (chisq .lt. 0.0 ) then
	   rave = -1.0
	   return
	end if
      call evstat (a, ycen, fwhm)
      xcen = a(3) * NRING
      cont = a(1)
      width = fwhm * NRING * wscale

c  loop to measure centers and radius

      iter = 0
      conver = .false.
      oldstep = 1000.0
      ndiv = 0

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
         iyhi = min (nint(ryc + rmax), imgy)
         do iy=iylo,iyhi
            ry = iy - ryc
            rysq = ry * ry
            dx = sqrt (max(1.0,(rmaxsq - rysq)))
            ixlo = max (nint(rxc - dx), 1)
            ixhi = min (nint(rxc + dx), imgx)
            do ix=ixlo,ixhi
	         if (image(ix,iy) .gt. BADLO) then
                  rx = ix - rxc
                  rsq = rysq + rx * rx
                  ir = 1.0 + NRING * rsq / rmaxsq
                  if (ir .le. NRING) then
                     iseg = (atan2(ry,rx) + 3.14159)*NSEG/6.28319 + 1
                     iseg = max(min(iseg, NSEG), 1)      
                     isam(ir,iseg) = isam(ir,iseg) + 1
                     if (mod(isam(ir,iseg),isample) .eq. 0) then
                        numb(ir,iseg) = numb(ir,iseg) + 1
                        dat(numb(ir,iseg),ir,iseg) = image(ix,iy)
                     end if
                  end if
	         end if
            end do
         end do
 
         rave = 0.0
         rerr = 0.0
         dxc  = 0.0
         dyc  = 0.0

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
	   xc = xcen - ilo + 1.0
         
         do j=1,NSEG 
            do i=1,NRING
               call biwgt (dat(1,i,j), numb(i,j), y(i), dy(i))
            end do
            
	      if(icut .gt. 0) call lowpass (y, NRING, icut, iwide)

            call center (y(ilo), num, xc, width, cont)
            if (xc .lt. 0.6) then
 			 rave = -2.0
               return
            end if
	      rc = xc + ilo - 1.0

            rad(j) = rmax * sqrt ((rc - 0.5) / NRING)
            rave = rave + rad(j)
            rerr = rerr + rad(j) * rad(j)
            ang = 6.28319 * (j - 0.5) / NSEG
            dxc = dxc + rad(j) * cos (ang)
            dyc = dyc + rad(j) * sin (ang)
         end do

	   if (fixed) then
	      dxc = 0.0
	      dyc = 0.0
	   else
            dxc = 2.0 * dxc / NSEG
            dyc = 2.0 * dyc / NSEG
	   end if
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
         conver = step .lt. (converg * max(1.0, rerr))
         if (step .lt. oldstep) then
            oldstep = step
            ndiv = 0
         else
            ndiv = ndiv + 1
         end if

         iter = iter + 1

      end do

	if (.not. conver) rave = - rave

	return
      end


