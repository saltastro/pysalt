c------------------------------------------------------------------------------
      subroutine center (y, ndata, xcen, width, cont)
c------------------------------------------------------------------------------

c  This routine determines the center of an emission feature in a one
c  dimensional data set with uniform sample spacing, using the IRAF
c  CENTER1D algorithm with linear interpolation.
c
c  Inputs:
c    Y     the y values of the data points; real, ndata
c    NDATA the number of data points; integer
c    XCEN  initial estimate of the center; real
c    WIDTH estimated width of the feature; real
c    CONT  continuum estimate; real
c
c  Output:
c    XCEN  the improved center estimate; real
c

      parameter (CONVERGE=0.005, EPSILON=0.01, RAD=10.0)
      parameter (MMAX=1024, WMIN=3.0, ITMAX=50, MCHECK=3)

      real y(ndata)
      real dataa(MMAX), datai(MMAX)


c  determine data range

      ilo = max (1, int(xcen - width/2.0 - RAD))
      ihi = min (ndata, int(xcen + width/2.0 + RAD + 1.0))
      ndat = ihi - ilo + 1

      j = 0
      do i=ilo,ihi
         j = j + 1
         dataa(j) = max (0.0, y(i) - cont)
         datai(j) = dataa(j) * j
      end do
      xcen = xcen - ilo + 1.0

      icheck = 0
      dxlast = ndat

      do iter=1,ITMAX

         hwidth = width / 2.0
         hwidth = min(hwidth, xcen - 1.0)
         hwidth = min(hwidth, ndat - xcen - 1.0)
         if (hwidth .lt. WMIN) then
            xcen = -1.0
            return
         end if

         a = xcen - hwidth
         b = xcen - hwidth / 2.0
         call integrate (dataa, ndat, a, b, area1)
         call integrate (datai, ndat, a, b, area2)
         sum1 = (xcen - hwidth) * area1 - area2
         sum2 = - area1

         a = b
         b = xcen + hwidth / 2.0
         call integrate (dataa, ndat, a, b, area1)
         call integrate (datai, ndat, a, b, area2)
         sum1 = sum1 - xcen * area1 + area2
         sum2 = sum2 + area1

         a = b
         b = xcen + hwidth
         call integrate (dataa, ndat, a, b, area1)
         call integrate (datai, ndat, a, b, area2)
         sum1 = sum1 + (xcen + hwidth) * area1 - area2
         sum2 = sum2 - area1

         if (sum2 .eq. 0.0) then
            xcen = -2.0
            return
         end if

         dx = sum1 / abs(sum2)
         dxabs = abs (dx)
         dx = max (-1.0, min (1.0, dx))
         xcen = xcen + dx

         if (dxabs .lt. CONVERGE) then
            xcen = xcen + ilo - 1.0
            return
         else if (dxabs .gt. (dxlast + EPSILON)) then
            icheck = icheck + 1
            if (icheck .gt. MCHECK) then
               xcen = -3.0
               return
            end if
         else if (dxabs .gt. (dxlast - EPSILON)) then
            xcen = xcen - dx / 2.0
            icheck = 0
         else
            icheck = 0
            dxlast = dxabs
         end if

      end do

      xcen = -4.0

      return
      end

c------------------------------------------------------------------------------
      subroutine integrate (y, ndata, a, b, area)
c------------------------------------------------------------------------------

c  This routine calculates the integral under a function between specified
c  limits, using linear interpolation.
c
c  Inputs:
c    Y     the y values of the data points; real, ndata
c    NDATA the size of the data array; integer
c    A     the lower limit; real
c    B     the upper limit; real
c
c  Output:
c    AREA  the calculated integral; real
c
c  If the limits exceed the data range, the integral is returned as 0.0


      real y(ndata)


      area = 0.0

      if (a .eq. b) then
         return
      else if (a .lt. b) then
         xa = a
         xb = b
      else
         xa = b
         xb = a
      end if

      na = xa
      nb = xb
      if ((na .lt. 1) .or. (nb .ge.ndata)) return

      s = xa - na
      t = xb - nb
      area = 0.5 * ((t*(2.0-t)+1.0)*y(nb) - (s*(2.0-s)+1.0)*y(na)
     $       + t*t*y(nb+1) - s*s*y(na+1))

      if (na .ne. nb) then
         do i=na,nb-1
            area = area + y(i)
         end do
      end if

      if (a .gt. b) area = - area

      return
      end

