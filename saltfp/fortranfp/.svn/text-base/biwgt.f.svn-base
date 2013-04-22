c----------------------------------------------------------------------------- 
      subroutine biwgt (x, n, xbl, xbs)
c----------------------------------------------------------------------------- 

c  This subroutine calculates the Biweight estimator of the location and
c  scale of a data set.  This is a robust estimator, and is preferred to
c  the median by all rational folks.  The Biweight function is given by:
c
c                                  u((1-u*u)**2)     abs(u) <= 1
c                         f(u) =
c                                  0                 abs(u) >  1
c
c  where u is defined by
c
c                         u = (X(I) - Med) / (c * Mad)  .
c
c  Med, Mad, and c are the median, the median absolute deviation from
c  the median, and the tuning constant respectively. The tuning constant
c  is a parameter which is chosen depending on the sample size and the 
c  specific function being used for the scale estimate.  Here the tuning
c  constant is set to 6.0 for calculation of the location, and to 9.0 for
c  calculation of the scale, as recommended by Tukey.
c
c  The input data are provided in the single precision real array X, of 
c  length N.  The location and scale are returned in the real variables 
c  XBL and XBS, respectively.  On input errors, the scale XBS is returned
c  with negative values:
c
c        for N < 1     XBS = -2.0  and  XBL = 0.0
c        for N = 1     XBS = -1.0  and  XBL = X(1)
c
c  WARNING: This version of the Biweight routine has the side effect of
c  sorting the input array X, and returning the sorted data in place.



      real*8  sum1, sum2, t0, t1, cmad, cmadsq
      real    x(n)


      if (n .lt. 1) then
         xbl =  0.0
         xbs = -2.0
         return
      else if (n .eq. 1) then
         xbl =  x(1)
         xbs = -1.0
         return
      else if (n .eq. 2) then
         xbl = 0.5 * (x(1) + x(2))
         xbs = abs(x(1) - xbl)
         return
      end if

      call medmad (x, n, xmed, xmad)

      if (xmad .eq. 0.0) then
         xbl = xmed
         xbs = xmad
         return
      end if

      xbl = xmed
      delta = xmad
      cmad = 6.0 * xmad
      cmadsq = cmad * cmad
      conv = 0.10 * xmad / sqrt(n - 1.0)
      iter = 0

      do while ((abs(delta) .ge. conv) .and. (iter .lt. 5))

         sum1 = 0.0
         sum2 = 0.0

         do i=1,n
            t0 = x(i) - xbl
            if (abs(t0) .lt. cmad) then
               t1 = cmadsq - t0 * t0
               t1 = t1 * t1
               sum1 = sum1  + (t0 * t1)
               sum2 = sum2 + t1
            endif
         end do
         
         delta = sum1 / sum2
         xbl = xbl + delta
         iter = iter + 1

      end do

      sum1 = 0.0
      sum2 = 0.0
      cmad = 9.0 * xmad
      cmadsq = cmad * cmad
      
      do i=1,n
         t0 = x(i) - xbl
         if (abs(t0) .lt. cmad) then
            t0 = t0 * t0
            t1 = cmadsq - t0
            sum1 = sum1 + (t0 * t1 * t1 * t1 * t1)
            sum2 = sum2 + t1 * (cmadsq - 5.0 * t0)
         endif
      end do

      xbs = n * sqrt(sum1/(n-1.0)) / abs(sum2)

      return
      end

c------------------------------------------------------------------------------
      subroutine medmad (x, n, xmed, xmad)
c------------------------------------------------------------------------------

c  This routine calculates the median and the median absolute deviation of
c  a data set.  The data are supplied in the array X, of length N. The
c  median is returned in XMED, and the median absolute deviation in XMAD.
c
c  If N is odd, the median is the value from the data set that has equal
c  numbers of values greater than it and less than it.  If N is even, the
c  median is the mean of the two central values of the data set.
c
c  The median absolute deviation, as the name implies, is the median of
c  the distribution of absolute values of the deviation of the data about
c  the median.
c
c  On input errors, XMAD is returned with negative values:
c       for N < 1     XMAD = -2.0  and  XMED undefined
c       for N = 1     XMAD = -1.0  and  XMED = X(1)
c
c  WARNING: This routine has the side effect of sorting the input array X.



      real x(n)
      logical even


      if (n .lt. 1) then
         xmad = -2.0
         return
      else if (n .eq. 1) then
         xmed = x(1)
         xmad = -1.0
         return
      else if (n .eq. 2) then
         xmed = 0.5 * (x(1) + x(2))
         xmad = abs(x(1) - xmed)
         return
      end if

      even = mod(n,2) .eq. 0
      n2 = n / 2

c  sort the data

      call sort (n,x)

c  calculate the median

      if (even) then
         xmed = 0.5 * (x(n2) + x(n2+1))
      else
         xmed = x(n2+1)
      endif

c  calculate the mad

      i = n2
      j = n2 + 2
      if (even) j = n2 + 1
      xi = xmed - x(i)
      xj = x(j) - xmed

      if (n .eq. 3) then
         xmad = min (xi, xj)
         return
      end if

      do k=1,n2
         if (xi .lt. xj) then
            xmad = xi
            if (i .gt. 1) then
               i = i - 1
               xi = xmed - x(i)
            else
               j = j + 1
               xj = x(j) - xmed
            end if
         else
            xmad = xj
            if (j .lt. n) then
               j = j + 1
               xj = x(j) - xmed
            else
               i = i - 1
               xi = xmed - x(i)
            end if
         end if
      end do

      if (even) then
         if (xi .lt. xj) then
            xmad = 0.5 * (xmad + xi)
         else
            xmad = 0.5 * (xmad + xj)
         end if
      end if

      return
      end

c------------------------------------------------------------------------------
      subroutine sort (n, x)
c------------------------------------------------------------------------------
 
c  This routine performs a heapsort of the data in array X, of length N.
c  The data are sorted into ascending order, and returned in the array X.
c
c  Stolen (unabashedly) from NUMERICAL RECIPES, with small n speed-ups
c  added



      real x(n)

      if (n .le. 3) then
         if (n .le. 1) return
         if (x(1) .gt. x(2)) then
            temp = x(1)
            x(1) = x(2)
            x(2) = temp
         end if
         if (n .eq. 2) return
         if (x(2) .gt. x(3)) then
            temp = x(2)
            x(2) = x(3)
            x(3) = temp
         end if
         if (x(1) .gt. x(2)) then
            temp = x(1)
            x(1) = x(2)
            x(2) = temp
         end if
         return
      end if

      m = n / 2 + 1
      k = n

 10   if (m .gt. 1) then
         m = m - 1
         temp = x(m)
      else
         temp = x(k)
         x(k) = x(1)
         k = k - 1
         if (k .eq. 1) then
            x(1) = temp
            return
         endif
      endif

      i = m
      j = m + m

 20   if (j .le. k) then
         if (j .lt. k) then
            if (x(j) .lt. x(j+1)) j = j + 1
         endif
         if (temp .lt. x(j)) then
            x(i) = x(j)
            i = j
            j = j + j
         else
            j = k + 1
         endif
         go to 20
      endif

      x(i) = temp
      go to 10

      end

