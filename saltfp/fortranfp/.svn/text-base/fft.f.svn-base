C     Fast Fourier Transform routines from Numerical Recipes

      subroutine realft (data,n,iflag)

      real*8 wr, wi, wpr, wpi, wtemp, theta
      real data(2*n+2)

      theta = 3.14159265358979d0 / dfloat(n)
      wr = 1.0d0
      wi = 0.0d0
      c1 = 0.5
      if (iflag .gt. 0) then
         c2 = -0.5
         call four1 (data, n, iflag)
         data(2*n+1) = data(1)
         data(2*n+2) = data(2)
      else
         c2 = 0.5
         theta = -theta
         data(2*n+1) = data(2)
         data(2*n+2) = 0.0
         data(2) = 0.0
      end if
      wpr = -2.0d0 * dsin(0.5d0 * theta) ** 2
      wpi = dsin(theta)
      n2p3 = 2 * n + 3
      do i=1,n/2+1
         i1 = 2 * i - 1
         i2 = i1 + 1
         i3 = n2p3 - i2
         i4 = i3 + 1
         wrs = sngl(wr)
         wis = sngl(wi)
         h1r =  c1 * (data(i1) + data(i3))
         h1i =  c1 * (data(i2) - data(i4))
         h2r = -c2 * (data(i2) + data(i4))
         h2i =  c2 * (data(i1) - data(i3))
         data(i1) =  h1r + wrs * h2r - wis * h2i
         data(i2) =  h1i + wrs * h2i + wis * h2r
         data(i3) =  h1r - wrs * h2r + wis * h2i
         data(i4) = -h1i + wrs * h2i + wis * h2r
         wtemp = wr
         wr = wr * wpr - wi * wpi + wr
         wi = wi * wpr + wtemp * wpi + wi
      end do
      if (iflag .gt. 0) then
         data(2) = data(2*n+1)
      else
         call four1 (data, n, iflag)
      end if

      return
      end


      subroutine four1 (data, nn, iflag)

      real*8 wr, wi, wpr, wpi, wtemp, theta
      real data(2*nn)

      n = 2 * nn
      j = 1
      do i=1, n, 2
         if (j .gt. i) then
            tempr = data(j)
            tempi = data(j+1)
            data(j) = data(i)
            data(j+1) = data(i+1)
            data(i) = tempr
            data(i+1) = tempi
         end if
         m = n / 2
         do while ((m .ge. 2) .and. (j .gt. m))
            j = j - m
            m = m / 2
         end do
         j = j + m
      end do
      mmax = 2
      do while (n .gt. mmax)
         istep = 2 * mmax
         theta = 6.28318530717959d0 / (iflag * mmax)
         wpr = -2.0d0 * dsin(0.5d0 * theta) ** 2
         wpi = dsin(theta)
         wr = 1.0d0
         wi = 0.0d0
         do m=1,mmax,2
            wrs = sngl (wr)
            wis = sngl (wi)
            do i=m, n, istep
               j = i + mmax
               tempr = wrs * data(j)   - wis * data(j+1)
               tempi = wrs * data(j+1) + wis * data(j)
               data(j)   = data(i)   - tempr
               data(j+1) = data(i+1) - tempi
               data(i)   = data(i)   + tempr
               data(i+1) = data(i+1) + tempi
            end do
            wtemp = wr
            wr = wr * wpr - wi    * wpi + wr
            wi = wi * wpr + wtemp * wpi + wi
         end do
         mmax = istep
      end do

      return
      end

