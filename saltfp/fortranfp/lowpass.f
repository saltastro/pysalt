c------------------------------------------------------------------------------
      subroutine lowpass (data, num, icut, iwide)
c------------------------------------------------------------------------------

c  This routine applies a low-pass filter to a data set.  The filter
c  function smoothly tapers from 1.0 at frequency (icut-iwide) to 0.0
c  at frequency (icut+iwide).
c
c  Inputs:  data        array of data; real, size numb+2
c           num         number of points; integer, must be a power of 2
c           icut        filter frequency; integer
c           iwide       filter cutoff width; integer
c
c  Outputs: data        the filtered data
c


      real data(num+2)
 
      num2 = num / 2
      loc = min (max (icut, 2), num2)
      ilo = max (loc - iwide, 2)
      ihi = min (loc + iwide, num2)
      step = 3.14159 / (ihi - ilo + 2)

      afit = data(1)
      bfit = (data(num) - data(1)) / (num - 1)
      do i=1,num
         data(i) = data(i) - afit - bfit * (i - 1)
      end do

      call realft (data, num2, 1)

      loc = ilo * 2 - 1
      arg = step
      do i=ilo,ihi
         filt = 0.5 * (1.0 + cos(arg))
         arg = arg + step
         data(loc) = data(loc) * filt
         data(loc+1) = data(loc+1) * filt
         loc = loc + 2
      end do

      if (ihi .lt. num2) then
         do i=ihi+1,num2
            data(loc) = 0.0
            data(loc+1) = 0.0
            loc = loc + 2
         end do
      end if

      data(2) = 0.0

      call realft (data, num2, -1)
            
      do i=1,num
         data(i) = (data(i) / num2) + afit + bfit * (i - 1)
      end do

      return
      end


