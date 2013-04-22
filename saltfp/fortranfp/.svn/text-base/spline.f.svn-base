c     Cubic Spline Interpolation from Numerical Recipes

      subroutine spline (x,y,n,y2)

      real x(n),y(n),y2(n),u(100)

c  use only a natural cubic spline

      y2(1) = 0.0
      y2(n) = 0.0
      u(1)  = 0.0

      do i=2,n-1
         sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
         p = sig * y2(i-1) + 2.0
         y2(i) = (sig - 1.0) / p
         u(i) = (6.0 * ((y(i+1) - y(i)) / (x(i+1) - x(i))
     $          - (y(i) - y(i-1)) / (x(i) - x(i-1)))
     $          / (x(i+1) - x(i-1)) - sig * u(i-1)) / p
      end do
      do i=n-1,1,-1
         y2(i) = y2(i) * y2(i+1) + u(i)
      end do

      return
      end


      subroutine splint (x,y,y2,n,xx,yy)

      real x(n),y(n),y2(n)

      data klo,khi /1,2/
      
      if ((x(klo) .gt. xx) .or. (x(khi) .lt. xx)) then
         klo = 1
         khi = n
         do while ((khi - klo) .gt. 1)
            k = (klo + khi) / 2
            if (x(k) .gt. xx) then
               khi = k
            else
               klo = k
            end if
         end do
      end if

      h = x(khi) - x(klo)
      a = (x(khi) - xx) / h
      b = (xx - x(klo)) / h
      yy = a*y(klo) + b*y(khi) + 
     $     ((a*a*a-a)*y2(klo) + (b*b*b-b)*y2(khi)) * h*h / 6.0

      return
      end
