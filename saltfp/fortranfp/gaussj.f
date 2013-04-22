c------------------------------------------------------------------------------
      subroutine gaussj (a, b, n, np, trouble)
c------------------------------------------------------------------------------

c  This routine performs a matrix inversion and solution of the corresponding
c  set of linear equations using the Gauss-Jordan method with maximum pivot
c  strategy.  Adapted from the code in Numerical Recipes (this version does
c  not solve for multiple right-hand sides).
c
c  Inputs:
c    A  the matrix of coefficients; real, np x np
c    B  the right-hand side of the equations; real, np
c    N  the number of equations; integer
c    NP the dimension of the arrays A and B; integer
c
c  Outputs:
c    A  the inverse of the matrix
c    B  the solution vector
c    TROUBLE false for normal termination, true if an error occurs; logical


      parameter (NMAX=50)

      real    a(np,np), b(np)
      integer ipiv(NMAX), indxr(NMAX), indxc(NMAX)
      logical trouble


      trouble = .false.

      if (n .gt. NMAX) then
         trouble = .true.
         return
      end if

      do i=1,n
         ipiv(i) = 0
      end do

      do i=1,n
         temp = 0.0
         do j=1,n
            if (ipiv(j) .ne. 1) then
               do k=1,n
                  if (ipiv(k) .eq. 0) then
                     if (abs(a(j,k)) .ge. temp) then
                        temp = abs(a(j,k))
                        irow = j
                        icol = k
                     end if
                  else if (ipiv(k) .gt. 1) then
                     trouble = .true.
                     return
                  end if
               end do
            end if
         end do

         ipiv(icol) = ipiv(icol) + 1
         if (irow .ne. icol) then
            do j=1,n
               temp = a(irow,j)
               a(irow,j) = a(icol,j)
               a(icol,j) = temp
            end do
            temp = b(irow)
            b(irow) = b(icol)
            b(icol) = temp
         end if

         indxr(i) = irow
         indxc(i) = icol
         if (a(icol,icol) .eq. 0.0) then
            trouble = .true.
            return
         end if

         temp = 1.0 / a(icol,icol)
         a(icol,icol) = 1.0
         do j=1,n
            a(icol,j) = a(icol,j) * temp
         end do
         b(icol) = b(icol) * temp

         do j=1,n
            if (j .ne. icol) then
               temp = a(j,icol)
               a(j,icol) = 0.0
               do k=1,n
                  a(j,k) = a(j,k) - a(icol,k) * temp
               end do
               b(j) = b(j) - b(icol) * temp
            end if
         end do
      end do

      do i=n,1,-1
         if (indxr(i) .ne. indxc(i)) then
            do j=1,n
               temp = a(j,indxr(i))
               a(j,indxr(i)) = a(j,indxc(i))
               a(j,indxc(i)) = temp
            end do
         end if
      end do

      return
      end

