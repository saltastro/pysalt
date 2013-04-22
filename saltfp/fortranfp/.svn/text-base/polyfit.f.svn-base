c------------------------------------------------------------------------------
      subroutine polyfit (x, y, n, fit, m, trouble)
c------------------------------------------------------------------------------

c  This routine fits a polynomial to a data set.  Note that M=1 is a constant, 
c  M=2 is a linear fit, etc.  All data points receive equal weight.
c  Uncertainties of the polynomial coefficients are not estimated.  
c  Adapted from the code in Numerical Recipes.
c
c  Inputs:
c    X  the vector of independent vlaues; real, n
c    Y  the vector of dependent values; real, n
c    N  the number data points in X and Y; integer
c    M  the order of polynomial to be fit; integer
c
c  Outputs:
c    FIT the polynomial coefficients
c    TROUBLE false for normal termination, true if an error occurs; logical



      parameter (MMAX=20)

      real    x(n), y(n), fit(m)
	real	covar(MMAX,MMAX), p(MMAX)
      logical trouble


      trouble = .false.

      if (m .gt. min(MMAX,n)) then
         trouble = .true.
         return
      end if

	do j=1,m
		fit(j) = 0.0
		do i=1,m
			covar(i,j) = 0.0
		end do
	end do


	do i=1,n
		p(1)=1.0
		if (m .gt. 1) then
			do j=2,m
				p(j) = p(j-1) * x(i)
			end do
		end if
		do j=1,m
			fit(j) = fit(j) + p(j) * y(i)
			do k=1,j
				covar(j,k) = covar(j,k) + p(j) * p(k)
			end do
		end do
	end do

	if (m .gt. 1) then
		do j=2,m
			do k=1,j-1
				covar(k,j) = covar(j,k)
			end do
		end do
	end if

	call gaussj (covar, fit, m, MMAX, trouble)

	return
	end


	function poly (x, fit, m)

	real fit(m)

	poly = fit(m)
	if (m .gt. 1) then
		do i=m-1,1,-1
			poly = poly * x + fit(i)
		end do
	end if

	return
	end
