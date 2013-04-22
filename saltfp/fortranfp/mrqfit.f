c------------------------------------------------------------------------------
      subroutine mrqfit (x, y, sig, ndata, a, afit, npar, 
     $     covar, chisq, funk)
c------------------------------------------------------------------------------

c  This routine determines the least-squares fit of a non-linear parametric
c  function to a set of data.
c
c  Inputs:
c    X     the x values of the data points; real, ndata
c    Y     the y values of the data points; real, ndata
c    SIG   the uncertainties of the y values of the data points; real, ndata
c    NDATA the number of data points; integer
c    A     initial estimates of the fit parameters; real, npar
c    AFIT  flags specifying which parameters to vary; logical, npar
c    FUNK  a subroutine to evaluate the function and its derivatives; external
c    NPAR  the number of parameters in the fit
c
c  Outputs:
c    A     the fitted parameter values
c    COVAR the covarience matrix of the fit; real, npar x npar
c    CHISQ chi squared for the fit; real


c	external funk

      real    CONVERGE
      integer MMAX
      parameter (MMAX=20, CONVERGE=0.001)

      real    x(ndata), y(ndata), sig(ndata)
      real    a(npar), covar(npar,npar)
      real    ocovar(MMAX,MMAX), oa(MMAX), da(MMAX), oda(MMAX)
      logical afit(npar), trouble, conver
      real funk

      nfit = 0
      do i=1,npar
         if (afit(i)) nfit = nfit + 1
      end do

      if (nfit .eq. 0) then
         chisq = 0.0
         return
      else if (nfit .gt. MMAX) then
         chisq = -1.0
         return
      end if

      alamda = 0.001
      call mrqcof  (x, y, sig, ndata, a, afit, npar, nfit,
     $     covar, da, npar, chisq, funk)
      conver = .false.

      do while (.not. conver)
         
         ochisq = chisq
         i = 0
         do k=1,npar
            if (afit(k)) then
               i = i + 1
               oa(i) = a(k)
               oda(i) = da(i)
               do j=1,nfit
                  ocovar(i,j) = covar(i,j)
               end do 
               covar(i,i) = covar(i,i) * (1.0 + alamda)
            end if
         end do
         
         call gaussj (covar, da, nfit, npar, trouble)
         if (trouble) then
            chisq = -2.0
            return
         end if
         
         i = 0
         do j=1,npar
            if (afit(j)) then
               i = i + 1
               a(j) = a(j) + da(i)
            end if
         end do
         
         call mrqcof (x, y, sig, ndata, a, afit, npar, nfit,
     $                covar, da, npar, chisq, funk)

         if (chisq .le. ochisq) then
            if ((ochisq - chisq) .le. CONVERGE * ochisq) then
               conver = .true.
            else
               alamda = 0.1 * alamda
            end if
         else
            alamda = 10.0 * alamda
            chisq = ochisq
            i = 0
            do k=1,npar
               if (afit(k)) then
                  i = i + 1
                  a(k) = oa(i)
                  da(i) = oda(i)
                  do j=1,nfit
                     covar(i,j) = ocovar(i,j)
                  end do
               end if
            end do
         end if

      end do

      call gaussj (covar, da, nfit, npar, trouble)
      if (trouble) then
         chisq = -3.0
         return
      end if
         
      if (nfit .lt. npar) then
         do i=nfit+1,npar
            do j=1,i
               covar(i,j) = 0.0
               covar(j,i) = 0.0
            end do
         end do

         k = nfit
         do j=npar,1,-1
            if (afit(j)) then
               do i=1,npar
                  swap       = covar(i,k)
                  covar(i,k) = covar(i,j)
                  covar(i,j) = swap
               end do
               do i=1,npar
                  swap       = covar(k,i)
                  covar(k,i) = covar(j,i)
                  covar(j,i) = swap
               end do
               k = k - 1
            end if
         end do
      end if

      return
      end

c------------------------------------------------------------------------------
      subroutine mrqcof (x, y, sig, ndata, a, afit, npar, nfit,
     $                    alpha, beta, nalp, chisq, funk)
c------------------------------------------------------------------------------

c  This subroutine calculates the coefficients and constants of the normal
c  equations for the improvements to the parameters in a non-linear
c  least-squares fit.  In addition, the chi-squared for the current fit
c  is calculated.
c
c  Inputs:
c    X     the x values of the data points; real, ndata
c    Y     the y values of the data points; real, ndata
c    SIG   the uncertainties of the y values of the data points; real, ndata
c    NDATA the number of data points; integer
c    A     current estimates of the fit parameters; real, npar
c    AFIT  flags specifying which parameters to vary; logical, npar
c    NPAR  the number of parameters in the fit
c    NFIT  the number of parameters to vary in the fit
c    NALP  the dimension of the alpha and beta arrays
c    FUNK  a subroutine to evaluate the function and its derivatives; external
c
c  Outputs:
c    ALPHA the coefficient matrix of the fit; real, npar x npar
c    BETA  the constant vector of the fit; real, npar
c    CHISQ chi squared for the fit; real


      parameter (MMAX=20)

      real     x(ndata), y(ndata), sig(ndata)
      real     a(npar), alpha(nalp,nalp), beta(nalp)
      real     dyda(MMAX)
      logical  afit(npar)
      external funk

      do i=1,nfit
         do j=1,i
            alpha(i,j) = 0.0
         end do
         beta(i) = 0.0
      end do

      chisq = 0.0
      do i=1,ndata
         call funk (x(i), a, npar, yfit, dyda)
         wgt = 1.0 / (sig(i) * sig(i))
         dy = y(i) - yfit

         j = 0
         do l=1,npar
            if (afit(l)) then
               j = j + 1
               tmp = dyda(l) * wgt
               k = 0
               do m=1,l
                  if (afit(m)) then
                     k = k + 1
                     alpha(j,k) = alpha(j,k) + tmp * dyda(m)
                  end if
               end do
               beta(j) = beta(j) + dy * tmp
            end if
         end do
         chisq = chisq + dy * dy * wgt
      end do

      do i=2,nfit
         do j=1,i-1
            alpha(j,i) = alpha(i,j)
         end do
      end do

      return
      end

