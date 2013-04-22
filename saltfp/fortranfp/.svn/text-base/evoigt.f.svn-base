c------------------------------------------------------------------------------
      subroutine evfit (wave,flux,sigma,num,fit,flag,errs,chisq)
c------------------------------------------------------------------------------

c  This program fits a Voigt emission line profile to a spectrum.
c
c  Inputs:  wave(i)     array of wavelengths of the points
c           flux(i)     array of fluxes of the points
c           sigma(i)    array of flux errors of the points
c           num         number of points
c           fit(j)      initial guesses for the parameters
c           flag(j)     flags for parameters to fit
c
c  Outputs: fit(j)      fitted parameter values
c           errs(j)      error estimates for the parameters
c           chisq       reduced chi square for the fit
c
c  If an error occurs during fitting, chisq is returned with a negative value.


      real flux(num), wave(num), sigma(num), fit(5), errs(5)
      real covar(5,5)
      logical flag(5)
      external evoigt


      dof = float(num)
      do i=1,5
         if (flag(i)) dof = dof - 1.0
      end do

      if (dof .lt. 0.0) then
         chisq = -1.0
         return
      end if
         
      call mrqfit (wave,flux,sigma,num,fit,flag,5,
     $        covar,chisq,evoigt)      

      if (chisq .lt. 0.0) return
      do i=1,5
         if (flag(i)) then
            errs(i) = sqrt(abs(covar(i,i)))
         else
            errs(i) = 0.0
         end if
      end do
      if (dof .gt. 0.0) then
         chisq = chisq / dof
      else 
         chisq = 0.0
      end if

      return
      end


c------------------------------------------------------------------------------
      subroutine evinit (wave, flux, num, fit)
c------------------------------------------------------------------------------

c  This routine makes reasonable initial guesses for an emission line
c  Voigt profile fit.
c
c  Inputs:  wave(i)     array of wavelengths of the points
c           flux(i)     array of fluxes of the points
c           num         number of points
c
c  Output:  fit(j)      initial guesses for the parameters


      real wave(num), flux(num), fit(5)


      imax = 1
      wmin = wave(1)
      wmax = wave(1)
      do i=2,num
         if (flux(i) .gt. flux(imax)) imax = i
         wmin = min(wave(i), wmin)
         wmax = max(wave(i), wmax)
      end do
      fit(3) = wave(imax)

      dwmax = abs(wave(1) - fit(3))
      fit(1) = flux(1)
      do i=2,num
         dw = abs(wave(i) - fit(3))
         if (dw .gt. dwmax) then
            fit(1) = flux(i)
            dwmax = dw
         end if
      end do

      fhalf = (flux(imax) + fit(1)) / 2.0
      dwmax = 0.0
      do i=1,num
         if (flux(i) .gt. fhalf) then
            dw = abs(wave(i) - fit(3))
            if (dw .gt. dwmax) then
               dwmax = dw
            end if
         end if
      end do
      if (dwmax .eq. 0.0) dwmax = (wmax - wmin) / num
      fit(4) = dwmax / 3.532
      fit(5) = dwmax / 1.766

      y = fit(5) / (1.414 * fit(4))
      call voi (0.0, y, v, dvdx, dvdy)
      fit(2) = fit(4) * (flux(imax) - fit(1)) / v

      return
      end


c------------------------------------------------------------------------------
      subroutine evoigt (wave, a, npar, vgt, dvda)
c------------------------------------------------------------------------------

c  This subroutine evaluates a Voigt emission profile and its 
c  derivatives at a specified wavelength.  The parameters are:
c
c          a(1)  continuum level
c          a(2)  total line strength (area under Voigt)
c          a(3)  center wavelength
c          a(4)  Gaussian scale (sigma)
c          a(5)  Lorentzian scale
c
c  The Lorentzian scale (a(5)) may not be negative, so negative values
c  are set to 0.0 (a pure Gaussian).  The Gaussian scale (a(4)) must be 
c  greater than 0.0, so if it isn't, it is set to 0.0001*a(5) (essentially
c  a pure Lorentzian).  If both scales are out of bounds, the function is
c  not evaluated, and vgt is set to -666.0 as an error indicator.  
c  NPAR (included for compatibility with mrqfit) must be 5.


      real a(npar), dvda(npar)

      a(5) = max (0.0, a(5))
      if (a(4) .le. 0.0) a(4) = 0.0001 * a(5)
      if (a(4) .eq. 0.0) then
         vgt = -666.0
         return
      end if

      r2a4 = 1.0 / (1.4142136 * a(4))
      a2a4 = a(2) / a(4)
      x = (wave - a(3)) * r2a4
      y = a(5) * r2a4

      call voi(x, y, v, dvdx, dvdy)

      vgt = a(1) + a2a4 * v

      dvda(1) = 1.0
      dvda(2) = v / a(4)
      dvda(3) = - a2a4 * r2a4 * dvdx
      dvda(4) = - a2a4 * (v + x * dvdx + y * dvdy) / a(4)
      dvda(5) = a2a4 * r2a4 * dvdy

      return
      end


c------------------------------------------------------------------------------
      subroutine evstat (a, vcen, fwhm)
c------------------------------------------------------------------------------

c  This subroutine evaluates the central intensity and full width at
c  half maximum of a Voigt emission profile, given the Voigt parameters.
c  The parameters are:
c
c          a(1)  continuum level
c          a(2)  total line strength (area under Voigt)
c          a(3)  center wavelength
c          a(4)  Gaussian scale (sigma)
c          a(5)  Lorentzian scale


      real a(5)
 
      sl = max(0.0, a(5))
      sg = max(0.0001 * sl, a(4))
      if (sg .eq. 0.0) then
         vcen = -666.0
         fwhm = -666.0
         return
      end if

      r2sg = 1.0 / (1.4142136 * sg)
      a2sg = a(2) / sg
      y = sl * r2sg

      call voi (0.0, y, vcen, dvdx, dvdy)
      
      xlo = 0.0
      vlo = vcen

      xhi =(2.0 * (sl + sg)) * r2sg
      call voi (xhi, y, vhi, dvdx, dvdy)
      do while (vhi .gt. 0.5 * vlo)
         xhi = xhi * 2.0
         call voi (xhi, y, vhi, dvdx, dvdy)
      end do

      do while ((xhi - xlo) .gt. (0.001 * xhi))
         xnew = (xlo + xhi) * 0.5
         call voi (xnew, y, vnew, dvdx, dvdy)

         if (vnew .gt. 0.5 * vcen) then
            xlo = xnew
            vlo = vnew
         else
            xhi = xnew
            vhi = vnew
         end if
      end do

      vcen = vcen * a2sg
      fwhm = (xlo + xhi) / r2sg

      return
      end


c------------------------------------------------------------------------------
      subroutine voi (x,y,v,dvdx,dvdy)
c------------------------------------------------------------------------------

c     This routine computes the normalized Voigt function and its partial 
c     derivatives.  The integral of this function is the value of the 
c     Gaussian sigma of the Voigt function.
c
c     The computation is valid for  y > 0, and the maximum relative error
c     for the Voigt function is < 2.E-6
c
c     Subroutine adapted from:
c     J. Humlicek, J. Quant. Spectrosc. Radiat. Transfer 21, 309 (1979)
c     by Philip L. Varghese        <820315.0112>
c
c     Modified to compute derivatives by Ted Williams 5/21/92


      real t(6), c(6), s(6)

      data t/0.314240376E+00, 0.947788391E+00, 0.159768264E+01,
     *       0.227950708E+01, 0.302063703E+01, 0.388972490E+01/
      data c/0.403621095E+00,-0.299993210E+00, 0.500980825E-02,
     *       0.399820281E-02,-0.956712138E-04, 0.199809468E-07/
      data s/0.555821146E+00, 0.922164680E-01,-0.619762811E-01,
     *       0.248076921E-02, 0.366661062E-04,-0.250346637E-06/

      v = 0.0
      dvdx = 0.0
      dvdy = 0.0

      y1 = y + 1.5
      y2 = y1 * y1

C  branch to region I or II depending on values of X and Y.

      if ((y.gt.0.85).or.(abs(x).lt.(18.1*y+1.65))) then

C  calculations for region I

         do i=1,6
            rm = x - t(i)
            dm = 1.0 / (rm * rm + y2)
            d1 = y1 * dm
            d2 = rm * dm
            d1d2 = d1 * d2
            rp = x + t(i)
            dp = 1.0 / (rp * rp + y2)
            d3 = y1 * dp
            d4 = rp * dp
            d3d4 = d3 * d4
            v  = v + c(i) * (d1 + d3) - s(i) * (d2 - d4)
            dvdx = dvdx - 2.0 *(c(i) * (d1d2 + d3d4)
     $             + s(i) * (0.5*(dm-dp) - d2*d2 + d4*d4))
            dvdy = dvdy + 2.0 * (s(i) * (d1d2 - d3d4)
     $             + c(i) * (0.5*(dm+dp) - d1*d1 - d3*d3))
         end do
      else

C  calculations for region II

         if (abs(x).lt.12.0) then
            v = exp(-x*x) / 2.50662827
            dvdx = -2.0 * x * v / 2.50662827
         end if

         y3 = y + 3.0
         do i=1,6
            rm = x - t(i)
            sm = rm * rm
            dm = 1.0 / (sm + y2)
            d1 = y1 * dm
            d2 = rm * dm
            tm = 1.0 / (sm + 2.25)
            d2tm = d2 * tm
            dmtm = dm * tm
            fm =  tm * (s(i)*y3*d2 + c(i)*(rm*d2 - 1.5*d1))
            rp = x + t(i)
            sp = rp * rp
            dp = 1.0 / (sp + y2)
            d3 = y1 * dp
            d4 = rp * dp
            tp = 1.0 / (sp + 2.25)
            d4tp = d4 * tp
            dptp = dp * tp
            fp = tp * (-s(i)*y3*d4 + c(i)*(rp*d4 - 1.5*d3))
            v = v + y * (fm + fp)
            dvdx = dvdx + y * (s(i)*y3*(dmtm - dptp)
     $           + 2.0 * (c(i)*(d2tm + d4tp)
     $           - fm*(d2 + rm*tm) - fp*(d4 + rp*tp)))
            dvdy = dvdy + fm + fp + y * (s(i)*(d2tm - d4tp)
     $           - 1.5*c(i)*(dmtm + dptp)
     $           - 2.0*(d1*fm + d3*fp))
         end do
      end if

      return
      end


