"""A general module with tools for use with the saltfp package"""

def fpfunc(z, r, t, coef=None):
   """A functional form fitting the Fabry Perot parameterization.  The 
      FP parameterization is given by:

      $\lambda = \\frac{A+Bz+Cz^2+Dz^3}{(1+r/F)^{0.5}}+Et$

      Parameters
      ----------
      z: float or ndarray
         z position of the etalon
      r: float or ndarray
         r position in the image
      t: float or ndarray
         r position in the image
      coef: list or ndarray
         An array with the coefficients for the FP equations.  The
         coefficients should be given as [A B C D E F] 

      Returns
      -------
      w:  float or ndarray
         wavelength at position r
   """
   if len(coef)!=6: raise Exception('Not enough FP Coefficients')
   if coef[5]==0:  raise Exception('F must not be zero')

   #calcuate the value according to the above equation
   w=coef[0]+coef[1]*z+coef[2]*z**2+coef[3]*z**3+coef[4]*t
   w=w/(1+(r/coef[5])**2)**0.5
   return w

