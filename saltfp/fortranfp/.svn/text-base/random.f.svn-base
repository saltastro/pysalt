c------------------------------------------------------------------------------
      function ran1 (idum)
c------------------------------------------------------------------------------

c  This function returns a unifor random deviate between 0.0 and 1.0.
c  Set the argument, IDUM, to a negative value to initialize or reinitialize
c  the sequence.  Adapted from Numerical Recipes.

      implicit none

      real RM1, RM2
      integer M1, M2, M3, IA1, IA2, IA3, IC1, IC2, IC3
      parameter (M1=259200, IA1=7141, IC1=54773, RM1=1.0/M1)
      parameter (M2=134456, IA2=8121, IC2=28411, RM2=1.0/M2)
      parameter (M3=243000, IA3=4561, IC3=51349)

      integer idum, iff, ix1, ix2, ix3, j
      real ran1, r(97)

      data iff /0/

      if ((idum .lt. 0) .or. (iff .eq. 0)) then
         iff = 1
         ix1 = mod(ic1-idum, m1)
         ix1 = mod(ia1*ix1+ic1, m1)
         ix2 = mod(ix1, m2)
         ix1 = mod(ia1*ix1+ic1, m1)
         ix3 = mod(ix1, m3)
         do j=1,97
            ix1 = mod(ia1*ix1+ic1, m1)
            ix2 = mod(ia2*ix2+ic2, m2)
            r(j) = (float(ix1) + float(ix2) * rm2) * rm1
         end do
         idum = 1
      end if

      ix1 = mod(ia1*ix1+ic1, m1)
      ix2 = mod(ia2*ix2+ic2, m2)
      ix3 = mod(ia3*ix3+ic3, m3)

      j = 1 + (97 * ix3) / m3
      if (j .lt. 1) j = 1
      if (j .gt. 97) j = 97
      ran1 = r(j)
      r(j) = (float(ix1) + float(ix2) * rm2) * rm1

      return
      end

c------------------------------------------------------------------------------
      function gasdev (idum)
c------------------------------------------------------------------------------

c  This function returns a normally distributed deviate with zero mean and
c  unit variance.  Set the argument, IDUM, to a negative value to initialize 
c  or reinitialize the sequence.  Adapted from Numerical Recipes.

      implicit none

      real gasdev, v1, v2, r, fac, gset, ran1
      integer idum, iset

      data iset /0/

      if (iset .eq. 0) then
         r = 2.0
         do while (r .ge. 1.0)
            v1 = 2.0 * ran1(idum) - 1.0
            v2 = 2.0 * ran1(idum) - 1.0
            r = v1 * v1 + v2 * v2
         end do
         fac = sqrt(-2.0 * log(r) / r)
         gset = v1 * fac
         gasdev = v2 * fac
         iset = 1
      else
         gasdev = gset
         iset = 0
      end if

      return
      end

