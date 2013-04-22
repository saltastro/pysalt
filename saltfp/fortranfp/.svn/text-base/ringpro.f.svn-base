c------------------------------------------------------------------------------
      subroutine ringpro (image, imgx, imgy, big, x, y, dy,
     >     axc,ayc,arad,rxc,ryc,icut,iwide,itmax,converg,wscale,rlo,rhi)
c------------------------------------------------------------------------------

c  This routine extracts and returns the radial profile of a PFIS FP 
c  calibration ring.  The geometrical parameters of the ring (aperture
c  center and radius, ring center) are specified in the PFPPAR common
c  block.  The ring image is contained in the real array IMAGE of dimension
c  IMGX by IMGY.  The logical parameter BIG controls the extraction mode: if 
c  true, the largest ring, with partial outermost arcs, is extracted; if 
c  false, the largest complete ring is extracted.  The ring profile is 
c  extracted in constant area annuli, with radii varying as r**(-2).  If the
c  parameter ICUT in PFPPAR is non-zero, the extracted profile is smoothed 
c  with a low-pass filter.  The profile is returned in the arrays X, Y, and DY.

c 220610 - nsl: Removed common block and added parameter parsing through the subroutine directly. This is to make the routine compatible with the salt pyraf package. 

      parameter (NRING=512, MAXPTS=1024)
      parameter (BADLO=-600.0)

      real image(imgx,imgy)
      real x(NRING), y(NRING+2), dy(NRING)
      real dat(MAXPTS,NRING)
      integer numb(NRING), isam(NRING)
      logical big

c      common /pfppar/ axc, ayc, arad, rxc, ryc, icut, iwide, itmax,
c     $                converg, wscale, rlo, rhi


	offset = sqrt((axc - rxc)**2 + (ayc - ryc)**2)
	if (big) then
	   rmax = arad + offset
	else
         rmax = arad - offset
	end if
      rmaxsq = rmax * rmax

      arsq = arad * arad
      isample = 1.10 + 3.14159 * arsq  / (MAXPTS * NRING)

      do i=1,NRING
         isam(i) = 0
         numb(i) = 0
      end do

      iylo = max(nint(ryc - rmax), 1)
      iyhi = min(nint(ryc + rmax), imgy)
      do iy=iylo,iyhi
         ry = iy - ryc
         rysq = ry * ry
         dx = sqrt(max(0.0,(rmaxsq - rysq)))
         ixlo = max (nint(rxc - dx), 1)
         ixhi = min (nint(rxc + dx), imgx)
         do ix=ixlo,ixhi
	      if (image(ix,iy) .gt. BADLO) then
               rx = ix - rxc
               rsq = rysq + rx * rx
               ir = 1.0 + NRING * rsq / rmaxsq
               if (ir .le. NRING) then
                  isam(ir) = isam(ir) + 1
                  if (mod(isam(ir),isample) .eq. 0) then
                     numb(ir) = numb(ir) + 1
                     dat(numb(ir),ir) = image(ix,iy)
                  end if
               end if
	      end if
         end do
      end do
 
      do i=1,NRING
	   x(i) = sqrt((i - 0.5) * rmaxsq / NRING)
         call biwgt (dat(1,i),numb(i),y(i),dy(i))
      end do

	if (icut .gt. 0) call lowpass (y, NRING, icut, iwide)

c        print*, axc, ayc, arad, rxc, ryc, filter,icut, iwide,plot,itmax,
c     >       converg, wscale,rlo,rhi,fixed


	return

      end 


