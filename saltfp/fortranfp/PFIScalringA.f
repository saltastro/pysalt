	program calring

c  This routine measures the radius and center of a PFIS FP calibration ring.


	use dfport

      parameter (NRING=512, MAXPTS=1024)

	integer*2 image(2148,2051)
      real dat(MAXPTS,NRING)
      real x(NRING), y(NRING), dy(NRING), fit(NRING)
      real a(5), da(5),  xv(2), yv(2)
      integer numb(NRING), isam(NRING)
      logical flag(5)
	character*80 name, string, comment
	character*1 answ


	call pgbegin (0, "/WS", 1, 1)
	call pgask (.false.)

! position to the directory

	write (*,1000)
1000	format ("Enter directory path name: "$)
	read (*,1001) name
1001	format (a)
	istat = chdir(name)
	if (istat .ne. 0) then
		write (*,1002)
1002		format ("Error changing directory")
		stop
	end if

	open (10, file="calring.dat", access='APPEND')


!  read the input image file

100	continue
	write (*,1004) 
1004	format ("Enter PFIS image filename: "$)
	read (*,1001) name
	istat = 0
	call ftopen (11, name, 1, iblk, istat)
	if (istat .ne. 0) then
		write (*,1005) istat
1005		format ("*** Error opening file, status = ",i)
		go to 101
	end if
	call ftgkys (11, "CCDSUM", string, comment, istat)
	if (istat .eq. 0) then
		read (string,*) ixbin, iybin
	else
		write (*,1010)
1010		format ("Error: cannot determine binning")
		istat = 0
		call ftclos (11, istat)
		go to 101
	end if
	call ftgkys (11, "UTC-OBS", string, comment, istat)
	if (istat .eq. 0) then
		read (string,1011) ihr, imn, isec
1011		format (i2,1x,i2,1x,i2)
		ut = ihr + imn/60.0 + isec/3600.0
	else
		write (*,1012)
1012		format ("Error: cannot determine time")
		istat = 0
		call ftclos (11, istat)
		go to 101
	end if
	call ftgpvi (11, 0, 1, 4405548, 0, image(1,1), ext, istat)
	call ftclos (11, istat)

	call calringf (image,2148,2051,2,rave,rerr,rxc,ryc) 

c  write output data

      write (*, 1006) rave, rerr, rxc, ryc
1006	format ("Radius: ",f7.2, " +/- ", f6.2, 5x, "X center: ",
     $     f7.2, 5x, "Y center: ", f7.2)
	i=lnblnk(name)
      write (10, 1007) rave, rerr, rxc, ryc, z, ut, wave, name(1:i)
1007	format (1x, f7.2, 1x, f6.2, 2(1x, f7.2), 1x, f6.3, 1x, 
     $     f7.4, 1x, f8.3, 1x, " 0", 1x, a)

c  extract ring profile
	axc     = 1068.
	ayc     = 1017.
	arad    = 958.
      arsq = arad * arad
      isample = 1.10 + 3.14159 * arsq  / (MAXPTS * NRING)


      do i=1,5
         flag(i) = .true.
      end do
      do i=1,NRING
         x(i) = (i - 0.5) / NRING
      end do
      
      do i=1,NRING
		isam(i) = 0
		numb(i) = 0
      end do

      rmax = arad - sqrt((axc - rxc)**2 + (ayc - ryc)**2)
      rmaxsq = rmax * rmax
      iylo = max(nint(ryc - rmax), 1)
      iyhi = min(nint(ryc + rmax), 2051)

      do iy=iylo,iyhi
		ry = iy - ryc
		rysq = ry * ry
		dx = sqrt(max(1.0,rmaxsq - rysq))
		ixlo = max (nint(rxc - dx), 1)
		ixhi = min (nint(rxc + dx), 2148)
		do ix=ixlo,ixhi
			if ((ix .le. 512) .or. 
     1			((ix .ge. 563) .and. (ix .le. 1586)) .or.
	2			(ix .ge. 1637)) then
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
         call biwgt (dat(1,i),numb(i),y(i),dy(i))
      end do


!  plot the profile

	pmin = y(1) - dy(1)
	pmax = y(1) + dy(1)
	do i=2,NRING
		pmin = min(pmin,(y(i)-dy(i)))
		pmax = max(pmax,(y(i)+dy(i)))
	end do	
	call pgsci(1)	
	call pgenv (0.0, 1.0, pmin, pmax, 0, 1)
	call pgline (NRING, x, y)

	xv(1) = 1.0 + NRING * Rave * rave / rmaxsq
	xv(2) = xv(1)
	yv(1) = pmin
	yv(2) = pmax
	call pgline (2, xv, yv)
	call pgupdt ()



101	continue



	write (*,1021)
1021	format ("Do another image? "$)
	read (*,1001) answ
	if (answ .eq. "y" .or. answ .eq. "Y") go to 100



c  ring completed - loop for more



c  done - close up shop

      close (10)
	call pgend ()

	stop
	end
