	Subroutine QuickCal (idata,ibin,   )

c *******************************************************************
c
c	QuickCal
c
c	This program measures the radius of a calibration ring in an RSS
c	FP image, for use in on-line quick-look wavelength zero point 
c	determination.  Processing steps include:
c		1) Assemble the central 4 amplifier segments into a single
c		   image with approximately correct geometry.
c		2) Measure and subtract the overscan bias.
c		3) Divide by the specified flat field.
c		4) Fit the ring.
c		5) Calculate the zero point.
c
c	This version of the code supports binning factors 1-9 in the 
c	original readout format.  Only square pixels are allowed.
c
c *******************************************************************

		

	real over(100)
	integer*2 idata(27663888)

	integer npre(9), ndata(9), nover(9), ngap1(9), ngap2(9)
	integer nline(9), noskp(9)
	logical flag
	character*80 flatname


c	geometry values - for original readout scheme
	data nline /4102, 2051, 1368, 1026,  821,  684,  586,  513,  456/
	data ndata /1024,  512,  341,  255,  204,  170,  145,  127,  113/
	data nover /  45,   23,   15,   11,    9,    8,    7,    6,    5/
	data noskp /   5,    2,    2,    2,    2,    1,    1,    1,    1/
	data npre  /  50,   25,   17,   13,   10,    9,    8,    7,    6/
	data ngap1 / 100,   50,   34,   25,    0,   18                  /
	data ngap2 / 100,   50,     ,     ,     ,     ,     ,     ,     / 

c	data npre   /   50,   25,   17,   12,    0,    9 /
c	data nover  /   50,   50,   50,   50,    0,   50 /
c	data ndata  / 1024,  512,  341,  256,    0,  170 /
c	data nyseg  / 4102, 2051, 1368, 1026,    0,  684 /
c	data ndata  / 1024,  512,  341,  256,    0,  170 /
c	data nblack /  100,   50,   34,   25,    0,   18 /
c	data numover / 45 / 

c	get information from auxilarily file

	open (20, file="/data/ccd/fabry/quicklook.dat", status=old, iostat=ierr)
	if (ierr .ne. 0) then
	  close (20)
	  ier = 1
	  return
	end if
	read (20,'(a)') flatname
	read (20,*) irmin, irmax
	read (20,*) wave
	read (20,*) CalA, CalB, CalF
	close (20)

	nseg = npre(ibin) + ndata(ibin) + nover(ibin) + noskp(ibin)

c	prepare the raw image

	iloc = 1
	jloc = 1

	do j=1,nline(ibin)

c	skip first segment
	  iloc = iloc + nseg

c	process second segment
	  do i=1,nover(ibin)
	    over(i) = idata(iloc)
	    iloc = iloc + 1
	  end do
	  call biwgt (over, nover(ibin), bias, sigma)
	  ibias = nint(bias)

	  iloc = iloc + noskp(ibin)
	  do i=1,ndata(ibin)
	    idata(jloc) = idata(iloc) - ibias
	    iloc = iloc + 1
	    jloc = jloc +1
	  end do

	  do i=1,ngap1(ibin)
	    idata(jloc) = -666
	    jloc = jloc + 1
	  end do

	  iloc = iloc + npre(ibin)

c	process third segment
	  kloc = iloc + ndata(ibin) + npre(ibin) + noskp(ibin)
	  do i=1,nover(ibin)
	    over(i) = idata(kloc)
	    kloc = kloc + 1
	  end do
	  call biwgt (over, nover(ibin), bias, sigma)
	  ibias = nint(ibias)

	  iloc = iloc + npre(ibin)
	  do i=1,ndata(ibin)
	    idata(jloc) = idata(iloc) - ibias
	    iloc = iloc + 1
	    jloc = jloc +1
	  end do

	  iloc = iloc + nover(ibin) + noskp(ibin)

c	process fourth segment
	  do i=1,nover(ibin)
	    over(i) = idata(iloc)
	    iloc = iloc + 1
	  end do
	  call biwgt (over, nover(ibin), bias, sigma)
	  ibias = nint(bias)

	  iloc = iloc + noskp(ibin)
	  do i=1,ndata(ibin)
	    idata(jloc) = idata(iloc) - ibias
	    iloc = iloc + 1
	    jloc = jloc +1
	  end do

	  do i=1,ngap2(ibin)
	    idata(jloc) = -666
	    jloc = jloc + 1
	  end do

	  iloc = iloc + npre(ibin)

c	process fifth segment
	  kloc = iloc + ndata(ibin) + npre(ibin) + noskp(ibin)
	  do i=1,nover(ibin)
	    over(i) = idata(kloc)
	    kloc = kloc + 1
	  end do
	  call biwgt (over, nover(ibin), bias, sigma)
	  ibias = nint(bias)
	  	
	  iloc = iloc + npre(ibin)
	  do i=1,ndata(ibin)
	    idata(jloc) = idata(iloc) - ibias
	    iloc = iloc + 1
	    jloc = jloc +1
	  end do

	  iloc = iloc + nover(ibin) + noskp(ibin)

c	skip sixth segment
	  iloc = iloc + nseg

	end do


c	divide by flat

c	if (flatname(1:lnblnk(flatname)) .ne. "null") then
c	  call ftopen (20, flatfile, 0, ibsize, istat)
c	read it line by line
c	    call ftgpvi (20, 0, iloc, npix, 0, array, flag, istat)
c	flatten the image
c	  call ftclos (20, istat)
c	end if

c	fit the ring
	call ringcen (...)

c	calculate the zero point

	CalA = wave * sqrt(1 - (rad/CalF)**2) - CalB * Z

	ier = 0

	return
	end
