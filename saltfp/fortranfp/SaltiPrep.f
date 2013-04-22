c *******************************************************************
c
c	SaltiPrep
c
c	This program prepares a clean single-extension FITS format image 
c	for later analysis.  Processing steps include:
c		1) Assemble the 4 amplifier segments into a single image.
c		2) Measure and subtract the overscan bias.
c
c
c *******************************************************************

		
	use dfport

	real over(100)
	integer idata(17211992), dat(4610648)
	integer naxes(2)
	integer npre(9), ndata(9), nover(9), ngap(9)
	integer nline(9), noskp(9)
	logical lsimple, ext
	character*80 name, indir, oudir, string, comment
	character*1 answ

c	geometry values - for original readout scheme

	data nline /4102, 2051, 1368, 1026,  821,  684,  586,  513,  456/
	data ndata /1024,  512,  341,  255,  204,  170,  145,  127,  113/
	data nover /  45,   23,   15,   11,    9,    8,    7,    6,    5/
	data noskp /   5,    2,    2,    2,    2,    1,    1,    1,    1/
	data npre  /  50,   25,   17,   13,   10,    9,    8,    7,    6/
	data ngap  / 100,   50,   34,   26,   22,   17,   15,   13,   12/


	write (*,998) 
998	format ("Enter input directory path name: "$)
	read (*,1001) indir
	istat = chdir(indir)
	if (istat .ne. 0) then
		write (*,999)
999		format ("Error changing directory")
		stop
	end if

	write (*,997) 
997	format ("Enter output directory path name: "$)
	read (*,1001) oudir
	istat = chdir(oudir)
	if (istat .ne. 0) then
		write (*,999)
		stop
	end if

100	continue

	write (*,1000) 
1000	format ("Enter input image filename: "$)
	read (*,1001) name
1001	format (a)
	istat = chdir(indir)
	call ftnopn (10, name, 0, istat)
	if (istat .ne. 0) then
		write (*,1011) istat
1011		format ("Error opening file, status = "i)
		call ftclos (10, istat)
		stop
	end if

c	get the binning factor
	call ftgkys (10, "CCDSUM", string, comment, istat)
	if (istat.ne. 0) then
		write (*,1060) istat
1060		format ("Error reading CCDSUM keyword, status = "i)
		call ftclos (10, istat)
		stop
	end if
	read (string, *) ibin, jbin
	if (ibin .ne. jbin) then
		write (*,1061) ibin, jbin
1061		format ("Error: unequal binning factors: ",i2,2x,i2)
		call ftclos (10, istat)
		stop
	end if
	if (ibin .gt. 9) then
	  write (*,1062) ibin
1062	  format ("Error: unsupported binning factor: "i2)
	  call ftclos (10, istat)
	  stop
	end if

c	calculate geometry
	ilinesize = npre(ibin) + ndata(ibin) + nover(ibin) + noskp(ibin)
	isegsize = nline(ibin) * ilinesize
	iskip = ndata(ibin) + npre(ibin) + noskp(ibin)
	jskip = npre(ibin) + nover(ibin) + noskp(ibin)
	kskip = 3 * ndata(ibin) + ngap(ibin)

	write (*,1002) 
1002	format ("Enter output image filename: "$)
	read (*,1001) name
	istat = chdir(oudir)
	call ftinit (11, name, 1, istat)
	if (istat .ne. 0) then
		write (*,1011) istat
		call ftclos (10, istat)
		stop
	end if
	call ftcphd (10, 11, istat)
	if (istat .ne. 0) then
		write (*,1012) istat
1012		format ("Error copying header, status = "i)
		call ftclos (10, istat)
		call ftclos (11, istat)
		stop
	end if


	do iseg = 1,4

c	read the file extent (2-5)
	  call ftmahd (10, iseg+1, itype, istat)	
	  if (istat .ne. 0) then
		write (*,1013) iseg, istat
1013		format ("Error moving to extent ",i1,", status = "i)
		call ftclos (10, istat)
		call ftclos (11, istat)
		stop
	  end if
	  call ftgpvj (10, 0, 1, isegsize, 0, dat, ext, istat)
	  if (istat .ne. 0) then
		write (*,1014) iseg, istat
1014		format ("Error reading extent ",i1,", status = "i)
		call ftclos (10, istat)
		call ftclos (11, istat)
		stop
	  end if

	  if (iseg .eq. 1) then
	    iloc = npre(ibin) + ndata(ibin) + noskp(ibin)
	    jloc = 1 + npre(ibin)
	    kloc = 1
	  else if (iseg .eq. 2) then
	    iloc = 1
	    jloc = 1 + nover(ibin) + noskp(ibin)
	    kloc = 1 + ndata(ibin)
	  else if (iseg .eq. 3) then
	    iloc = npre(ibin) + ndata(ibin) + noskp(ibin)
	    jloc = 1 + npre(ibin)
	    kloc = 1 + 2 * ndata(ibin) + ngap(ibin)
	  else if (iseg .eq. 4) then
	    iloc = 1
	    jloc = 1 + nover(ibin) + noskp(ibin)
	    kloc = 1 + 3 * ndata(ibin) + ngap(ibin)
	  endif

	  do j=1,nline(ibin)

c	calculate overscan
	    do i=1,nover(ibin)
		  over(i) = dat(iloc)
		  iloc = iloc + 1
		end do
	    call biwgt (over, nover(ibin), bias, sigma)
	    ibias = nint(bias)
 	    iloc = iloc + iskip
	
c	place data and subtract bias
	    do i=1,ndata(ibin)
	      idata(kloc) = dat(jloc) - ibias
	      jloc = jloc + 1
	      kloc = kloc + 1
	    end do
	    jloc = jloc + jskip
	    kloc = kloc + kskip

	  end do

	end do

c	mask unexposed areas
	loc = 1 + 2 * ndata(ibin)
	locinc = 4 * ndata(ibin)
	do j=1,nline(ibin)
	  do i=1,ngap(ibin)
	    idata(loc) = -666
	    loc = loc + 1
	  end do
	  loc = loc + locinc
	end do


	naxes(1) = 4*ndata(ibin) + ngap(ibin)
	naxes(2) = nline(ibin)
	ndatasize = naxes(1) * naxes(2)
	call ftrsim (11, 32, 2, naxes, istat)
	call ftpprj (11, 0, 1, ndatasize, idata(1), istat)
	call ftclos (11, istat)
	call ftclos (10, istat)

	write (*,1003)
1003	format ("Do another (y/n)? "$)
	read (*,1001) answ
	if (answ .eq. "Y" .or. answ .eq. "y") go to 100

	stop
	end
