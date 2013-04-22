c *******************************************************************
c
c	FPprep
c
c	This program prepares a clean single-extension FITS format image 
c	for later Fabry-Perot analysis.  Processing steps include:
c		1) Assemble the central 4 amplifier segments into a single
c		   image.
c		2) Measure and subtract the overscan bias.
c		3) Divide by the specified flat field.
c		4) Mask non-illuminated region.
c
c	This version of the code does not implement pixel interpolation,
c	and thus only supports binning factors 1,2,3,4,and 6!
c
c *******************************************************************

		
	use dfport

	real over(100), bias(4102)
	integer idata(17622192), dat(4610648)
	integer naxes(2)
	integer npre(6), nxseg(6), nover(6), nyseg(6)
	integer ndata(6), nblack(6)
	logical lsimple, ext
	character*80 name, indir, oudir, string, comment
	character*1 answ

c	geometry values - these may need to be changed as readout scheme evolves
	data npre   /   50,   25,   17,   12,    0,    9 /
	data nover  /   50,   25,   12,    6,    0,    8 /
	data nxseg  / 1024,  512,  341,  256,    0,  170 /
	data nyseg  / 4102, 2051, 1368, 1026,    0,  684 /
	data ndata  / 1024,  512,  341,  256,    0,  170 /
	data nblack /  100,   50,   34,   25,    0,   18 /
c	data numover / 45 / 
	data numover / 23 /

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
1000	format ("Enter PFIS image filename: "$)
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
	if ((ibin .eq. 5) .or. (ibin .gt. 6)) then
	  write (*,1062) ibin
1062	  format ("Error: unsupported binning factor: "i)
	  call ftclos (10, istat)
	  stop
	end if

c	calculate geometry
	ixsegsize = npre(ibin) + nxseg(ibin) + nover(ibin)
	iysegsize = nyseg(ibin)
	isegsize = ixsegsize * iysegsize

	

c	get the image size
c	call ftgkys (10, "DETSIZE", string, comment, istat)
c	if (istat.ne. 0) then
c		write (*,1063) istat
c1063		format ("Error reading DETSIZE keyword, status = "i)
c		call ftclos (10, istat)
c		stop
c	end if
c	read (string, 1064) ixlo,ixhi,iylo,iyhi
c1064	format (4(1x,i))


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

	joc = 1
	do iseg = 1,4

c	read the file extent (3-6)
	  call ftmahd (10, iseg+2, itype, istat)	
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

c	calculate overscan vector

	  loc = 1
	  locinc = npre(ibin) + nxseg(ibin) + nover(ibin) - numover
	  if ((iseg .eq. 2) .or. (iseg .eq. 4)) then
	    loc = loc + locinc
	  end if
	  do i=1,iysegsize
	    do j=1,numover
		  over(j) = dat(loc)
		  loc = loc + 1
		end do
	    call biwgt (over, numover, bias(i), sigma)
 	    loc = loc + locinc
	  end do
c	future overscan vector smoothing here

c	place data and subtract bias

	  locinc = 3 * ndata(ibin) + 2 * nblack(ibin)
	  jocinc = npre(ibin) + nover(ibin)
	  if (iseg .eq. 1) then
	    loc = 1
	    joc = 1 + nover(ibin)
	  end if
	  if (iseg .eq. 2) then
	    loc = 1 + ndata(ibin) + nblack(ibin)
	    joc = 1 + npre(ibin)
	  end if
	  if (iseg .eq. 3) then
	    loc = 1 + 2 * ndata(ibin) + nblack(ibin) 
	    joc = 1 + nover(ibin)
	  end if
	  if (iseg .eq. 4) then
	    loc = 1 + 3 * ndata(ibin) + 2 * nblack(ibin)
	    joc = 1 + npre(ibin)
	  end if

	  do j=1,iysegsize
	    do k=1,nxseg(ibin)
	      idata(loc) = dat(joc) - bias(j)
	      loc = loc + 1
	      joc = joc + 1
	    end do
	    loc = loc + locinc
	    joc = joc + jocinc
	  end do
	end do

c	mask unexposed areas
	loc = 1 + ndata(ibin)
	locinc = 4 * ndata(ibin) + 2 * nblack(ibin) - nblack(ibin)
	joc = 2 * ndata (ibin) + nblack(ibin)
	do i=1,iysegsize
	  do j=1,nblack(ibin)
	    idata(loc) = -666
	    idata(loc+joc) = -666
	    loc = loc + 1
	  end do
	  loc = loc + locinc
	end do
c do circle mask here


c	divide by flat











	naxes(1) = 4*nxseg(ibin) + 2*nblack(ibin)
	naxes(2) = nyseg(ibin)
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
