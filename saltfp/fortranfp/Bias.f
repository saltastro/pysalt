c *******************************************************************
c
c	RSS_Prep
c
c	This program prepares a clean single-extension FITS format image 
c	for later analysis.  Processing steps include:
c		1) Assemble the 6 amplifier segments into a single image.
c		2) Measure and subtract the overscan bias.
c
c
c *******************************************************************
	program Bias

		
	use dfport


	real over(100), y(4102), sigma(4102), x(4102)
	real ylo(4102), yhi(4102), ys(4102), fit(20)  
	real xs(4102)
	integer dat(4610648)
	integer naxes(2)
	integer npre(9), ndata(9), nover(9), ngap(9)
	integer nline(9), noskp(9)
	logical lsimple, ext, trouble
	character*80 name, indir, oudir, string, comment
	character*1 answ

c	geometry values - for original readout scheme

	data nline /4102, 2051, 1368, 1026,  821,  684,  586,  513,  456/
	data ndata /1024,  512,  341,  255,  204,  170,  145,  127,  113/
	data nover /  45,   23,   15,   11,    9,    8,    7,    6,    5/
	data noskp /   5,    2,    2,    2,    2,    1,    1,    1,    1/
	data npre  /  50,   25,   17,   13,   10,    9,    8,    7,    6/
	data ngap  / 100,   50,   34,   26,   22,   17,   15,   13,   12/

	call pgbegin (0, "?", 1, 1)

	write (*,998) 
998	format ("Enter input directory path name: "$)
	read (*,1001) indir
	istat = chdir(indir)
	if (istat .ne. 0) then
		write (*,999)
999		format ("Error changing directory")
		stop
	end if

c	write (*,997) 
c997	format ("Enter output directory path name: "$)
c	read (*,1001) oudir
c	istat = chdir(oudir)
c	if (istat .ne. 0) then
c		write (*,999)
c		stop
c	end if

100	continue

	write (*,1000) 
1000	format ("Enter input image filename: "$)
	read (*,1001) name
1001	format (a)
c	istat = chdir(indir)
	call ftnopn (10, name, 0, istat)
	if (istat .ne. 0) then
		write (*,1011) istat
1011		format ("Error opening file, status = ",i)
		call ftclos (10, istat)
		stop
	end if

c	get the binning factor
	call ftgkys (10, "CCDSUM", string, comment, istat)
	if (istat.ne. 0) then
		write (*,1060) istat
1060		format ("Error reading CCDSUM keyword, status = ",i)
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
1062	  format ("Error: unsupported binning factor: ",i2)
	  call ftclos (10, istat)
	  stop
	end if
	write (*,1063) ibin
1063	format ("Binning Factor: ",i1)

c	calculate geometry
	ilinesize = npre(ibin) + ndata(ibin) + nover(ibin) + noskp(ibin)
	isegsize = nline(ibin) * ilinesize
	iskip = npre(ibin) + ndata(ibin) + noskp(ibin)
	jskip = npre(ibin) + nover(ibin) + noskp(ibin)
	kskip = 5 * ndata(ibin) + 2 * ngap(ibin)

c	write (*,1002) 
c1002	format ("Enter output image filename: "$)
c	read (*,1001) name
c	istat = chdir(oudir)
c	call ftinit (11, name, 1, istat)
c	if (istat .ne. 0) then
c		write (*,1011) istat
c		call ftclos (10, istat)
c		stop
c	end if
c	call ftcphd (10, 11, istat)
c	if (istat .ne. 0) then
c		write (*,1012) istat
c1012		format ("Error copying header, status = ",i)
c		call ftclos (10, istat)
c		call ftclos (11, istat)
c		stop
c	end if


	do iseg = 1,6

c	read the file extent
	  call ftmahd (10, iseg+1, itype, istat)	
	  if (istat .ne. 0) then
		write (*,1013) iseg, istat
1013		format ("Error moving to extent ",i1,", status = ",i)
		call ftclos (10, istat)
		call ftclos (11, istat)
		stop
	  end if
	  call ftgpvj (10, 0, 1, isegsize, 0, dat, ext, istat)
	  if (istat .ne. 0) then
		write (*,1014) iseg, istat
1014		format ("Error reading extent ",i1,", status = ",i)
		call ftclos (10, istat)
		call ftclos (11, istat)
		stop
	  end if

	  if (iseg .eq. 1) then
		iloc = 1 + npre(ibin) + ndata(ibin) + noskp(ibin)
	    jloc = 1 + npre(ibin)
	    kloc = 1
	  else if (iseg .eq. 2) then
	    iloc = 1
	    jloc = 1 + nover(ibin) + noskp(ibin)
	    kloc = 1 + ndata(ibin)
	  else if (iseg .eq. 3) then
		iloc = 1 + npre(ibin) + ndata(ibin) + noskp(ibin)
	    jloc = 1 + npre(ibin)
	    kloc = 1 + 2 * ndata(ibin) + ngap(ibin)
	  else if (iseg .eq. 4) then
	    iloc = 1
	    jloc = 1 + nover(ibin) + noskp(ibin)
	    kloc = 1 + 3 * ndata(ibin) + ngap(ibin)
	  else if (iseg .eq. 5) then
		iloc = 1 + npre(ibin) + ndata(ibin) + noskp(ibin)
	    jloc = 1 + npre(ibin)
	    kloc = 1 + 4 * ndata(ibin) + 2 * ngap(ibin)
	  else if (iseg .eq. 6) then
	    iloc = 1
	    jloc = 1 + nover(ibin) + noskp(ibin)
	    kloc = 1 + 5 * ndata(ibin) + 2 * ngap(ibin)
	  endif

	  do j=1,nline(ibin)

	    x(j) = j

c	calculate overscan
	    do i=1,nover(ibin)
		  over(i) = dat(iloc)
		  iloc = iloc + 1
		end do
	    call biwgt (over, nover(ibin), y(j), sigma(j))
c	    ibias = nint(bias)
 	    iloc = iloc + iskip
	
c	place data and subtract bias
c	    do i=1,ndata(ibin)
c	      idata(kloc) = dat(jloc) - ibias
c	      jloc = jloc + 1
c	      kloc = kloc + 1
c	    end do
c	    jloc = jloc + jskip
c	    kloc = kloc + kskip

	  end do

c	plot and fit bias vector here

	  xmin = x(1)
	  xmax = x(1)
	  ylo(1) = y(1) - sigma(1)
	  yhi(1) = y(1) + sigma(1)
	  ymin = ylo(1)
	  ymax = yhi(1)
	  do j=2,nline(ibin)
	    xmin = min (xmin, x(j))
	    xmax = max (xmax, x(j))
	    ylo(j) = y(j) - sigma(j)
	    yhi(j) = y(j) + sigma(j)
	    ymin = min (ymin, ylo(j))
	    ymax = max (ymax, yhi(j))
	  end do
	  xmin = xmin - 0.05 * (xmax - xmin)
	  xmax = xmax + 0.05 * (xmax - xmin)
	  ymin = ymin - 0.05 * (ymax - ymin)
	  ymax = ymax + 0.05 * (ymax - ymin)
	  call pgenv (xmin, xmax, ymin, ymax, 0, 1)
	  call pglabel ("row", "signal", "Bias")
c	  call pgpoint (nline(ibin), x, y, 0)
c	  call pgerry (nline(ibin), x, ylo, yhi, 1.0)
	  call pgline (nline(ibin), x, y)

c123	  do j=1,nline(ibin)
c	    ys(j) = y(j)
c	  end do

c	write (*,'("Enter box halfwidth: ",$)')
c	read (*,*) iw
c	if (iw .eq.0) go to 456

c	iw = nline(ibin)/50
c	sum = y(1)
c	num = 1
c	loc = 2
c	joc = 1
c
c	do j=1,iw
c	ys(j) = sum / num
c	sum = sum + y(loc)+ y(loc+1)
c	loc = loc + 2
c	num = num + 2
c	end do
c
c	do j=iw+1,nline(ibin)-iw-1
c	ys(j) = sum / num
c	sum = sum - y(joc) + y(loc)
c	joc = joc + 1
c	loc = loc + 1
c	end do
c
c	do j=nline(ibin)-iw,nline(ibin)-1
c	ys(j) = sum / num
c	sum = sum - y(joc) - y(joc+1)
c	joc = joc + 2
c	num = num - 2
c	end do
c	ys(nline(ibin)) = sum / num
c	write (*,*) num

c	ys(j) = ys(j) + y(j+k)
c	call lowpass (ys,2048,icut,2)

	write (*,'("Enter polynomial order: "$)')
	read (*,*) iord
	do i=1,nline(ibin)
		xs(i) = x(i) / nline(ibin)
	end do
	call polyfit (xs,y,nline(ibin), fit, iord, trouble)
	do i=1,nline(ibin)
		ys(i) = poly(xs(i), fit, iord)
	end do
			
	call pgsci(2)
	call pgline(nline(ibin),x,ys)
	call pgsci(1)
c	go to 123
c456	continue

	end do

c	mask unexposed areas
c	loc = 1
c	locinc = 2 * ndata(ibin)
c	do j=1,nline(ibin)
c	  loc = loc + locinc
c	  do i=1,ngap(ibin)
c	    idata(loc) = -666
c	    loc = loc + 1
c	  end do
c	  loc = loc + locinc
c	  do i=1,ngap(ibin)
c	    idata(loc) = -666
c	    loc = loc + 1
c	  end do
c	  loc = loc + locinc
c	end do


c	naxes(1) = 6 * ndata(ibin) + 2 * ngap(ibin)
c	naxes(2) = nline(ibin)
c	ndatasize = naxes(1) * naxes(2)
c	call ftrsim (11, 32, 2, naxes, istat)
c	call ftpprj (11, 0, 1, ndatasize, idata(1), istat)
c	call ftclos (11, istat)
	call ftclos (10, istat)

	write (*,1003)
1003	format ("Do another (y/n)? "$)
	read (*,1001) answ
	if (answ .eq. "Y" .or. answ .eq. "y") go to 100
	call pgend ()

	stop

	end program Bias
