c *******************************************************************
c
c	RSS_Prep
c
c	This program prepares a clean single-extension FITS format image 
c	for later analysis.  Processing steps include:
c		1) Assemble the 6 amplifier segments into a single image.
c		2) Measure and subtract the overscan bias.
c	The overscan is fit with polynomial to smooth the noise.  Fitting
c	can be done in slit mode, which fits three segments (vertically),
c	of non-slit mode, which fits the entire overscan vector.
c	The output image is converted to real pixel type.
c
c *******************************************************************
	subroutine RSS_prep(answ,indir,listname,prefix, plottype)

		
c	use dfport
c	use dflib

	real over(100)
	integer dat(4610648)
	real idata(26023088)
	integer naxes(2), statusr
	integer npre(9), ndata(9), nover(9), ngap(9)
	integer nline(9), noskp(9)
	real x(4102), xs(4102), y(4102), ys(4102), ysig(4102)
	real fit1(20), fit2(20), fit3(20)
	logical lsimple, ext, trouble, slit
	character*80 name,indir,oudir,string,comment,listname,outname
	character*1 answ
	character*4 plottype
	integer ibeg, iend
	character*10 prefix

c	geometry values - for original readout scheme

	data nline /4102, 2051, 1368, 1026,  821,  684,  586,  513,  456/
	data ndata /1024,  512,  341,  255,  204,  170,  145,  127,  113/
	data nover /  45,   23,   15,   11,    9,    8,    7,    6,    5/
	data noskp /   5,    2,    2,    2,    2,    1,    1,    1,    1/
	data npre  /  50,   25,   17,   13,   10,    9,    8,    7,    6/
	data ngap  / 100,   50,   34,   26,   22,   17,   15,   13,   12/

	call pgbegin (0, plottype, 2, 3)
c	call pgbegin (0, "?", 2, 3)
	call pgask (.false.)

	slit = .false.
c	write (*,996)
c996	format ("Slit mode (y/n)? " )
c	read (*,1001) answ
	if (answ .eq. "y" .or. answ .eq. "Y") slit = .true.

c	write (*,998) 
c998	format ("Enter input directory path name: " )
c, 	read (*,1001) indir
	istat = chdir(indir)
	if (istat .ne. 0) then
		write (*,999)
999		format ("Error changing directory")
		stop
	end if

c	write (*,997) 
997	format ("Enter output directory path name: " )
c	read (*,1001) oudir
c	istat = chdir(oudir)
c	if (istat .ne. 0) then
c		write (*,999)
c		stop
c	end if

100	continue


c	write (*,1000) 
c1000	format ("Enter input filelist  filename: ")
c	read (*,1001) listname
1001	format (a)
c	write (*,1002) 
c1002	format ("Enter prefix for reduced files: ")
c	read (*,1001) prefix

	open(1, file=listname, status='unknown')
	statusr = 0
	name=' '
	do while (statusr .eq. 0)
	   read(1, 1001, iostat=statusr) name
	   if (statusr .eq. 0) then
	      istat=0

c	istat = chdir(indir)
	      call ftnopn (10, name, 0, istat)
	      if (istat .ne. 0) then
		 write (*,1011) istat
 1011		 format ("Error opening file status = ",i3)
		 call ftclos (10, istat)
		 stop
	      end if

c	get the binning factor
	      call ftgkys (10, "CCDSUM", string, comment, istat)
	      if (istat.ne. 0) then
		 write (*,1060) istat
 1060		 format ("Error reading CCDSUM keyword status = ",i3)
		 call ftclos (10, istat)
		 stop
	      end if
	      read (string, *) ibin, jbin
	      if (ibin .ne. jbin) then
		 write (*,1061) ibin, jbin
 1061		 format ("Error: unequal binning factors: ",i2,2x,i2)
		 call ftclos (10, istat)
		 stop
	      end if
	      if (ibin .gt. 9) then
		 write (*,1062) ibin
 1062		 format ("Error: unsupported binning factor: ",i2)
		 call ftclos (10, istat)
		 stop
	      end if
	      write (*,1063) ibin
 1063	      format ("Binning Factor: ",i1)
	      
c	calculate geometry
	      ilinesize = npre(ibin)+ndata(ibin)+nover(ibin)+noskp(ibin)
	      isegsize = nline(ibin) * ilinesize
	      iskip = npre(ibin) + ndata(ibin) + noskp(ibin)
	      jskip = npre(ibin) + nover(ibin) + noskp(ibin)
	      kskip = 5 * ndata(ibin) + 2 * ngap(ibin)
	      
	      if (slit) then
		 write (*,1070)
 1070		 format ("Enter slit bottom and top and poly order: ")
		 read (*,*) ibot, itop, iord
		 iov1 = 1
		 nov1 = ibot - iov1
		 iord1 = iord
		 iov2 = ibot
		 nov2 = itop - ibot + 1
		 iord2 = iord
		 iov3 = itop + 1
		 nov3 = nline(ibin) - itop
		 iord3 = iord
	      else
		 iov1 = 1
		 nov1 = nline(ibin)
		 iord1 = 5
	      end if
	      
c	      write (*,1002) 
c 1002	      format ("Enter output image filename: ")
c	      read (*,1001) name
c	istat = chdir(oudir)
	      call strlim (prefix, ibeg, iend)
	      outname(ibeg:iend) = prefix(ibeg:iend) 
	      outname(iend+1:) = name
	      call ftinit (11, outname, 1, istat)
	      if (istat .ne. 0) then
		 write (*,1011) istat
		 call ftclos (10, istat)
		 stop
	      end if
	      call ftcphd (10, 11, istat)
	      if (istat .ne. 0) then
		 write (*,1012) istat
 1012		 format ("Error copying header status = ",i3)
		 call ftclos (10, istat)
		 call ftclos (11, istat)
		 stop
	      end if
	      
	      
	      do iseg = 1,6
		 
c	read the file extent
		 call ftmahd (10, iseg+1, itype, istat)	
		 if (istat .ne. 0) then
		    write (*,1013) iseg, istat
 1013		    format ("Error moving to extent ",i1,"status = ",i3)
		    call ftclos (10, istat)
		    call ftclos (11, istat)
		    stop
		 end if
		 call ftgpvj (10, 0, 1, isegsize, 0, dat, ext, istat)
		 if (istat .ne. 0) then
		    write (*,1014) iseg, istat
 1014		    format ("Error reading extent ",i1,"status = ",i3)
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
		    xs(j) = x(j) / nline(ibin)
c	calculate overscan
		    do i=1,nover(ibin)
		       over(i) = dat(iloc)
		       iloc = iloc + 1
		    end do
		    call biwgt (over, nover(ibin), y(j), ysig(j))
		    iloc = iloc + iskip
		 end do
		 
		 call polyfit(xs(iov1), y(iov1), nov1, fit1, iord1, trouble)
		 if (slit) then
		    call polyfit(xs(iov2), y(iov2), nov2, fit2, iord2, trouble)
		    call polyfit(xs(iov3), y(iov3), nov3, fit3, iord3, trouble)
		 end if
		 
		 do j=iov1,nov1
		    ys(j) = poly(xs(j), fit1, iord1)
		 end do
		 if (slit) then
		    do j=iov2,iov2+nov2+1
		       ys(j) = poly(xs(j), fit2, iord2)
		    end do
		    do j=iov3,iov3+nov3-1
		       ys(j) = poly(xs(j), fit3, iord3)
		    end do
		 end if
		 
		 xmin = x(1)
		 xmax = x(1)
		 ymin = y(1)
		 ymax = y(1)
		 do j=2,nline(ibin)
		    xmin = min (xmin, x(j))
		    xmax = max (xmax, x(j))
		    ymin = min (ymin, y(j))
		    ymin = min (ymin, ys(j))
		    ymax = max (ymax, y(j))
		    ymax = max (ymax, ys(j))
		 end do
		 ymin = ymin - 0.5 * (ymax - ymin)
		 ymax = ymax + 0.5 * (ymax - ymin)
		 call pgenv(xmin,xmax,ymin,ymax,0,1)
		 call pglabel ("row", "signal", "Bias")
		 call pgline (nline(ibin), x, y)
		 call pgsci (2)
		 call pgline (nline(ibin), x, ys)
		 call pgsci (1)
		 
		 do j=1,nline(ibin)
		    
c	place data and subtract bias
		    do i=1,ndata(ibin)
		       idata(kloc) = dat(jloc) - ys(j)
		       jloc = jloc + 1
		       kloc = kloc + 1
		    end do
		    jloc = jloc + jskip
		    kloc = kloc + kskip
		    
		 end do
		 
	      end do
	      
c	mask unexposed areas
	      loc = 1
	      locinc = 2 * ndata(ibin)
	      do j=1,nline(ibin)
		 loc = loc + locinc
		 do i=1,ngap(ibin)
		    idata(loc) = -666
		    loc = loc + 1
		 end do
		 loc = loc + locinc
		 do i=1,ngap(ibin)
		    idata(loc) = -666
		    loc = loc + 1
		 end do
		 loc = loc + locinc
	      end do
	      
	      
	      naxes(1) = 6 * ndata(ibin) + 2 * ngap(ibin)
	      naxes(2) = nline(ibin)
	      ndatasize = naxes(1) * naxes(2)
	      call ftrsim (11, -32, 2, naxes, istat)
	      call ftppre (11, 0, 1, ndatasize, idata(1), istat)
	      call ftclos (11, istat)
	      call ftclos (10, istat)
	   endif
	enddo
	close(1)
c	      write (*,1003)
c 1003	      format ("Do another (y/n)? ")
c	      read (*,1001) answ
c	      if (answ .eq. "Y" .or. answ .eq. "y") go to 100
	 
	write (*, '("Processing completed")')
	call pgend ()

	end subroutine RSS_prep
	
