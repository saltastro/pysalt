c------------------------------------------------------------------------------
      subroutine masking(logfile,axc, ayc,arad,mask,radi,rado,value,
     >     dir,infile,answer,outfile,prefix,suffix)
c------------------------------------------------------------------------------

c  This program masks the pixels outside the circular entrance aperture in
c  a Fabry-Perot image.  The mask value may be a specified constant or the
c  biweight of the pixel values in a specified annulus of the image.
c
c  The input image(s) are specified by name or in a list (list names begin
c  with the @ character); wildcards are not permitted.  The output image(s)
c  may replace the input image(s), or be new images with names specified
c  or derived by appending a prefix or a suffix to the input name.
c
c  Parameters specifying the aperture geometry and mask scheme may optionally
c  be specified in the rfp parameter file.  The actions of the program
c  are recorded in the "rfp.log" file.

c 12.05.10 replaced all references to imfort routines with fitsio routines
c replaced isize with naxes
c replaced calls to getrfp to getpfp
c 18.06.10 made compatible with pysalt package 

      implicit none

      integer MAXSIZE, MAXSAMP, ISHORT, IREAL, IRW, IRO
c      parameter (MAXSIZE=2048, MAXSAMP=5000)
c      upadated by nic (11.05.10)
      parameter (MAXSIZE=3172, MAXSAMP=5000)
      parameter	(ISHORT=3, IREAL=6, IRW=3, IRO=1)

      real         indata(MAXSIZE), maskdata(MAXSAMP)
      real         outdata(MAXSIZE*MAXSIZE)
      real         axc, ayc, arad, radi, rado
      real         aradsq, aysq, rsq, value, risq, rosq, tmp
      integer      mask, nargs, img, jmg, loc, ier, isam, jsam, iout
      integer      line, ipix, ixlo, ixhi, iylo, iyhi, idim, itype
      integer      ibeg, iend, psbeg, psend, obeg, oend, opix
      integer      lnblnk,cleariargc
      integer      isize(7), naxes(2), nkeys, nspace
      logical      ilist, olist, inopen, outopen, batch, flag
      character*40 infile, outfile, flist, param,logfile
      character*40 presuf, prefix, suffix, dir
      character*24 datetime, comment
      character*1  bell, answer
      integer      blocksize, readwrite, inunit, outunit, keyval
      integer      group, nullvalj, npixels, firstpix, bitpix
      logical      anynull, simple, extend
      real         nullvale
      integer      nelements, naxis


      bell    = char(7)
      inopen  = .false.
      outopen = .false.
      ilist   = .false.
      olist   = .false.

      write (*, '(/"Fabry-Perot Image Mask")')

c      nargs = iargc ()
c      batch = nargs .gt. 0
c      if (nargs .gt. 2) then
c         write (*, '(a,"Too many command line arguments: ",i2)') 
c     $         bell, nargs
c         stop
c      end if


c nsl added, to change to dir containing data and record all logs there

      ier = 1
      if (lnblnk(dir) .gt. 0) then
         ier = chdir(dir)
         if (ier .ne. 0) write (*, '("Error changing directory")')
      else
         ier = 0
      end if


c  start the log file

      open (10, file=logfile, access="append", iostat=ier)
      if (ier .ne. 0) then
         write (*, '(a,"Error opening log file")') bell
         stop
      end if
      write (10, '(/"Fabry-Perot Image Mask")')
      call fdate (datetime)
      write (10, '(a24)') datetime


c  get the aperture data

c      call getpfp ("axc", param)
c      read (param, *, iostat=ier) axc
c      do while (ier .ne. 0)
c         write (*, '("Enter aperture x center: "$)')
c         read (*, *, iostat=ier) axc
c      end do

c      call getpfp ("ayc", param)
c      read (param, *, iostat=ier) ayc
c      do while (ier .ne. 0)
c         write (*, '("Enter aperture y center: "$)')
c         read (*, *, iostat=ier) ayc
c      end do

c      call getpfp ("arad", param)
c      read (param, *, iostat=ier) arad
c      do while (ier .ne. 0)
c         write (*, '("Enter aperture radius: "$)')
c         read (*, *, iostat=ier) arad
c      end do

      aradsq = arad * arad
      write (10, '("Aperture xc, yc, radius: ", 3f10.3)') 
     $      axc, ayc, arad


c  get the mask value scheme

c      mask = 0

c      call getpfp ("mask.value", param)
c      read (param, *, iostat=ier) value
c      if (ier .eq. 0) mask = 1

c      call getpfp ("mask.region", param)
c      read (param, *, iostat=ier) radi, rado
c      if (ier .eq. 0) mask = 2

      do while (mask .eq. 0)
         write (*, '("Mask method: constant (c) or region (r)? "$)')
         read (*, '(a)') answer
         if ((answer .eq. "c") .or. (answer .eq. "C")) then
            ier = 1
            do while (ier .ne. 0)
               write (*, '("Enter mask constant: "$)')
               read (*, *, iostat=ier) value
            end do
            mask = 1
         else if ((answer .eq. "r") .or. (answer .eq. "R")) then
            ier = 1
            do while (ier .ne. 0)
               write (*, '("Enter inner and outer radii: "$)')
               read (*, *, iostat=ier) radi, rado
            end do
            mask = 2
         else
            write (*, '(a,"Please answer c or r")') bell
         end if
      end do

      if (mask .eq. 1) then
         write (10, '("Mask value: ", f10.3)') value
      else
         radi = abs (radi)
         rado = abs (rado)
         if (radi .gt. rado) then
            tmp = radi
            radi = rado
            rado = tmp
         end if
         rado = min(rado, arad)
         if (radi .ge. rado) radi = rado - 5.0
         risq = radi * radi
         rosq = rado * rado
         isam = 1.1 + 3.14159 * (rosq - risq) / MAXSAMP
         write (10, '("Mask region: ", 2f10.3,i5)') radi, rado, isam
      end if


c  get input file list

c 100  if (batch) then
c         call getarg (1, infile)
c      else
c         infile = " "
c         do while (lnblnk(infile) .le. 0)
c            write (*, '("Input image or list: "$)')
c            read (*, '(a)') infile
c         end do
c      end if

      call strlim (infile, ibeg, iend)
      if (infile(ibeg:ibeg) .eq. "@") then
         flist = infile(ibeg+1:iend)
         open (11, file=flist, status="old", iostat=ier)
         if (ier .ne. 0) then
            write (*, '(a,"Error opening input image list ",a)') 
     $           bell, flist(1:lnblnk(flist))

c            if (batch) go to 400
c            go to 100
         end if
         ilist = .true.
      else
         ilist = .false.
      end if


c  get the output method

c 200  if (nargs .eq. 2) then
c         call getarg (2, outfile)
c         iout = 1
c      else
c         iout = -1
c      end if

c      do while (iout .lt. 0)
c         write (*, '("Specify output method (new, prefix, suffix ",
c     $     "): " $)')
c         read (*, '(a)') answer


         if ((answer .eq. "n") .or. (answer .eq. "N")) then
c            outfile = " "
c            do while (lnblnk(outfile) .le. 0)
c               write (*, '("Output image or list: "$)')
c               read (*, '(a)') outfile
c            end do
            iout = 1

         else if ((answer .eq. "p") .or. (answer .eq. "P")) then
            presuf = " "
c            do while (lnblnk(presuf) .le. 0)
c               write (*, '("Enter prefix: "$)')
c               read (*, '(a)') presuf
c            end do
            presuf = prefix
            call strlim (presuf, psbeg, psend)
            iout = 2
 
         else if ((answer .eq. "s") .or. (answer .eq. "S")) then
            presuf = " "
c            do while (lnblnk(presuf) .le. 0)
c               write (*, '("Enter suffix: "$)')
c               read (*, '(a)') presuf
c            end do
            presuf = suffix
            call strlim (presuf, psbeg, psend)
            iout = 3

c         else if ((answer .eq. "o") .or. (answer .eq. "O")) then
c            iout = 4
            
         else
            write (*,'(a,"Please answer n, p, or s")')  bell
            
         end if
         
c      end do
 
      if (iout .eq. 1) then
         call strlim (outfile, obeg, oend)
         if (outfile(obeg:obeg) .eq. "@") then
            flist = outfile(obeg+1:oend)
            open (12, file=flist, status="old", iostat=ier)
            if (ier .ne. 0) then
               write (*, '(a,"Error opening output image list ",a)')
     $               bell, flist(1:lnblnk(flist))
c               if (batch) go to 400
c               go to 200
            end if
            if (.not. ilist) then
               ier = 1
               write (*, '(a,"Error: input single, output list")')
     $              bell
c               if (batch) go to 400
c               go to 200
            end if
            olist = .true.
         else
            if (ilist) then
               ier = 1
               write (*, '(a,"Error: input list, output single")')
     $              bell
c               if (batch) go to 400
c               go to 200
            end if
            olist = .false.
         end if
      end if

c  processing loop

c  open the input image

 300  if (ilist) then
         read (11, '(a)', end=400) infile
         call strlim (infile, ibeg, iend)
      end if

c converting to fitsio from iraf routines
C  The STATUS (ier) parameter must always be initialized.
      ier=0
      blocksize=1
C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(inunit,ier)

      if (iout .eq. 4) then
c         call  imopen (infile, IRW, img, ier)
         readwrite=1            ! read/write access 
      else
c         call  imopen (infile, IRO, img, ier)
         readwrite=0            ! read access 
      end if
      call ftopen(inunit,infile,readwrite,blocksize,ier)

      if (ier .ne. 0) then
         write (*, '(a,"Error opening input image ",a)') 
     $          bell, infile(ibeg:iend)
         go to 400
      end if
      inopen = .true.

      call ftgknj(inunit,'NAXIS',1,2,naxes,idim,ier)
c      call imgsiz (img, isize, idim, itype, ier)


      if (ier .ne. 0) then
         write (*, '(a,"Error reading input image header ",a)')
     $          bell, infile(ibeg:iend)
         go to 400
      end if

      if (idim .ne. 2) then
         ier = 1
         write (*, '(a,"Wrong dimension image ",a)')
     $          bell, infile(ibeg:iend)
         go to 400
      end if

      call FTGKYJ(inunit,'BITPIX', keyval,comment,ier)
      bitpix=keyval

c      if ((itype .ne. ISHORT) .and. (itype .ne. IREAL)) then
      if ((bitpix .ne. 16) .and. (bitpix .ne. -32)) then
         ier = 1
         write (*, '(a,"Wrong type image ",a)')  
     $          bell, infile(ibeg:iend)
         go to 400
      end if

c isize(1) before
c      print*, naxes(1), maxsize
      if (naxes(1) .gt. MAXSIZE) then
         ier = 1
         write (*, '(a,"Image too large ",a)')
     $          bell, infile(ibeg:iend)
         go to 400
      end if


c  open the output image

      if (olist) then
         ier = 1
         read (12, '(a)', end=400) outfile

      else if (iout .eq. 2) then
         loc = index (infile, "/")
         if (loc .eq. 0) then
            outfile = presuf(psbeg:psend) // infile(ibeg:iend)
         else
            outfile = infile(ibeg:loc) // presuf(psbeg:psend) // 
     $           infile(loc+1:iend)
         end if

      else if (iout .eq. 3) then
         loc = index (infile, ".")
         if (loc .eq. 0) then
            outfile = infile(ibeg:iend) // presuf(psbeg:psend) 
         else
            outfile = infile(ibeg:loc-1) // presuf(psbeg:psend) // 
     $           infile(loc:iend)
         end if
      end if

c      if (iout .eq. 4) then
c         jmg = img
c         outfile = infile
c         obeg = ibeg
c         oend = iend
c      else

      call strlim (outfile, obeg, oend)
c         call imcrea (outfile, isize, 2, IREAL, ier)

C  Delete the file if it already exists, so we can then recreate it.
C  The deletefile subroutine is listed at the end of this file.
      ier=0         
      call deletefile(outfile,ier)
      
c Get an unused logical unit number
      ier=0
      call ftgiou(outunit,ier)

C  Create the new empty FITS file.  The blocksize parameter is a
C  historical artifact and the value is ignored by FITSIO.

      blocksize=1
      ier=0
      
      call ftinit(outunit,outfile,blocksize,ier)
      
      if (ier .ne. 0) then
         write (*, '(a,"Error creating output image ",a)')
     $        bell, outfile(obeg:oend)
         go to 400
      end if
             
c         call imopen (outfile, IRW, jmg, ier)
c         if (ier .ne. 0) then
c            write (*, '(a,"Error opening output image ",a)')
c     $           bell, outfile(obeg:oend)
c            go to 400
c         end if

      outopen = .true.

c Copy header from input image to output image.
c         call imhcpy (img, jmg, ier)
      call FTCPHD(inunit, outunit, ier)

      if (ier .ne. 0) then
         write (*, '(a,"Error writing output image header ",a)')
     $        bell, outfile(obeg:oend)
         go to 400
      end if
c     end if
     
c  determine the mask value

C  Initialize variables for image reading

      npixels=naxes(1) ! no in x direction
      firstpix=1
      nullvalj=-999
      nullvale=-999.0
      group=1

      if (mask .eq. 2) then
         loc = 0
         jsam = 0
         ixlo = max(1,        nint(axc - rado - 1.0))
         ixhi = min(naxes(1), nint(axc + rado + 1.0))
         iylo = max(1,        nint(ayc - rado - 1.0))
         iyhi = min(naxes(2), nint(ayc + rado + 1.0))
         do line=iylo,iyhi
            ier=0
            aysq = (line - ayc)**2
c            call imgl2r (img, indata, line, ier)
            firstpix= 1 + (line-1)*(npixels)
            call ftgpve(inunit,group,firstpix,npixels,nullvale,
     &           indata,anynull,ier)
            if (ier .ne. 0) then
               write (*, '(a,"Error reading input image, line ",i4)')
     $              bell, line
               go to 400
            end if
  
            do ipix=ixlo,ixhi
               rsq = aysq + (ipix - axc)**2
               if ((rsq .ge. risq) .and. (rsq .le. rosq)) then
                  jsam = jsam + 1
                  if (mod(jsam,isam) .eq. 0) then
                     loc = loc + 1
                     maskdata(loc) = indata(ipix)
                  end if
               end if
            end do
         end do
         call biwgt (maskdata, loc, value, tmp)
      end if

c      call imakwr (jmg, "MASK", value, "image mask value", ier)
      ier=0
      call ftpkye(outunit,'MASK',value,2,'Image mask value',ier)

      if (ier .ne. 0) then
         write (*, '("Error writing output image header ",a)')
     $        bell, outfile(obeg:oend)
         go to 400
      end if

      write (*, '(a," --> ",a,5x,f10.3)')
     $       infile(ibeg:iend), outfile(obeg:oend), value

      write (10, '(a," --> ",a,5x,f10.3)')
     $       infile(ibeg:iend), outfile(obeg:oend), value

c  mask the image

      npixels=naxes(1) ! no. in x direction
      firstpix=1
      nullvalj=-999
      nullvale=-999.0
      group=1
      opix=0

      do line=1,naxes(2)
         aysq = (line - ayc)**2
         firstpix= 1 + (line-1)*(npixels)
         if (bitpix .eq. 16) then
            call ftgpvj(inunit,group,firstpix,npixels,nullvalj,
     &           indata,anynull,ier)
         else if (bitpix .eq. -32) then
            call ftgpve(inunit,group,firstpix,npixels,nullvale,
     &           indata,anynull,ier)            
         endif

c         call imgl2r (img, indata, line, ier)

         if (ier .ne. 0) then
            write (*, '(a,"Error reading input image, line ",i4)') 
     $           bell, line
            go to 400
         end if
  
         do ipix=1,naxes(1)
            rsq = aysq + (ipix - axc)**2
            if (rsq .gt. aradsq) indata(ipix) = value
            opix=  ipix + (line-1)*(npixels)
            outdata(opix)= indata(ipix)
         end do

c         call impl2r (jmg, indata, line, ier)
         group=1
         firstpix=1
         nelements=naxes(1)*naxes(2)
         
         call ftppre(outunit,group,firstpix,nelements,outdata,ier)

         if (ier .ne. 0) then
            write (*, '(a)') "Error writing output image"
            go to 400
         end if
 
      end do

c  done - close the image(s) and loop for more
 
      call ftclos (inunit, ier)
      call ftfiou (inunit, ier)
      if (ier .ne. 0) then
         write (*, '(a,"Error closing input image")') bell
         go to 400
      end if
      inopen = .false.
 
      if (iout .ne. 4) then
         call ftclos (outunit, ier)
         call ftfiou (outunit, ier)
         if (ier .ne. 0) then
            write (*, '(a,"Error closing output image")') bell
            go to 400
         end if
      outopen = .false.
      end if

      if (ilist) go to 300


c  all done

 400  if (ier .ne. 0) then
         write (*,  '("Error termination")')
         write (10, '("Error termination")')
         if (inopen)  call ftclos (inunit,ier)
         if (outopen)  call ftclos (outunit,ier)
      else
         write (*,  '("Successful completion")')
         write (10, '("Successful completion")')
      end if
      if (ilist) close (11)
      if (olist) close (12)

c      if (.not. batch) then
c         ier = 1
c         do while (ier .ne. 0)
c            write (*, '("Mask more images (y/n)? "$)')
c            read (*, '(a)') answer
c            call yesno (answer, flag, ier)
c         end do
c         if (flag) go to 100
c      end if

      close (10)

c     stop
      return
      end subroutine masking
