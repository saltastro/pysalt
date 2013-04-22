c------------------------------------------------------------------------------
      program mask
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


      implicit none

      integer MAXSIZE, MAXSAMP, ISHORT, IREAL, IRW, IRO
      parameter (MAXSIZE=2048, MAXSAMP=5000)
      parameter	(ISHORT=3, IREAL=6, IRW=3, IRO=1)

      real         data(MAXSIZE), maskdata(MAXSAMP)
      real         axc, ayc, arad, radi, rado
      real         aradsq, aysq, rsq, value, risq, rosq, tmp
      integer      mask, nargs, img, jmg, loc, ier, isam, jsam, iout
      integer      line, ipix, ixlo, ixhi, iylo, iyhi, idim, itype
      integer      ibeg, iend, psbeg, psend, obeg, oend
      integer      lnblnk, rindex, iargc
      integer      isize(7)
      logical      ilist, olist, inopen, outopen, batch, flag
      character*80 infile, outfile, flist, param
      character*40 presuf
      character*24 datetime
      character*1  bell, answer


      bell    = char(7)
      inopen  = .false.
      outopen = .false.
      ilist   = .false.
      olist   = .false.

      write (*, '(/"Fabry-Perot Image Mask")')

      nargs = iargc ()
      batch = nargs .gt. 0
      if (nargs .gt. 2) then
         write (*, '(a,"Too many command line arguments: ",i2)') 
     $         bell, nargs
         stop
      end if


c  start the log file

      open (10, file="rfp.log", access="append", iostat=ier)
      if (ier .ne. 0) then
         write (*, '(a,"Error opening log file")') bell
         stop
      end if
      write (10, '(/"Fabry-Perot Image Mask")')
      call fdate (datetime)
      write (10, '(a24)') datetime


c  get the aperture data

      call getrfp ("axc", param)
      read (param, *, iostat=ier) axc
      do while (ier .ne. 0)
         write (*, '("Enter aperture x center: "$)')
         read (*, *, iostat=ier) axc
      end do

      call getrfp ("ayc", param)
      read (param, *, iostat=ier) ayc
      do while (ier .ne. 0)
         write (*, '("Enter aperture y center: "$)')
         read (*, *, iostat=ier) ayc
      end do

      call getrfp ("arad", param)
      read (param, *, iostat=ier) arad
      do while (ier .ne. 0)
         write (*, '("Enter aperture radius: "$)')
         read (*, *, iostat=ier) arad
      end do

      aradsq = arad * arad
      write (10, '("Aperture xc, yc, radius: ", 3f10.3)') 
     $      axc, ayc, arad


c  get the mask value scheme

      mask = 0

      call getrfp ("mask.value", param)
      read (param, *, iostat=ier) value
      if (ier .eq. 0) mask = 1

      call getrfp ("mask.region", param)
      read (param, *, iostat=ier) radi, rado
      if (ier .eq. 0) mask = 2

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
         rado = amin0 (rado, arad)
         if (radi .ge. rado) radi = rado - 5.0
         risq = radi * radi
         rosq = rado * rado
         isam = 1.1 + 3.14159 * (rosq - risq) / MAXSAMP
         write (10, '("Mask region: ", 2f10.3,i5)') radi, rado, isam
      end if


c  get input file list

 100  if (batch) then
         call getarg (1, infile)
      else
         infile = " "
         do while (lnblnk(infile) .le. 0)
            write (*, '("Input image or list: "$)')
            read (*, '(a)') infile
         end do
      end if

      call strlim (infile, ibeg, iend)
      if (infile(ibeg:ibeg) .eq. "@") then
         flist = infile(ibeg+1:iend)
         open (11, file=flist, status="old", iostat=ier)
         if (ier .ne. 0) then
            write (*, '(a,"Error opening input image list ",a)') 
     $           bell, flist(1:lnblnk(flist))
            if (batch) go to 400
            go to 100
         end if
         ilist = .true.
      else
         ilist = .false.
      end if


c  get the output method

 200  if (nargs .eq. 2) then
         call getarg (2, outfile)
         iout = 1
      else
         iout = -1
      end if

      do while (iout .lt. 0)
         write (*, '("Specify output method (new, prefix, suffix, ",
     $     "overwrite): "$)')
         read (*, '(a)') answer

         if ((answer .eq. "n") .or. (answer .eq. "N")) then
            outfile = " "
            do while (lnblnk(outfile) .le. 0)
               write (*, '("Output image or list: "$)')
               read (*, '(a)') outfile
            end do
            iout = 1

         else if ((answer .eq. "p") .or. (answer .eq. "P")) then
            presuf = " "
            do while (lnblnk(presuf) .le. 0)
               write (*, '("Enter prefix: "$)')
               read (*, '(a)') presuf
            end do
            call strlim (presuf, psbeg, psend)
            iout = 2
 
         else if ((answer .eq. "s") .or. (answer .eq. "S")) then
            presuf = " "
            do while (lnblnk(presuf) .le. 0)
               write (*, '("Enter suffix: "$)')
               read (*, '(a)') presuf
            end do
            call strlim (presuf, psbeg, psend)
            iout = 3

         else if ((answer .eq. "o") .or. (answer .eq. "O")) then
            iout = 4
            
         else
            write (*,'(a,"Please answer n, p, s, or o")')  bell
            
         end if
      end do
 
      if (iout .eq. 1) then
         call strlim (outfile, obeg, oend)
         if (outfile(obeg:obeg) .eq. "@") then
            flist = outfile(obeg+1:oend)
            open (12, file=flist, status="old", iostat=ier)
            if (ier .ne. 0) then
               write (*, '(a,"Error opening output image list ",a)')
     $               bell, flist(1:lnblnk(flist))
               if (batch) go to 400
               go to 200
            end if
            if (.not. ilist) then
               ier = 1
               write (*, '(a,"Error: input single, output list")')
     $              bell
               if (batch) go to 400
               go to 200
            end if
            olist = .true.
         else
            if (ilist) then
               ier = 1
               write (*, '(a,"Error: input list, output single")')
     $              bell
               if (batch) go to 400
               go to 200
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

      if (iout .eq. 4) then
         call  imopen (infile, IRW, img, ier)
      else
         call  imopen (infile, IRO, img, ier)
      end if
      if (ier .ne. 0) then
         write (*, '(a,"Error opening input image ",a)') 
     $          bell, infile(ibeg:iend)
         go to 400
      end if
      inopen = .true.

      call imgsiz (img, isize, idim, itype, ier)
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

      if ((itype .ne. ISHORT) .and. (itype .ne. IREAL)) then
         ier = 1
         write (*, '(a,"Wrong type image ",a)')  
     $          bell, infile(ibeg:iend)
         go to 400
      end if

      if (isize(1) .gt. MAXSIZE) then
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
         loc = rindex (infile, "/")
         if (loc .eq. 0) then
            outfile = presuf(psbeg:psend) // infile(ibeg:iend)
         else
            outfile = infile(ibeg:loc) // presuf(psbeg:psend) // 
     $           infile(loc+1:iend)
         end if

      else if (iout .eq. 3) then
         loc = rindex (infile, ".")
         if (loc .eq. 0) then
            outfile = infile(ibeg:iend) // presuf(psbeg:psend) 
         else
            outfile = infile(ibeg:loc-1) // presuf(psbeg:psend) // 
     $           infile(loc:iend)
         end if
      end if

      if (iout .eq. 4) then
         jmg = img
         outfile = infile
         obeg = ibeg
         oend = iend
      else
         call strlim (outfile, obeg, oend)
         call imcrea (outfile, isize, 2, IREAL, ier)
         if (ier .ne. 0) then
            write (*, '(a,"Error creating output image ",a)')
     $           bell, outfile(obeg:oend)
            go to 400
         end if
         
         call imopen (outfile, IRW, jmg, ier)
         if (ier .ne. 0) then
            write (*, '(a,"Error opening output image ",a)')
     $           bell, outfile(obeg:oend)
            go to 400
         end if

         outopen = .true.

         call imhcpy (img, jmg, ier)
         if (ier .ne. 0) then
            write (*, '(a,"Error writing output image header ",a)')
     $           bell, outfile(obeg:oend)
            go to 400
         end if
      end if


c  determine the mask value

      if (mask .eq. 2) then
         loc = 0
         jsam = 0
         ixlo = imax0(1,        nint(axc - rado - 1.0))
         ixhi = imin0(isize(1), nint(axc + rado + 1.0))
         iylo = imax0(1,        nint(ayc - rado - 1.0))
         iyhi = imin0(isize(2), nint(ayc + rado + 1.0))
         do line=iylo,iyhi
            aysq = (line - ayc)**2
            call imgl2r (img, data, line, ier)
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
                     maskdata(loc) = data(ipix)
                  end if
               end if
            end do
         end do

         call biwgt (maskdata, loc, value, tmp)
      end if

      call imakwr (jmg, "MASK", value, "image mask value", ier)
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

      do line=1,isize(2)
 
         aysq = (line - ayc)**2
 
         call imgl2r (img, data, line, ier)
         if (ier .ne. 0) then
            write (*, '(a,"Error reading input image, line ",i4)') 
     $           bell, line
            go to 400
         end if
  
         do ipix=1,isize(1)
            rsq = aysq + (ipix - axc)**2
            if (rsq .gt. aradsq) data(ipix) = value
         end do

         call impl2r (jmg, data, line, ier)
         if (ier .ne. 0) then
            write (*, '(a,"Error writing output image, line ",i4)')
     $           bell, line
            go to 400
         end if
 
      end do
 

c  done - close the image(s) and loop for more
 
      call imclos (img, ier)
      if (ier .ne. 0) then
         write (*, '(a,"Error closing input image")') bell
         go to 400
      end if
      inopen = .false.
 
      if (iout .ne. 4) then
         call imclos (jmg, ier)
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
         if (inopen)  call imclos (img,ier)
         if (outopen)  call imclos (jmg,ier)
      else
         write (*,  '("Successful completion")')
         write (10, '("Successful completion")')
      end if
      if (ilist) close (11)
      if (olist) close (12)

      if (.not. batch) then
         ier = 1
         do while (ier .ne. 0)
            write (*, '("Mask more images (y/n)? "$)')
            read (*, '(a)') answer
            call yesno (answer, flag, ier)
         end do
         if (flag) go to 100
      end if

      close (10)

      stop
      end
