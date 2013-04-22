      program delitt

	use dfport
	use dflib


      real data(6344,4102), bias(6344)
	integer isize(7)
	logical anyf
      character*80 biasfile, infile, outfile
      character*24 datetime
      character*1 bell



      write (*, '(/"RSS DeLitt")')


c  start the log file

      open (10, file="rssfp.log", access="append", iostat=ier)
      if (ier .ne. 0) then
         write (*, '(a1,"Error opening log file")') bell
         stop
      end if
      write (10, '(/"RSS DeLitt")')
      call fdate (datetime)
      write (10, '(a24)') datetime


c  get the bias image

 
      write (*, '("Enter name of bias image: "$)')
      read (*, '(a)') biasfile

      call ftnopn (11, biasfile, 0, ier)
      if (ier .ne. 0) then
         write (*, '(a1,"Error opening input image ",a)') 
     $        bell, biasfile(1:lnblnk(biasfile))
	   go to 666
      end if


	do i=1,nx
	   bias(i) = 0.0
	end do




c  get the output image or list
 
      write (*, '("Enter name of image to process: "$)')
      read (*, '(a)') infile


	do while (lnblnk(infile) .gt. 0)

         call ftnopn (12, infile, IRO, ier)
         if (ier .ne. 0) then
            write (*, '(a1,"Error opening input image ",a)') 
     $           bell, infile(1:lnblnk(infile))
 	      go to 666
         end if
	
	   call ftgipr (12, 7, bitpix, idim, isize, ier)
         if (ier .ne. 0) then
            write (*, '(a1,"Error reading input image header ",a)')  
     $           bell, infile(1:lnblnk(infile))
            call ftclos (12, ier)
	      go to 666
         end if
         nx = isize(1)
         ny = isize(2)
 
	   iloc = 1
         do line=1,ny
	      call ftgpve (12, 0, iloc, nx, 0, data(1,line), anyf, ier)
            if (ier .ne. 0) then
               write (*, '(a1,"Error reading image line ",i4)')
     $              bell, line
               go to 666
            end if
	      iloc = iloc + nx
         end do
  
	   call ftclos (12, ier)

	   write (*, '("Enter Littrow center x,y: ",$)')
	   read (*,*) xg, yg
	   write (*, '("Enter slit bottom and top: ",$)')
	   read (*,*) ibot, itop
	   lbot = 2.0 * yg - itop
	   ltop = 2.0 * yg - ibot


	   do iy=lbot,ltop
	      do ix = 1,10
	         back(ix) = data(xg-10-ix,iy)
	         back(ix+10) = data (xg +10+ix,iy)
	      end do
	      call biweight (back,20,gbk(iy),sig)
	      do ix=1,21
	         ghost(ix,iy) = data(xg+ix-11,iy) - gbk(iy)
	         data(xg+ix-11,iy) = gbk(iy)
	      end do
	   end do

	   do iy=ibot,itop
	      

	      

         write (*, '("Enter name of output image: "$)')
         read (*, '(a)') outfile
         call ftinit (13, outfile, 1, ier)
         if (ier .ne. 0) then
            write (*, '(a1,"Error creating output image ",a)') 
     $        bell, outfile(1:lnblnk(outfile))
            go to 666
         end if

         call ftcphd (12, 13, ier)
         if (ier .ne. 0) then
            write (*, '(a1,"Error writing output header ",a)')
     $        bell, outfile(1:lnblnk(outfile))
            go to 666
         end if

         call ftrsim (13, -32, 2, isize, ier)
         if (ier .ne. 0) then
            write (*, '(a1,"Error writing output header ",a)')
     $        bell, outfile(1:lnblnk(outfile))
            go to 666
         end if

         write (10, '("Output image: ",a)') outfile(1:lnblnk(outfile))


	      call ftppre (13, 0, iloc, nx, data, ier)
            if (ier .ne. 0) then
               write (*, '(a1,"Error writing output image line ",i4)')
     $           bell, line
               go to 666
            end if




	   call ftclos(12, ier)
	   call ftclos (13, ier)

	   write (10, '(a," --> ", a)') infile(1:lnblnk(infile)),
     $      outfile(1:lnblnk(outfile))

         write (*, '("Enter name of image to process: "$)')
         read (*, '(a)') infile

	end do

	write (10, '("Processing complete")')
	close (10)

	stop

666	write (10, '("Error Termination!")')
	close (10)
	call ftclos (12, ier)
	call ftclos (13, ier)
	stop

      end program delitt
 
