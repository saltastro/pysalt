      subroutine deletefile(filename,status)

C  A simple little routine to delete a FITS file

      integer status,unit,blocksize
      character*(*) filename

C  Simply return if status is greater than zero
      if (status .gt. 0)return


C  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)

C  Try to open the file, to see if it exists
      write(*, '(a)') filename
      call ftopen(unit,filename,1,blocksize,status)
      if (status .eq. 0) then
         print*, 'file exists already - deleting'
C         file was opened;  so now delete it 
          call ftdelt(unit,status)
      else if (status .eq. 103) then
c         print*, 'file doesnt exist'
C         file doesn't exist, so just reset status to zero and clear errors
          status=0
          call ftcmsg
c          call ftdelt(unit,status) ! nic added
      else
c         print*, 'some other error', status
C         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

C  Free the unit number for later reuse
      call ftfiou(unit, status)
      end
