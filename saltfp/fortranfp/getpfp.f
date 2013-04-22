c------------------------------------------------------------------------------
      subroutine getpfp (key, par)
c------------------------------------------------------------------------------

c  This routine retrieves entries from the PFP parameter file.  This file
c  is "pfp.par" in the current directory or specified by the environment
c  variable PFP_PAR.  If the file is located, it is searched for the
c  keyword string passed in KEY.  If the string is found, the remaining 
c  characters on its line are returned in the string PAR.  All failures
c  return an empty string.  All characters of a line following a "%" character
c  are considered to be comments and are ignored.  Blank lines are ignored.

c
c	use DFPORT

      character*80  line, filename
      character*(*) key, par


      par = ""

      call strlim (key, kbeg, kend)
      if (kend .le. 0) return


c  Open the parameter file

      call getlu (lu)
      if (lu .eq. 0) return

      open (lu, file="pfp.par", status="old", iostat=ier)
      if (ier .ne. 0) then
         call getenv ("PFP_PAR", filename)
         if (lnblnk(filename) .ne. 0) then
            open (lu, file=filename, status="old", iostat=ier)
         end if
      end if
      if (ier .ne. 0) return


c  Read and parse lines from the parameter file

 10   read (lu, '(a)', end=20) line

      call strlim (line, lbeg, lend)
      if (lend .le. 0) go to 10

      icom = index(line, "%")
      if (icom .eq. lbeg) go to 10
      if (icom .gt. 0) lend = icom - 1

      iend = index(line(lbeg:lend), " ")
      if (iend .le. 0) go to 10
      iend = lbeg + iend - 2
      if (line(lbeg:iend) .ne. key(kbeg:kend)) go to 10


c  Keyword found - return parameters

      iend = iend + 1
      call strlim (line(iend:lend), jbeg, jend)
      jbeg = iend + jbeg - 1
      jend = iend + jend - 1
      par = line(jbeg:jend)
      close (lu)
      return


c  Keyword not found

 20   close (lu)
      return
 
      end

c------------------------------------------------------------------------------
      subroutine yesno (string, flag, error)
c------------------------------------------------------------------------------

c  This routine interprets a string as a yes/no answer.  If the first 
c  non-blank character of the string is "y" or "Y" the flag is set true;
c  if it is "n" or "N", the flag is set false; in either case, the error
c  is set to 0.  If the string is null or starts with any other character,
c  the error is set to 1 and the flag is unchanged.


      integer       error
      logical       flag
      character*(*) string
      character*1   char

      call strlim (string, ibeg, iend)
      if (iend .eq. 0) then
         error = 1
      else
         char = string(ibeg:ibeg)
         if ((char .eq. "y") .or. (char .eq. "Y")) then
            flag = .true.
            error = 0
         else if ((char .eq. "n") .or. (char .eq. "N")) then
            flag = .false.
            error = 0
         else
            error = 1
         end if
      end if

      return
      end

c------------------------------------------------------------------------------
      subroutine strlim (string, ibeg, iend)
c------------------------------------------------------------------------------

c  This routine finds the first and last non-blank characters in a string,
c  returning their positions in ibeg and iend, respectively.  If the string
c  is null or contains only blanks, ibeg and iend are returned as 0.


      character*(*) string

      iend = lnblnk(string)
      if (iend .gt. 0) then
         ibeg = 1
         do while (string(ibeg:ibeg) .eq. " ")
            ibeg = ibeg + 1
         end do
      else
         ibeg = 0
      end if

      return
      end

c------------------------------------------------------------------------------
      subroutine getlu (lu)
c------------------------------------------------------------------------------

c  This routine returns the first available logical unit for Fortran I/O.
c  If there is no available unit in the range 20:99, LU is returned as 0.


      logical op

      lu = 20
      inquire (unit=lu, opened=op)

      do while (op .and. (lu .lt. 99))
         lu = lu + 1
         inquire (unit=lu, opened=op)
      end do

      if (op) lu = 0

      return
      end
