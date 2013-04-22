c------------------------------------------------------------------------------
      subroutine getrfp (key, par)
c------------------------------------------------------------------------------

c  This routine retrieves entries from the RFP parameter file.  This file
c  is "rfp.par" in the current directory or specified by the environment
c  variable RFP_PAR.  If the file is located, it is searched for the
c  keyword string passed in KEY.  If the string is found, the remaining 
c  characters on its line are returned in the string PAR.  All failures
c  return an empty string.  All characters of a line following a "%" character
c  are considered to be comments and are ignored.  Blank lines are ignored.


      implicit none

      integer       lu, ier, icom, kbeg, kend, lbeg, lend
      integer       jbeg, jend, iend, lnblnk
      logical       op
      character*80  line, filename
      character*(*) key, par


      par = ""

      call strlim (key, kbeg, kend)
      if (kend .le. 0) return


c  Open the parameter file

      lu = 20
      inquire (unit=lu, opened=op)
      do while (op)
         lu = lu + 1
         inquire (unit=lu, opened=op)
      end do
      open (lu, file="rfp.par", status="old", iostat=ier)
      if (ier .ne. 0) then
         call getenv ("RFP_PAR", filename)
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

      implicit none

      integer       error, ibeg, iend
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

      implicit none

      integer ibeg, iend, lnblnk
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
