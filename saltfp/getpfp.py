#!/usr/bin/env python

# Author                 Version      Date
# -----------------------------------------------
# Nicola Loaring (SAAO)    1.0          2 June 2009

# getpfp.py, converted from Ted's Fortran getpfp.f
# this is just for reading the parameter files

#  This routine retrieves entries from the PFP parameter file.  This file
#  is "pfp.par" in the current directory or specified by the environment
#  variable PFP_PAR.  If the file is located, it is searched for the
#  keyword string passed in KEY.  If the string is found, the remaining 
#  characters on its line are returned in the string PAR.  All failures
#  return an empty string.  All characters of a line following a "%" character
#  are considered to be comments and are ignored.  Blank lines are ignored.


#      character*80  line, filename
#      character*(*) key, par

def getpfp (filename, key):
 par = ""
 flag = 0
# strlim (key, kbeg, kend)
# if (kend <= 0) return


#  Open the parameter file, create a file object

 infile= open(filename, "r")
 line="dummy start"

 while line != "":
     line= infile.readline()
     if line.find(key) != -1:
         cols=line.split()
         nocols= len(cols)
         flag=1
# for either mask.region or mask.value
         if nocols == 1:
             par = cols[0]
         elif nocols == 2:
             par = cols[1]
         elif nocols == 3:
             par = [' ', ' ']
             par[0] = cols[1] 
             par[1] = cols[2] 
        

# print par, flag, 'value and flag' 
 return [par,flag]

#c------------------------------------------------------------------------------
#      subroutine yesno (string, flag, error)
#c------------------------------------------------------------------------------
#
#c  This routine interprets a string as a yes/no answer.  If the first 
#c  non-blank character of the string is "y" or "Y" the flag is set true;
#c  if it is "n" or "N", the flag is set false; in either case, the error
#c  is set to 0.  If the string is null or starts with any other character,
#c  the error is set to 1 and the flag is unchanged.
#
#
#      integer       error
#      logical       flag
#      character*(*) string
#      character*1   char
#
#      call strlim (string, ibeg, iend)
#      if (iend .eq. 0) then
#         error = 1
#      else
#         char = string(ibeg:ibeg)
#         if ((char .eq. "y") .or. (char .eq. "Y")) then
#            flag = .true.
#            error = 0
#         else if ((char .eq. "n") .or. (char .eq. "N")) then
#            flag = .false.
#            error = 0
#         else
#            error = 1
#         end if
#      end if
#
#      return
#      end
#
#c------------------------------------------------------------------------------
#      subroutine strlim (string, ibeg, iend)
#c------------------------------------------------------------------------------
#
#c  This routine finds the first and last non-blank characters in a string,
#c  returning their positions in ibeg and iend, respectively.  If the string
#c  is null or contains only blanks, ibeg and iend are returned as 0.
#
#
#      character*(*) string
#
#      iend = lnblnk(string)
#      if (iend .gt. 0) then
#         ibeg = 1
#         do while (string(ibeg:ibeg) .eq. " ")
#            ibeg = ibeg + 1
#         end do
#      else
#         ibeg = 0
#      end if
#
#      return
#      end
#

