#! /usr/bin/env python
"""
saltarith is a copy of IRAF's imarith routine
but designed to deal with SALT's multi extension
data in a smart way. 

Author                 Version      Date
-----------------------------------------------
A Eigenbrot (Madison)      1.0       16  July 2012
                           1.1       05  Aug  2012
S.M. Crawford (SAAO)       1.2       13 Nov   2012

TODO
-----------------------------------------------


Updates
-----------------------------------------------
15 Aug 2012  Removed support for using the mean, as it can
             only lead to Bad Things. Fits header is now
             updated just like all other pysalt routines.

13 Nov   2012 -Added paremter file
              -Added to PySALT package


"""
from __future__ import with_statement

import time, numpy
from pyraf import iraf

import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging, history

from salterror import SaltError

debug=True

def saltarith(operand1, op, operand2, result, outpref, divzero=0, clobber=False,     \
             logfile='salt.log',verbose=True):


   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',operand1)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', result, outpref, infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')

       #let's keep track of whether operand2 is an image or not
       is_image = False

       #load in operand2, or, if it's not an image, assume it's a number
       try: operand2struct = float(operand2)
       except ValueError:
          operand2struct = saltio.openfits(operand2)
          is_image = True

       #open the input image files
       for infile, outfile in zip(infiles,outfiles):
           struct = saltio.openfits(infile)

           #do some math!
           outstruct = arith(struct, op, operand2struct, is_image, divzero)
           try:
               pass 
           except Exception,e:
               msg='Unable to do math %s because %s' % (infile, e)
               raise SaltError(msg)

           #update header stuff
           fname, hist = history(level=1, wrap=False)
           saltkey.housekeeping(struct[0],'SARITH', 'Some arithmatic was performed',hist)

           #write it. close it.
           saltio.writefits(outstruct,outfile,clobber=clobber)
           saltio.closefits(struct)

           #output the information
           log.message('imarith: %s %s %s %s' % (infile, op, operand2, outfile), with_header=False, with_stdout=verbose)

       #close the operand2 image
       if is_image: saltio.closefits(operand2struct)

# -----------------------------------------------------------
# Actually do the math

def arith(struct, op, opstruct, is_image, divzero):
   """performs op on struct with opstruct as the argument

      return struct
   """
   # Determine the number of extensions
   nextend=len(struct)

   #do the math
   for i in range(nextend):
      if struct[i].name!='PRIMARY' or len(struct)==1:
         #screw variance frames. Pay me more
         #actually do the math
         if is_image: 
            struct[i].data = eval('struct[i].data'+str(op)+'opstruct[i].data')
            
            if op == '/':
               #find where division by zero is going to happen
               zidx = numpy.where(opstruct[i].data == 0)
               struct[i].data[zidx] = divzero

         else: 
            struct[i].data = eval('struct[i].data'+str(op)+'opstruct')
            
            if op == '/' and opstruct == 0:
               struct[i].data = numpy.ones(struct[i].data.shape)*divzero

   return struct
            
# -----------------------------------------------------------
# main code

if not iraf.deftask('saltarith'):
   parfile = iraf.osfn("saltred$saltarith.par")
   t = iraf.IrafTaskFactory(taskname="saltarith",value=parfile,function=saltarith, pkgname='saltred')
