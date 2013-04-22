#!/usr/bin/env python

# Author                 Version      Date
# -----------------------------------------------
# S M Crawford (SAA0)    0.1          17 Apr 2010
#
# interidentify is the interactive portion of the identify
# task.  It allows the user to look at the data and identify
# spectral features.   The user can also run and check
# the matching algorythms that can be run.
#

import saltprint, saltkey, saltio, salttime

import numpy
import time

import spectools as st
from PySpectrograph import WavelengthSolution

import Tkinter as Tk
from pylab import *
from matplotlib.widgets import Cursor, SpanSelector, Slider, CheckButtons
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg



class InterIdentify:
 
 
 def __init__(self, xarr, specarr, slines, sfluxes, ws, xdiff=20, function='poly', order=3, verbose=True):
    """Plot the input spectra and allow the user to interactively or automatically identify spectral features
   
       returns status
    """
    #set up the variables
    maxcolumn=5
    self.xarr=xarr
    self.specarr=specarr  
    self.slines=slines
    self.sfluxes=sfluxes
    self.ws=ws
    self.xdiff=xdiff
    self.verbose=verbose
    self.function=function
    self.order=order
    self.niter=5
    self.sigma=5
    self.xlimit=5
    self.wlimit=10
    self.fpvalue=0.5*self.specarr.max()
 
    self.spmatch=None
    self.splines=None
    #empty lists that are the values for matched coordinates
    #self.xp=[1904.3030000000001, 1686.8838411900001, 1135.0326118600001, 2560.39483096]
    #self.wp=[4671.2259999999997, 4624.2757000000001, 4500.9772000000003, 4807.0190000000002]
    #self.fp=[10818.823630112836, 10818.823630112836, 10818.823630112836, 10818.823630112836]
    self.xp=[]
    self.wp=[]
    self.fp=[]

    #set up the GUI
    self.root= Tk.Tk()
    self.root.wm_title("InterIdentify")
    self.root.bind("<Destroy>", self.destroy)
    self.root.bind("?", self.help)
    self.root.bind("q", self.destroy)
    #self.root.bind("<Motion>", self.set_focus)
    #self.root.bind("<Button-1>", self.callback)

    #set up some variables
    self.xtext= Tk.StringVar(master=self.root )
    self.wtext= Tk.StringVar(master=self.root )
    self.xtext.set(1904.303)
    self.wtext.set(4671.226)
    self.waveaxis=Tk.IntVar(master=self.root)
    self.waveaxis.set(0)


    #set up the spectra figiure
    self.spfig=figure(figsize=(8,8),dpi=72)

    #plot the spectra
    self.plotspectra()     

    #set the rows for the data
    sprow=0
    tbrow=1
    fwrow=2
    oprow=3
    qurow=4
    
    #add spectra plot
    self.spcanvas = FigureCanvasTkAgg(self.spfig, master=self.root)
    self.spcanvas.show()

    self.spcanvas.get_tk_widget().grid(row = sprow, column = 0, columnspan = maxcolumn, sticky = 'news')
    self.spcanvas.mpl_connect('key_press_event',self.keypress)
    self.spcanvas.mpl_connect('motion_notify_event',self.set_focus)

    #add the toolbar
    self.toolbar=NavigationToolbar2TkAgg(self.spcanvas,self.root)
    self.toolbar.grid(row=tbrow, column = 0, columnspan = maxcolumn, sticky = 'news')

    #add the information panel
    self.inFrame = Tk.Frame(master=self.root)
    self.inFrame.grid(row=fwrow, column=0, columnspan=maxcolumn, sticky='nw')


    self.xLabel = Tk.Label(master=self.inFrame, fg='#000000',text='X: ', relief='solid')
    self.xLabel.grid(row=0, column=0, sticky='w')
    self.xEntry = Tk.Entry(master=self.inFrame, fg='#000000',textvariable=self.xtext, relief='solid')
    self.xEntry.grid(row=0, column=1, sticky='w')

    self.wLabel = Tk.Label(master=self.inFrame, fg='#000000',text='W: ', relief='solid')
    self.wLabel.grid(row=0, column=2, sticky='ew')
    self.wEntry = Tk.Entry(master=self.inFrame, fg='#000000',textvariable=self.wtext, relief='solid')
    self.wEntry.grid(row=0, column=3, sticky='ew')

    self.addbutton = Tk.Button(master=self.inFrame, text='Add', command=self.addmatches)
    self.addbutton.grid(row=0, column=4, sticky='e')

    self.fitbutton = Tk.Button(master=self.inFrame, text='Fit', command=self.fit)
    self.fitbutton.grid(row=0, column=5, sticky='e')

    self.updatebutton = Tk.Button(master=self.inFrame, text='Update', command=self.updateplot)
    self.updatebutton.grid(row=0, column=6, sticky='e')

    #add the options panel
    self.opFrame = Tk.Frame(master=self.root)
    self.opFrame.grid(row=oprow, column=0, columnspan=maxcolumn, sticky='nw')

    self.pixbutton= Tk.Radiobutton(master=self.opFrame, text='pixels', variable=self.waveaxis, value=0, command=self.pixorwave)
    self.pixbutton.grid(row=0, column=0, sticky='e')
    self.wavebutton= Tk.Radiobutton(master=self.opFrame, text='wavelength', variable=self.waveaxis, value=1, command=self.pixorwave)
    self.wavebutton.grid(row=0, column=1, sticky='e')

    #add the quit button
    self.quFrame = Tk.Frame(master=self.root)
    self.quFrame.grid(row=qurow, column=0, columnspan=maxcolumn, sticky='ew')
    self.exitbutton = Tk.Button(master=self.quFrame, text='Quit', command=self.exit)
    self.exitbutton.grid(row=0, column=3, sticky='ew')
	
    return 


 def runplotdata(self):

    Tk.mainloop()

    return 
    
 def destroy(self, e): 
    self.root.quit()
    return self.ws

 def exit(self):
    self.root.quit()
    return self.ws

 def set_focus(self, e):
    self.spcanvas.get_tk_widget().focus_set()
    #e.widget.focus()
    return

 def pixorwave(self):
    """Set the axis to either be in pixels or wavelength"""
    self.updateplot()
    return

 def help(self, e):
    """Print the help message and the key-bindings available to the user"""
    helpmessage="""
The following commands are available to the user:
? - Print this information           q - quit the viewer
    """
    print helpmessage
    return

 def keypress(self, e):
    """React to different key presses"""
    if e.key=='p' and e.xdata:
      #print out the value
      if self.ws:
         print e.xdata, self.ws.value(e.xdata)
      else:
         print e.xdata
    if e.key=='c' and e.xdata and e.ydata:
      #find the centroid of the line and insert into xtext
      try:
         print e.xdata, e.ydata
         cx = st.centroid(self.xarr, self.specarr, xc=e.xdata, xdiff=self.xdiff)
         print cx
      except Exception, ex:
         print ex
         cx=e.xdata
      self.xtext.set(cx)

    if e.key=='m' and e.xdata and e.ydata:
      #find the line and insert into xtext
      self.xtext.set(e.xdata)

    if e.key=='f': 
      #find the fit
      self.fit()
    if e.key=='b':
      #auto-idenitfy features
      self.autoidentify()
    if e.key=='z':
      #Assume the solution is correct and find the zeropoint
      #that best matches it from cross correllation
      self.zeropoint()
    if e.key=='e':
      #find closest feature from existing fit and line list
      #and match it
      self.matchline(e.xdata)
    if e.key=='l':
      #plot the features from existing list
      self.splines=True
      self.updateplot()
    if e.key=='i':
      #reset identified features
      pass
    if e.key=='d':
      #Delete feature 
      self.deletefeature(e.xdata)

    if e.key=='q':
      #quit
      self.root.quit()

 def zeropoint(self):
   """Assume a zero-point shift in the solution.  Find the best solution by zero-point
      shifting the results and finding the best cross-correlation and plot the result
   """
   #create the artificial spectra
   res=0.2
   dres=0.1
   lmax=self.specarr.max()
   wmin=self.ws.value(self.xarr.min())
   wmax=self.ws.value(self.xarr.max())
   print wmin, wmax, self.ws.coef
   mask=(self.slines>wmin)*(self.slines<wmax)
   sl=self.slines[mask]
   sf=self.sfluxes[mask]
   print sl, sf
   swarr, sfarr=st.makeartificial(sl, sf, lmax, res, dres)

   #Now find the best fitting coefficients for the wavelength solution
   ws=self.ws
   nws=WavelengthSolution.WavelengthSolution(ws.x_arr, ws.w_arr, order=self.order)
   nws.coef=ws.coef
 
   #create the range of coefficents
   dcoef=ws.coef*0.0
   #dcoef[-3]=ws.coef[-3]
   #dcoef[-2]=float(ws.coef[-2]*0.5)
   dcoef[-1]=10
   nstep=20
   print 'dcoef:', dcoef
   dlist=st.mod_coef(ws.coef, dcoef, 0, nstep)
   #print dlist

   #loop through them and deteremine the best cofficient
   cc_arr=np.zeros(len(dlist), dtype=float)
   for i in range(len(dlist)):
       #set the coeficient
       nws.coef=dlist[i]

       #set the wavelegnth coverage 
       warr=nws.value(self.xarr)

       #resample the artificial spectrum at the same wavelengths as the 
       asfarr=numpy.interp(warr, swarr, sfarr, left=0.0, right=0.0)

       #calculate the correlation value
       cc_arr[i]=st.ncor(self.specarr, asfarr)

   #plot the results
   nws.coef=dlist[cc_arr.argmax()]
   warr=nws.value(self.xarr)
   print ws.coef
   print nws.coef
   self.ws=WavelengthSolution.WavelengthSolution(self.xarr, warr)
   self.ws.fit()
   asfarr=numpy.interp(warr, swarr, sfarr, left=0.0, right=0.0)
   self.spcross=self.spplot.plot(self.xarr,asfarr, ls='-', color='#770000')
   self.spcanvas.draw()

   


 def autoidentify(self):
   """Auto-identify features.  If the wavelength solution is not already defined, then
      use the same program as the main one.  If the wavelength solution is defined,
      then use the initial fit to generate the better fit
   """
   if self.ws:
       #detect all of the lines in the image and determine their fluxes
       xp, xf=st.findpoints(self.xarr, self.specarr, self.sigma, self.niter)
       i_ord=xf.argsort()
       self.sppeaks=self.spplot.plot(xp, xf, ls='', marker='o', ms=5, color='#00FF00')
       self.spcanvas.draw()
 
       #find the wavelength matches
       wp=st.findmatch(self.xarr, self.specarr, xp, xf, self.slines, self.sfluxes, self.ws)
  
       #now loop through the list and if it is positive then add it to the matches
       for x, w in zip(xp, wp):
         if w>0:
           print x,w, self.ws.value(x)
           self.xp.append(x)
           self.wp.append(w)
           self.fp.append(self.fpvalue)
           self.spmatch=self.spplot.plot([x],[self.fpvalue], ls='', marker='|', ms=20, color='#770000')
   else:
       pass

   self.spcanvas.draw()


 def deletefeature(self, x):
   """Find the closest line in the matched sources and remove it"""
   #centroid the data to find the exact line
   try:
       cx = st.centroid(self.xarr, self.specarr, xc=x, xdiff=self.xdiff)
   except:
       cx=x

   #now find the closest line to it
   try:
       in_minw=abs(self.xp-cx).argmin()
       self.xp.__delitem__(in_minw)
       self.wp.__delitem__(in_minw)
       self.fp.__delitem__(in_minw)
   except Exception,e:
       message='WARNING--Cannot delete feature due to %s ' % e
       print message
       return

   self.updateplot()

 def matchline(self, x):
   """Find the closest line in the database and match it"""
   #centroid the data to find the exact line
   try:
       cx = st.centroid(self.xarr, self.specarr, xc=x, xdiff=self.xdiff)
   except:
       cx=x

   #find the wavelength of it
   if self.ws:
       w=self.ws.value(cx)
   else:
       message='WARNING--No solution has been found yet'
       print message
       return

   #now find the closest line to it
   try:
       in_minw=abs(self.slines-w).argmin()
       self.xp.append(cx)
       self.wp.append(self.slines[in_minw])
       self.fp.append(0.5*self.specarr.max())
   except Exception,e:
       message='WARNING--Cannot find a matching feature due to %s ' % e
       print message
       return

   #now plot the value
   self.spmatch=self.spplot.plot([cx], [self.fpvalue], ls='', marker='|', ms=20, color='#FF00FF')
   self.spcanvas.draw()

 def fit(self):
   """Fit the data and calculate the Wavelength Solution"""
   if self.xp and self.wp:
      self.ws=WavelengthSolution.WavelengthSolution(self.xp, self.wp, order=self.order)
      self.ws.fit()
      for x,w in zip(self.xp, self.wp):
       print x,w, self.ws.value(x)
      print self.ws.coef
	
 def addmatches(self):
    """Add a value to the list for matched values"""

    #First get the two values and convert to float
    try:
       x=float(self.xtext.get())
    except ValueError:
       print 'ERROR--Could not determine x-value'
       return
    try:
       w=float(self.wtext.get())
    except ValueError:
       print 'ERROR--Could not determine w-value'
       return

    #add the two values to the lists
    self.xp.append(x)
    self.wp.append(w)
    self.fp.append(self.fpvalue)

    #add the point to the plot
    if not self.spmatch: self.spmatch=True
    self.updateplot()

 def plotlines(self):
   """Plot the spectra lines using the existing values for the fit"""
   if not self.ws: 
       return

   #plot the values
   xl=self.ws.invvalue(self.slines)
   print self.ws.coef
   print self.slines, xl
   fl=xl*0.0+0.25*self.specarr.max()

   #update the plot and add the lines
   self.splines=self.spplot.plot(xl, fl, ls='', marker='|', ms=20, color='#00FF00')
 

 
    
 def updateplot(self):
   """Handle updating the light curve plot and the data array plot when the 
       data array image is changed
   """

   for i in range(len(self.xp)):
     pass#print self.xp[i], self.wp[i], self.fp[i]

   #update the canvas
   self.spfig.delaxes(self.spplot)
   self.plotspectra()
   if self.spmatch: 
       self.spmatch=self.spplot.plot(self.xp, self.fp, ls='', marker='|', ms=20, color='#FF0000')
   if self.splines: 
       self.plotlines()

   #update the canvas
   self.spcanvas.draw()


 def plotspectra(self): 
    """Plot the spectra

    """
    self.spplot = self.spfig.add_axes([0.15,0.10,0.8,0.80], autoscale_on=False, adjustable='datalim'  )
    self.spplot.hold(True)

    #plot the curve
    self.spcurve,=self.spplot.plot(self.xarr,self.specarr,linewidth=0.5,linestyle='-',marker='',color='b')


    #set the axis limits
    self.spplot.set_xlim(self.xarr.min(),self.xarr.max())
    self.spplot.set_ylim(self.specarr.min(),self.specarr.max())

    #set the labels
    self.spplot.set_ylabel('Counts')
    self.spplot.set_xlabel('Pixels')

    #set up the other wavelength solution
    if self.ws and self.waveaxis.get():
       self.spplot.set_xlabel('Wavelength')
       xt=self.spplot.get_xticks()
       xt=self.ws.value(xt)
       xtl=[]
       for k in xt: xtl.append('%i' % k)
       self.spplot.set_xticklabels(xtl)


