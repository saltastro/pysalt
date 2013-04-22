
import numpy as np
import numpy.lib.recfunctions as rfn
from xml.dom.minidom import Document

class SlitError(Exception):
    """Base class for exceptions from the slitlet class"""
    pass


class Slitlets:
    """Slitlets is a class describing slitlets.   It has properties
       which describe the slitlet including information that would be
       included by the user and the shape of the slitlet.

       name: string character description of the target
       ra : right acesion of the target in decimal degrees
       dec: declination of the target in decimal degrees
       equinox:  equinox of the ra/dec 
       priority:  priority of the target
       mag: magnitude of the target
       band: pass band the magnitude was measured in
       selecte: d_flag:  whether the slitlet is currently selected for a mask
       refstar_flag:  whether the slitlet is a reference star
    """
    def __init__(self, shape=(10), default_width=1.5, default_length=10.0):
        self.dnames=('name', 'targ_ra', 'targ_dec', 'equinox', 'mag', 'band', 'priority', 'width', 'len1', 'len2', 'tilt',
                     'slit_ra', 'slit_dec', 'slit_width', 'slit_len1', 'slit_len2', 'slit_tilt', 'slit_radius', 
                     'inmask_flag', 'collision_flag', 'refstar_flag', 'fov_flag', 'collision_id', 'flags'
                     )
        self.dformat=('S30', 'f4', 'f4', 'i4', 'f4', 'S1', 'f4', 'f4', 'f4',  'f4', 'f4',
                      'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
                      'i4', 'i4', 'i4', 'i4', 'S30', 'S128'
                      )
        self.dtype=zip(self.dnames, self.dformat)
        self.data=None
        self.default_width=default_width
        self.default_length=default_length

    def add_arrays(self, x,y):
        return rfn.stack_arrays((x,y), usemask=False)

    def findslitlet(self, value):
        """Find the slitlet object with a name given by value"""
        name=self.data['name']
        return np.where(name == value)[0][0]

    def findtarget(self, ra, dec):
        """Find the slitlet object with the closest ra nad dec"""
        dist = ((self.data['targ_ra'] - ra)**2 + (self.data['targ_dec'] - dec)**2)**0.5
        return dist.argmin()

    def addtomask(self, sid):
        """Add object with a value of sid to the mask"""
        self.data['inmask_flag'][sid] = 1

    def updatevalue(self, sid, value, column):
        """Given a slitlet, update a column with the appropriate value"""
        dformat=self.dformat[self.dnames.index(column)]
        value=self.parseValue(value, dformat)
        self.data[column][sid] = value

    def parseValue(self, value, dformat):
        if dformat.count('S'):
           return str(value)
        elif dformat.count('i'):
           return int(value)
        elif dformat.count('f'):
           return float(value)
        else:
            message='%s is not a supported  format' % dformat
            raise SlitError(message)

    def make_default(self, dformat):
        '''
        depending on the format of the data type, either return a empty string
        or a 0
        '''
        if dformat.count('S'):
           return ''
        return 0
    
    def create_default_slitlet(self):
        '''
        create a default slit using the class dtype
        '''
        
        slitlet = np.rec.array(tuple([self.make_default(i) for i in self.dformat]),dtype=self.dtype)
        return slitlet

    def add_slitlet(self,isrefstar=0):
        '''
        add a default slitlet to the rec array and set the inmask flag = 1 so
        that the new slitlet shows in the slit table or in the reference star
        table
        '''
        
        slitlet = self.create_default_slitlet()
        if isrefstar:
            slitlet['refstar_flag'] = 1
        slitlet['inmask_flag'] = 1
        self.data = self.add_arrays(self.data, slitlet)
    

    def update_slit_flags(self,i):
        '''
        write a string containing all the flags to be displayed in the flags
        column of the tables. this method only check the flags for one row
        in the rec array
        '''
        txt = ''
        if self.data[i]['inmask_flag']:
            txt = txt + 'inMask '

        if self.data[i]['fov_flag']:
            txt = txt + 'inFoV '
        else:
            txt = txt + 'outFoV '

        if self.data[i]['refstar_flag']:
            txt = txt + 'Ref '

        if self.data[i]['collision_flag']:
            txt = txt + 'Col ' + self.data[i]['collision_id']

        self.data[i]['flags'] = txt.strip()

    def update_flags(self):
        '''
        update the flags for all the objects
        '''

        for i in range(0,len(self.data)):
            self.update_slit_flags(i)

    def readascii(self, infile, form='short'):

        if form=='short':
            dnames=('name', 'targ_ra', 'targ_dec', 'equinox', 'mag', 'band', 'priority')
            dformat=('S30', 'f4', 'f4', 'i4', 'f4', 'S1', 'f4')
            object_arr=np.loadtxt(infile, dtype={'names': dnames, 'formats': dformat},
                         converters={1:ra_read, 2:dec_read})
            #determine the missing values
            mnames=[]
            mtypes=[]
            for i in range(len(self.dnames)):
                if self.dnames[i] not in dnames:
                   mnames.append(self.dnames[i])
                   mtypes.append(self.dformat[i])
            #set up the default values
            default_list=[np.zeros(len(object_arr))]*len(mnames)
            default_list[0]=default_list[0]+self.default_width
            default_list[1]=default_list[1]+0.5*self.default_length
            default_list[2]=default_list[2]+0.5*self.default_length
            object_arr=rfn.append_fields(object_arr, names=mnames, data=default_list, dtypes=mtypes,
                     fill_value=0, usemask=False)
        elif form=='long':
            dnames=('name', 'targ_ra', 'targ_dec', 'equinox', 'mag', 'band', 'priority', 'width', 'length', 'tilt')
            dformat=('S30', 'f4', 'f4', 'i4', 'f4', 'S1', 'f4', 'f4', 'f4', 'f4')
            object_arr=np.loadtxt(infile, dtype={'names': dnames, 'formats': dformat},
                         converters={1:ra_read, 2:dec_read})
            #determine the missing values
            mnames=[]
            mtypes=[]
            for i in range(len(self.dnames)):
                if self.dnames[i] not in dnames:
                   mnames.append(self.dnames[i])
                   mtypes.append(self.dformat[i])
            #set up the default values
            default_list=[np.zeros(len(object_arr))]*len(mnames)
            object_arr=rfn.append_fields(object_arr, names=mnames, data=default_list, dtypes=mtypes,
                     fill_value=0, usemask=False)
        else:
            message='This format is not supported'
            raise SlitError(message)

        #set objects that are preselected
        object_arr['inmask_flag'] = 1.0*(object_arr['priority'] >= 1.0)
        #set reference stars
        object_arr['refstar_flag'] = 1.0*(object_arr['priority'] == -1.0)
     
        #stack the data if it already exists
        if self.data is None:  
           self.data=object_arr
        else:
           self.data=self.add_arrays(self.data, object_arr)
        # total number of objects:
        self.nobjects=len(self.data)
        self.update_flags()

    def readxml(self,slits_dict,refstars_dict):
        '''
        read the dictionaries produced by slitmask.readmaskxml and populate the
        slitlets rec array.
        '''
        slitdnames=('name', 'targ_ra', 'targ_dec', 'mag', 'priority', 'width', 'len1','len2')
        slitdformat=('S30', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4','f4')
        slitdtype=zip(slitdnames, slitdformat)

        refstardnames=('name', 'targ_ra', 'targ_dec', 'mag')
        refstardformat=('S30', 'f4', 'f4', 'f4')
        refstardtype=zip(refstardnames, refstardformat)


        # map the dict values to the recarray names
        map = {'name':'id','targ_ra':'xce', 'targ_dec':'yce','mag':'mag',\
        'priority':'priority','width':'width','len1':'len1','len2':'len2'}

        slitlist = []
        for i in slits_dict.keys():
            tmp = [slits_dict[i][map[j]] for j in slitdnames]
            slitlist.append(tmp)

        slits_arr = np.rec.array(slitlist,slitdtype)
        #determine the missing values
        mnames=[]
        mtypes=[]
        for i in range(len(self.dnames)):
            if self.dnames[i] not in slitdnames:
               mnames.append(self.dnames[i])
               mtypes.append(self.dformat[i])
        #set up the default values
        default_list=[np.zeros(len(slits_arr))]*len(mnames)
        slits_arr=rfn.append_fields(slits_arr, names=mnames, data=default_list,\
            dtypes=mtypes, fill_value=0, usemask=False)

        refstarlist = []
        for i in refstars_dict.keys():
            tmp = [refstars_dict[i][map[j]] for j in refstardnames]
            refstarlist.append(tmp)

        refstars_arr = np.rec.array(refstarlist,refstardtype)
        mnames=[]
        mtypes=[]
        for i in range(len(self.dnames)):
            if self.dnames[i] not in refstardnames:
               mnames.append(self.dnames[i])
               mtypes.append(self.dformat[i])

        #set up the default values
        default_list=[np.zeros(len(refstars_arr))]*len(mnames)
        refstars_arr=rfn.append_fields(refstars_arr, names=mnames, data=default_list, dtypes=mtypes,
                 fill_value=0, usemask=False)
        refstars_arr['priority'] = -1
        refstars_arr['refstar_flag'] = 1

        
        object_arr = self.add_arrays(slits_arr, refstars_arr)
        #set objects that are preselected to be in the mask
        object_arr['inmask_flag'] = 1.0*(object_arr['priority'] >= 1.0)
        
        #stack the data if it already exists
        if self.data is None:
           self.data=object_arr
        else:
           self.data=self.add_arrays(self.data, object_arr)
        self.nobjects=len(self.data)
        self.update_flags()

        
    def asregion(self, i, pa=0.):
        """Create a ds9 region string from the shape of the object 
        *TODO*
        """
        ra=self.data['targ_ra'][i] 
        dec=self.data['targ_dec'][i]
        if self.data['refstar_flag'][i]:
            shape='circle'
            radius=3#self.data['radius'][i]
#            regstr='%s(%f,%f,%f")' % (shape, ra, dec, radius)
            regstr='%s(%f,%f,%f") # color={blue} tag = {ref} ' % (shape, ra, dec, radius)
        else:
            shape = 'box'
            name = self.data['name'][i]
            width = self.data['width'][i]
            length = self.data['len1'][i]+self.data['len2'][i]
            # needs to be posang+tilt ?
            tilt = self.data['tilt'][i]-pa#-90.0
#            regstr='%s(%f,%f,%f",%f",%f) # tag = {slit} tag = {%s} ' % (shape, ra, dec, width, length, tilt, name)
            color = 'green'
            if self.data['collision_flag'][i] == 1:
                color = 'red'
            print name, self.data['collision_flag'][i], self.data['collision_id'][i]
            regstr = '%s(%f,%f,%f",%f",%f) # color={%s} tag = {slit} tag = {%s}\n ' % (shape, ra, dec, width, length, tilt, color, name)
            #regstr+='%s(%f,%f,%f",%f",%f) # color={green} tag = {slit} tag = {%s} ' % (shape, ra, dec, 1000, length, tilt, name)
        return regstr

    def asregionspec(self, i, pa=0.):
        """Create a ds9 region string from the shape of the object 
        *TODO*
        """
        ra = self.data['targ_ra'][i]
        dec = self.data['targ_dec'][i]
        shape = 'box'
        name = self.data['name'][i]
        width = self.data['width'][i]
        length = self.data['len1'][i]+self.data['len2'][i]
        # needs to be posang+tilt ?
        tilt = self.data['tilt'][i]+pa#-90.0
#            regstr='%s(%f,%f,%f",%f",%f) # tag = {slit} tag = {%s} ' % (shape, ra, dec, width, length, tilt, name)
        #regstr='%s(%f,%f,%f",%f",%f) # color={red} tag = {slit} tag = {%s}\n ' % (shape, ra, dec, width, length, tilt, name)
        regstr = '%s(%f,%f,%f",%f",%f) # color={green} tag = {spec} tag = {%s} ' % (shape, ra, dec, 1000, length, tilt, name)
        return regstr


    def asxml(self, i):
        doc=Document()
        if self.data['refstar_flag'][i]:
            card=doc.createElement("refstar")
            card.setAttribute("id", "%s"%(str(self.data['name'][i])))
            card.setAttribute("xce", "%f"%(self.data['targ_ra'][i]))
            card.setAttribute("yce", "%f"%(self.data['targ_dec'][i]))
            card.setAttribute("width", "%f"%(self.data['width'][i]))
            card.setAttribute("length", "%f"%(self.data['len1'][i]+self.data['len2'][i]))
            card.setAttribute("priority", "%f"%(self.data['priority'][i]))
            card.setAttribute("mag", "%f"%(self.data['mag'][i]))
        else:
            card=doc.createElement("slit")
            card.setAttribute("id", "%s"%(str(self.data['name'][i])))
            card.setAttribute("xce", "%f"%(self.data['targ_ra'][i]))
            card.setAttribute("yce", "%f"%(self.data['targ_dec'][i]))
            card.setAttribute("width", "%f"%(self.data['width'][i]))
            card.setAttribute("length", "%f"%(self.data['len1'][i]+self.data['len2'][i]))
            card.setAttribute("priority", "%f"%(self.data['priority'][i]))
            card.setAttribute("mag", "%f"%(self.data['mag'][i]))
        return card
 

class slitshape:
    def __init__(self, ra, dec, width=1, length=5, angle=0):
        """Describe the shape for a typical slit.  It is a rectangle
           with a set height and width.  The RA and DEC are the center
           position of the slit and maybe different from the target

           ra--Central Right Ascension in decimal degrees of slit
           dec--Central Declination in decimal degress of slit
           width--width in arcseconds
           height--height in arcseconds
           angle--position angle in degrees with respect to the position angle
                  of the slitmask
        """

        self.ra=ra
        self.dec=dec
        self.width=width
        self.length=length
        self.angle=angle
        self.name='rectangle'

    def findslitlet(value):
        """For a value given for the slitlet name, return that slitlet"""
        try:
            x=np.where(self.data['name']=='8')[0][0]
            return self.data[x]
        except Exception,e:
            raise SlitError(e)


class starshape:
    def __init__(self, ra, dec, radius=1):
        """Describe the shape for the reference star.  It is a circle
           with a set size and radius

           ra--Right Ascension in decimal degrees
           dec--Declination in decimal degress
           radius--radius in arcseconds
        """

        self.ra=ra
        self.dec=dec
        self.radius=radius
        self.name='circle'


def sex2dec(x):
    x=x.split(':')
    if float(x[0])>=0:
        return float(x[0])+float(x[1])/60.0+float(x[2])/3600.0
    else:
        return -(abs(float(x[0]))+float(x[1])/60.0+float(x[2])/3600.0)

def ra_read(x):
    try:
       return float(x)
    except ValueError,e:
       try:
           return 15*sex2dec(x)
       except: 
           return None

def dec_read(x):
    try:
       return float(x)
    except ValueError,e:
       try:
           return sex2dec(x)
       except ValueError,e:
           return None
 


if __name__=='__main__':
   import sys
   s=Slitlets()
   s.readascii(sys.argv[1], form='short')
   print s.data[s.data['inmask_flag']==1]
   print s.data.dtype
   #object_arry, slit_dict=readasciicatalog(sys.argv[1], form='short')
   #print object_arry[0]
   #s=slit_dict['177']
   #print s.shape.name, s.selected  
   #print s.asregion()
