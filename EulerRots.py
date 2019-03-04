import pandas as pd
import numpy as np
import copy,os
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cmx
from mpl_toolkits.basemap import Basemap, shiftgrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
import RotKit_f as rotkit_f
import cartopy.crs as ccrs

#Imports a table for looking up plate codes - potentially useful when building new rotations?
platecodes=pd.read_table(os.path.join(__location__,'Datafiles/PlateCodes.txt'), header=None, names=['NumCode','LetterCode','Description'],index_col='NumCode')

#adds a new ellipse function to Basemap class
def ellipse(self, lon, lat, a, b, az, **kwargs):
    degrad=np.pi/180.
    ax = self._check_ax()
    ellipse_angles=np.arange(0,360,0.5)*degrad
    a,b=a*np.pi/180.,b*np.pi/180.
    ellipser=(a*b)/np.sqrt((a**2*np.sin(ellipse_angles)**2)+(b**2*np.cos(ellipse_angles)**2))
    ellipse_coords=sphere_point_along_bearing(lat,lon,ellipse_angles/degrad+az,ellipser/degrad)
    x, y = self(ellipse_coords[1],ellipse_coords[0])
    poly = Polygon(zip(x,y), **kwargs)
    ax.add_patch(poly)
    # Set axes limits to fit map region.
    self.set_axes_limits(ax=ax)
    return poly

Basemap.ellipse = ellipse

def ellipse_pars(lon, lat, a, b, az, **kwargs):
    """
    given a specified error ellipse will generate a polygon patch 
    that can then be added to a map: currently set up for use with Cartopy.
    """
    degrad=np.pi/180.
    ellipse_angles=np.arange(0,360,0.5)*degrad
    a,b=a*np.pi/180.,b*np.pi/180.
    ellipser=(a*b)/np.sqrt((a**2*np.sin(ellipse_angles)**2)+(b**2*np.cos(ellipse_angles)**2))
    ellipse_coords=sphere_point_along_bearing(lat,lon,ellipse_angles/degrad+az,ellipser/degrad)
    poly = Polygon(zip(ellipse_coords[1],ellipse_coords[0]), transform=ccrs.PlateCarree(),**kwargs)
    return poly

def get_chron_age(chron,timescale='CK95'):
	"""
	Given a chron name as a string with appended 'y', 'm' or 'o', will attempt to look up 
	and return the age of the end/young end, middle or beginning/old end of the designated
	chron. Currently two timescales are available: 'CK95' (Cande & Kent 1995, default) and
	'GTS12' (Geological Timescale with some astronomically tuned chron boundaries)
	"""
    agemodel=pd.read_table(os.path.join(__location__,'Datafiles/'+timescale+'.txt'))
    if chron[-2]=='n':
        select=agemodel[agemodel.Chron==chron[1:-1]].index[0]
        if chron[-1]=='y': age=agemodel.iloc[select].Young_Age
        elif chron[-1]=='o': age=agemodel.iloc[select].Old_Age
        elif chron[-1]=='m': age=agemodel.iloc[select].Young_Age+0.5*(agemodel.iloc[select].Old_Age-agemodel.iloc[select].Young_Age)
    else: #if reversed (unlisted) find the equivalent normal chron
        select=agemodel[agemodel.Chron==chron[1:-2]+'n'].index[0]
        if chron[-1]=='y': age=agemodel.iloc[select].Old_Age
        elif chron[-1]=='o': age=agemodel.iloc[select+1].Young_Age
        elif chron[-1]=='m': age=agemodel.iloc[select].Old_Age+0.5*(agemodel.iloc[select+1].Young_Age-agemodel.iloc[select].Old_Age)  
    return age    

def find_plate_from_name(text):
    """
    Returns (in the form of a pandas DataFrame) the numerical and text codes and descriptions
    of plates with the input string fragment 'text' within their descriptions. 
    Warning: can return multiple matches! (e.g., Pacific will return any matches with 'Pacific' 
    in the plate description, which is not just the Pacific Plate). 
    """
    return platecodes[platecodes.Description.str.contains(text)]

def find_plate_from_number(code):
	"""
	Given a numerical plate code, will return a pandas DataFrame row with the 
	matching plate code, text code, and description 
	"""
    return platecodes.loc[code]

def rotmat_to_pole(rot_matrix):
	"""
	Converts given a 3x3 rotation matrix rot_matrix into a Euler rotation (lat,long, angle)
	"""
    pole_lon=np.arctan((rot_matrix[0,2]-rot_matrix[2,0])/(rot_matrix[2,1]-rot_matrix[1,2]))*180/np.pi
    if rot_matrix[2,1]-rot_matrix[1,2]<0. : pole_lon=pole_lon+180
    if pole_lon>180: pole_lon=pole_lon-360
    if pole_lon<-180: pole_lon=pole_lon+360
    toss=np.sqrt((rot_matrix[2,1]-rot_matrix[1,2])**2+(rot_matrix[0,2]-rot_matrix[2,0])**2+(rot_matrix[1,0]-rot_matrix[0,1])**2)
    pole_lat=np.arcsin((rot_matrix[1,0]-rot_matrix[0,1])/toss)*180/np.pi
    temp=(rot_matrix[0,0]+rot_matrix[1,1]+rot_matrix[2,2]-1.0)
    pole_ang=np.arctan(toss/(rot_matrix[0,0]+rot_matrix[1,1]+rot_matrix[2,2]-1.0))*180/np.pi
    if temp<0: pole_ang=pole_ang+180
    return [pole_lat,pole_lon,pole_ang]

def sphere_ang_dist(lat1,long1,lat2,long2,degrees=True):
    """
    calculates the length of the great circle path between points (lat1, long1) 
    and (lat2, long2) using haversine formula. By default, assumes input coordinates 
    are in degrees and returns distance in degrees. If degrees=False, then assumes input 
    is in radians and returns distance in radians
    
    Source: http://www.movable-type.co.uk/scripts/latlong.html    
    """
    if degrees==True: degrad=np.pi/180.
    else: degrad=1.
    dlong=(long2-long1)*degrad
    dlat=(lat2-lat1)*degrad
    lat1=lat1*degrad
    lat2=lat2*degrad 
    a=np.sin(dlat/2)*np.sin(dlat/2)+np.cos(lat1)*np.cos(lat2)*np.sin(dlong/2)*np.sin(dlong/2)
    sepdist=2*(np.arctan2(np.sqrt(a),np.sqrt(1-a)))
    return sepdist/degrad
   
def sphere_bearing(lat1,long1,lat2,long2,degrees=True):
    """
    returns bearing of point (lat2, long2) wrt point (lat1,long1). By default, assumes 
    input coordinates are in degrees and returns distance in degrees. If degrees=False, 
    then assumes input is in radians and returns distance in radians.
    
    Source: http://www.movable-type.co.uk/scripts/latlong.html
    """
    if degrees==True: degrad=np.pi/180.
    else: degrad=1.
    dlong=(long2-long1)*degrad
    lat1=lat1*degrad
    lat2=lat2*degrad
    y = np.sin(dlong)*np.cos(lat2)
    x = np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(dlong)
    if degrees==True:
    	bearing=(np.arctan2(y,x))/degrad
    	if (bearing<0): bearing=bearing+360
    return bearing
    
def sphere_dist_bearing(lat1,long1,lat2,long2,degrees=True)
	"""
	returns a list containing angular distance and bearing of point (lat2, long2) wrt point
	(lat1,long1). By default, assumes input coordinates are in degrees and returns distance 
	in degrees. If degrees=False, then assumes input is in radians and returns distance in radians.
	"""
	return [sphere_ang_dist(lat1,long1,lat2,long2,degrees), sphere_bearing(lat1,long1,lat2,long2,degrees)]    
    

def sphere_point_along_bearing(lat1,long1,bearing,d,degrees=True):
   """
   returns lat and long of a point angular distance d at azimuth bearing from initial point
   lat1,long1.  By default, assumes input coordinates are in degrees and returns distance 
   in degrees. If degrees=False, then assumes input is in radians and returns distance in radians.
   """        
   if degrees==True: degrad=np.pi/180.
   else: degrad=1.
   lat1,long1,bearing,d=lat1*degrad,long1*degrad,bearing*degrad,d*degrad
   lat2=np.arcsin(np.sin(lat1)*np.cos(d)+np.cos(lat1)*np.sin(d)*np.cos(bearing))
   long2=long1+np.arctan2(np.sin(bearing)*np.sin(d)*np.cos(lat1),np.cos(d)-np.sin(lat1)*np.sin(lat2))
   return lat2/degrad,long2/degrad
                                                                                        
def dir2cart(d):
    """
    converts list or array of vector directions [D,I] to array of cartesian coordinates, in x,y,z
    Source: PMagPy (Lisa Tauxe) https://github.com/PmagPy/PmagPy
    """
    ints=np.ones(len(d)).transpose() # get an array of ones to plug into dec,inc pairs
    d=np.array(d)
    rad=np.pi/180.
    if len(d.shape)>1: # array of vectors
        decs,incs=d[:,0]*rad,d[:,1]*rad
        if d.shape[1]==3: ints=d[:,2] # take the given lengths
    else: # single vector
        decs,incs=np.array(d[0])*rad,np.array(d[1])*rad
        if len(d)==3: 
            ints=np.array(d[2])
        else:
            ints=np.array([1.])
    cart=np.array([ints*np.cos(decs)*np.cos(incs),ints*np.sin(decs)*np.cos(incs),ints*np.sin(incs)]).transpose()
    return cart

def cart2dir(cart):
    """
    converts list of cartesian [x,y,z] values to spherical coordinates [D,I]. 
    Source: PMagPy (Lisa Tauxe) https://github.com/PmagPy/PmagPy
    """
    cart=np.array(cart)
    rad=np.pi/180. # constant to convert degrees to radians
    if len(cart.shape)>1:
        Xs,Ys,Zs=cart[:,0],cart[:,1],cart[:,2]
    else: #single vector
        Xs,Ys,Zs=cart[0],cart[1],cart[2]
    Rs=np.sqrt(Xs**2+Ys**2+Zs**2) # calculate resultant vector length
    Decs=(np.arctan2(Ys,Xs)/rad)%360. # calculate declination taking care of correct quadrants (arctan2) and making modulo 360.
    try:
        Incs=np.arcsin(Zs/Rs)/rad # calculate inclination (converting to degrees) # 
    except:
        print 'trouble in cart2dir' # most likely division by zero somewhere
        return np.zeros(3)
        
    return np.array([Decs,Incs,Rs]).transpose() # return the directions list

def angle(D1,D2):
    """
    call to angle(D1,D2) returns array of angles between lists of two directions D1,D2 where D1 is for example, [[Dec1,Inc1],[Dec2,Inc2],etc.]
    Source: PMagPy (Lisa Tauxe) https://github.com/PmagPy/PmagPy
    """
    X1=dir2cart(np.array(D1)) # convert to cartesian from polar
    X2=dir2cart(np.array(D2))
    angles=[] # set up a list for angles
    for k in range(X1.shape[0]): # single vector
        angle= np.arccos(np.dot(X1[k],X2[k]))*180./np.pi # take the dot product
        angle=angle%360.
        angles.append(angle)
    return np.array(angles)

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
class EulerRotationModel(object):
    """
    A class which organises the sets of rotations associated with particular plates, initially loaded from an input rotation file.
    Required header/columns for file: 'MovingPlate','FixedPlate','Chron','EndAge','RotLat','RotLong','RotAng','Kappahat','a','b','c','d','e','f','Points','Segs','Plates','DOF','Source'
    where a,b,c,d,e,f are parameters from the covariance matrix; Points, Segs, Plates and DOF/Degrees of Freedom are parameters for the original Hellinger fit for the specified rotation.
    
    Attributes:
    
    - rotationsets: list of FiniteRotationSets for a particular plate pair 
    - FixedPlate (currently set to 'None' and not currently used for anything)
    
    """      
    def __init__(self, rotfile):
        self.FixedPlate='None'
        rot_data=pd.read_table(rotfile)
        rot_data['StartAge']=0
        rot_data['StartChron']='None'
        self.rotationsets=[]
        for i,row in rot_data.iterrows():
            if i==0:
                CurrentMoving=row.MovingPlate
                CurrentFixed=row.FixedPlate
                rotations=[]
            if row.MovingPlate==CurrentMoving and row.FixedPlate==CurrentFixed:
                rotations.append(EulerRotation(row)) 
            else:
                self.rotationsets.append(FiniteRotationSet(rotations,CurrentMoving,CurrentFixed))
                CurrentMoving=row.MovingPlate
                CurrentFixed=row.FixedPlate
                rotations=[EulerRotation(row)]
        self.rotationsets.append(FiniteRotationSet(rotations,CurrentMoving,CurrentFixed))
        
    def insert_rots(self,newrots):
        """
        given a list of FiniteRotationSets (either calculated or from another source), 
        adds them to the rotation model. 
        NB: currently no consistency checking, so can add conflicting rotations, which would be bad. 
        """
        #todo: check for conflicting rotations before adding in new ones.
        for rotationset in newrots:
            self.rotationsets.append(rotationset)
    
    def remove_rots(self,plate1,plate2=-1):
        """
        Remove FiniteRotationSets for the pair plate1, plate2. 
        If plate2 is not defined, it will removed any FiniteRotationSet that involves plate1.
        """
        newrotset=[]
        for rotationset in self.rotationsets:
            #if no plate 2 specified, removes any rotationset with fixed or moving plate with plate1 code.
            #or, if you want to be pedantic, doesn't add it to an updated rotationset list.
            if plate2==-1:
                if (rotationset.MovingPlate==plate1 or rotationset.FixedPlate==plate1) is False:
                    newrotset.append(rotationset)
            #if plate 2 is specified, looks for the plate1-plate2 (or plate2-plate1) pair.
            else:
                if (rotationset.MovingPlate==plate1 and rotationset.FixedPlate==plate2) is False and (rotationset.MovingPlate==plate2 and rotationset.FixedPlate==plate1) is False:
                    newrotset.append(rotationset)
        self.rotationsets=newrotset

    def find_pairs(self,plate):
		"""
		Searches for any FiniteRotationSets in which the specified plate is one of the plate pair;
		returns the pairs of platecodes for those FiniteRotationSets (with plate always listed first). 
		"""    
        circuits=[]
        pairs=[rotationset for rotationset in self.rotationsets if (rotationset.MovingPlate==plate)]
        for pair in pairs:
            circuits.append([pair.MovingPlate,pair.FixedPlate])
        pairs=[rotationset for rotationset in self.rotationsets if (rotationset.FixedPlate==plate)] 
        for pair in pairs:
            circuits.append([pair.FixedPlate,pair.MovingPlate])
        return circuits  
     
    def find_circuit(self,startplate,endplate,preferred=-1):
    	"""
    	Searches for ways to link startplate to endplate using rotations within the currently 
    	defined rotation model. The returned circuit is a list of platecodes that starts with 
    	startplate and ends with endplate. If multiple possible circuits exist, they will all be returned;
    	can filter by setting a preferred platecode, although this may not always restrict the list to just one! 
    	"""
        # should check to see if both startplate and endplate are included in the rotation model...
        circuits=self.find_pairs(startplate)
        plates_used=[startplate]+[pair[-1] for pair in circuits]
        while len(plates_used)<len(self.plates()): #technically, could check for a pair with endplate in it but if that exists this shouldn't be running...
            for circuit in circuits:
                for pair in self.find_pairs(circuit[-1]):
                    if pair[-1] not in plates_used: #because we don't want to go backwards!
                        plates_used.append(pair[-1])
                        circuits.append(circuit+[pair[-1]])
        found=[circuit for circuit in circuits if circuit[-1]==endplate]
        if preferred>0: found=[circuit for circuit in found if preferred in circuit[1:-1]]
        return found
            
    def get_rots(self,movingplate,fixedplate,ages=[],preferred=-1):
    	"""
    	Returns a FiniteRotationSet for the movingplate-fixedplate pair. If a list of ages is supplied,
    	then the rotationset will consist of the interpolated finite rotations for those ages; otherwise it will be all
    	available ages in the rotationset. If more than one circuit is possible, then the circuit can be forced through preferred.
    	
    	If multiple movingplate-fixedplate rotationsets already exist within the rotation model (hopefully unlikely), 
    	or no viable plate circuit can be found, an empty list is returned. 
    	If multiple viable plate circuits are found, it will use the first one in the list returned from find_circuit().
    	"""
        #Theoretically, this should return a single rotationset, because it either exists in its stated or inverted form 
        selected1=[rotationset for rotationset in self.rotationsets if (rotationset.MovingPlate==movingplate) & (rotationset.FixedPlate==fixedplate)]
        selected2=[rotationset.invert() for rotationset in self.rotationsets if (rotationset.MovingPlate==fixedplate) & (rotationset.FixedPlate==movingplate)]
        #more than one rotationset for a defined plate pair: shouldn't happen in a properly set up rotation model/file...
        if len(selected1)+len(selected2)>1: 
            print 'Uh-oh: multiple rotations fit parameters'
            rots_got=[]
        #no rotationset for defined plate pair: next want to see if can construct a circuit from other plate pairs in the model     
        elif selected1==[] and selected2==[]:
            found=self.find_circuit(movingplate,fixedplate,preferred)
            if found: 
                toadd=[]
                #if a circuit exists, now add the circuit to get the new parameters
                #note that if more than one circuit is found, currently the first one in the list is used.
                #could potentially have an interactive prompt?
                for plate1,plate2 in zip(reversed(found[0][:-1]),reversed(found[0][1:])):
                    selected1=[rotationset for rotationset in self.rotationsets if (rotationset.MovingPlate==plate1) & (rotationset.FixedPlate==plate2)]
                    selected2=[rotationset.invert() for rotationset in self.rotationsets if (rotationset.MovingPlate==plate2) & (rotationset.FixedPlate==plate1)]
                    if selected1: toadd.append(selected1[0])
                    else: toadd.append(selected2[0])
                for rotset in toadd[1:]: rots_got=rotset.addrots(toadd[0])
            else: 
                print 'Uh-oh: no rotations fit parameters'
                rots_got=[]
        #not sure if these warning messages are the best way to go around it, but if returning nothing should cause any automated process to bork.
        elif selected1: rots_got=selected1[0]
        else: rots_got=selected2[0]
        #now interpolate ages if some have been input
        if ages:
            rots_got=rots_got.interpolate(ages)
        return rots_got
        
    def synthetic_APWP(self,moving_plate,absolute_ref_frame,ages):
    	"""
    	Returns a pandas DataFrame that predicts the Apparent Polar Wander path that should have been generated 
    	by motion of moving_plate in absolute_ref_frame (i.e reconstructed position of geographic North Pole
    	in the moving_plate reference frame) for specified list of age points. Note that currently, there is no 
    	restriction on what reference frame is used, but this will only be meaningful if it is an absolute frame 
    	(e.g, hotspot frames 001/Atlantic or 003/Pacific)
    	"""
        if ages[0]==0.: ages[0]=0.01 #Gets a bit fussy for the 0 rotation.
        reconstruction_rots=self.get_rots(moving_plate,absolute_ref_frame,ages)
        #At the reconstruction age, the VGP is at the North Pole and then drifts away from it. so use the inverted rotations
        NPole=Point(pd.Series(['VGP',moving_plate,0.,0.,90,0.],index=['Name','PlateCode','FeatureAge','ReconstructionAge','Lat','Lon']))
        VGPs=[NPole.rotate(rotation) for rotation in reconstruction_rots.invert().rotations[1:]] #zero rotation: more trouble than it's worth?   
        return pd.DataFrame([['VGP-'+`age`,moving_plate,age,0.,point.LocPars.PointLat,point.LocPars.PointLong,point.LocPars.MaxError,point.LocPars.MinError,point.LocPars.MaxBearing] for point,age in zip(VGPs,ages)],
                                columns=['Name','PlateCode','FeatureAge','ReconstructionAge','Lat','Lon','MaxError','MinError','MaxBearing']) 

    def synthetic_APWP_flowline(self,moving_plate,absolute_ref_frame,ages,SetName='APWP',PlotColor='orange',PlotLevel=5):
    	"""
    	Create an APWP object that predicts the Apparent Polar Wander path that should have been generated 
    	by motion of moving_plate in absolute_ref_frame (i.e reconstructed position of geographic North Pole
    	in the moving_plate reference frame) for specified list of age points.
    	"""
        return APWP(self.synthetic_APWP(moving_plate,absolute_ref_frame,ages),absolute_ref_frame,SetName,PlotColor,PlotLevel=5)      

    def newagemodel(self, timescale='CK95'):
        """
        WARNING: experimental. *Should* overwrite old rotationsets with new ones where rotation ages timed to reversal
        ages have been recalibrated to the specified timescale. Current options: 'CK95', 'GTS12'  
        """
        self.rotationsets=([rotationset.newagemodel(timescale) for rotationset in self.rotationsets])                       
                                                                                                                                                          
    def summary(self):
    	"""
    	Returns a summary DataFrame for FiniteRotationSets within the current rotation model: 
    	plate pair codes and names, number of finite rotations, age range.
    	"""
        return pd.DataFrame([[item.MovingPlate,find_plate_from_number(item.MovingPlate).Description,item.FixedPlate,find_plate_from_number(item.FixedPlate).Description,
                            item.N,item.rotations[0].EndAge,item.rotations[-1].EndAge] for item in self.rotationsets],
                            columns=['MovingPlate','MovingPlateName','FixedPlate','FixedPlateName','N','Youngest','Oldest'])
    def plates(self):
        """
        Returns a list of unique plate IDs within the current rotation model.
        """
        return list(set([item.MovingPlate for item in self.rotationsets]+[item.FixedPlate for item in self.rotationsets]))

class FiniteRotationSet(object):
    """
    a collection of finite EulerRotations with common moving and fixed plates. Methods to easily list set parameters such as ages, chrons,
    generate interpolated rotations and stage rotations
    """
    def __init__(self, rotations,movingplate,fixedplate):
        #where rotations is a list of EulerRotation objects
        self.MovingPlate=movingplate
        self.FixedPlate=fixedplate
        self.rotations=rotations
        #next bit makes interpolation easier: adds in a zero age rotation if none exists
        if (self.rotations[0].EndAge==0 and self.rotations[0].StartAge==0)==False: #a possibly unneccesarilty more robust check
            newrot=copy.deepcopy(rotations[0]) #creates new rather than referenced copy
            newrot.EndAge=0.
            newrot.StartAge=0.
            newrot.EndChron='None'
            newrot.RotPars.RotAng=0.
            self.rotations.insert(0,newrot)
        self.N=len(rotations)
        
    def invert(self, inverting='plate'):
        """
        returns an inverted rotationset by sequentially calling the invert method on each rotation
        By default will also invert the plate codes, set second argument to 'time'
        to invert starting and ending ages for rotation instead.
        """
        return FiniteRotationSet([rotation.invert(inverting) for rotation in self.rotations],self.FixedPlate,self.MovingPlate)
        
    def interpolate(self, ages):
        """ 
        Returns a rotation set of interpolated poles and covariances, for specified list of target ages. 
        Searches for the two bracketing rotations for each target age, following method of Doubrovine and Tarduno (2008). 
        Uses fortran routine ibfrlib0.3.f provided by Pavel Doubrovine, wrapped for python using f2py
        compiled with python wrapped using f2py
        """
        interpolated=[]
        ages.sort() #just in case
        for age in ages:
            age=float(age) #also just in case
            poleages=[item.EndAge for item in self.rotations]
            if age>=min(poleages) and age<=max(poleages):
                #check provided age is within range of given poles
                match=[rotation for rotation in self.rotations if (rotation.EndAge==age)]
                if match:
                    #checks there is a point in interpolating. If a rotation exists with the given age, can just use that.
                    interpolated.append(match[0])
                else: 
                    #get bracketing rotations
                    A1=[rotation for rotation in self.rotations if (rotation.EndAge<age)][-1]
                    A2=[rotation for rotation in self.rotations if (rotation.EndAge>age)][0]
                    xi=(age-A1.EndAge)/(A2.EndAge-A1.EndAge)
                    #do the interpolation by calling finb2       
                    newrot,newcovmat=rotkit_f.finb2(A1.RotPars.tolist(),A1.make_covmat(),A2.RotPars.tolist(),A2.make_covmat(),xi)
                    ndata=(A1.HellingerInfo.Points+A2.HellingerInfo.Points)/2 #not sure if this is justified but need something..
                    #this is admittedly a little clunky, but it does the job.
                    interpolated.append(EulerRotation(pd.Series([A1.MovingPlate,A1.FixedPlate,A1.StartChron,'None',0.,age,
                    newrot[0],newrot[1],newrot[2],
                    1.0,newcovmat[0,0],newcovmat[0,1],newcovmat[0,2],newcovmat[1,1],newcovmat[1,2],newcovmat[2,2],
                    ndata,A1.HellingerInfo.Segs,A1.HellingerInfo.Plates,A1.HellingerInfo.DOF,'Interpolated'],
                    index=['MovingPlate','FixedPlate','StartChron','EndChron','StartAge','EndAge','RotLat','RotLong',
                    'RotAng','Kappahat','a','b','c','d','e','f','Points','Segs','Plates','DOF','Source'])))
        return FiniteRotationSet(interpolated,self.MovingPlate,self.FixedPlate)
                  
    def stagerots(self):
        srots=[]
        #get the stage rotations between each finite rotation in the RotationSet.
        for rot1,rot2 in zip(self.rotations[:-1],self.rotations[1:]):  
            srot=rotmat_to_pole(rot2.make_rotmat()*np.transpose(rot1.make_rotmat()))
            #calculation of covariance matrix for stage rotation between A1 and A2 = A1*(cov2-cov1)
            #from Appendix to Doubrovine & Tarduno 2008
            cov_matrix_s=rot1.make_rotmat()*(rot2.make_covmat()-rot1.make_covmat())
            #not sure that this is totally valid - may be better to just put in smallest number of points?
            points=(rot1.HellingerInfo.Points+rot2.HellingerInfo.Points)/2
            srots.append(EulerRotation(pd.Series([self.MovingPlate,self.FixedPlate,rot1.EndChron,rot2.EndChron,rot1.EndAge,rot2.EndAge,
                    srot[0],srot[1],srot[2],
                    1.0,cov_matrix_s[0,0],cov_matrix_s[0,1],cov_matrix_s[0,2],cov_matrix_s[1,1],cov_matrix_s[1,2],cov_matrix_s[2,2],
                    points,'NA','NA',10000,'Interpolated Stage Pole'],
                    index=['MovingPlate','FixedPlate','StartChron','EndChron','StartAge','EndAge','RotLat','RotLong',
                    'RotAng','Kappahat','a','b','c','d','e','f','Points','Segs','Plates','DOF','Source'])))    
        return StageRotationSet(srots,self.MovingPlate,self.FixedPlate)
        
    def addrots(self,other_rotset,ages=[],sense='to'):
        #adds this rotationset to another one at specified age steps: if no ages provided will interpolate to all finite rotation ages 
        # Default sense is that you are adding 'to' other_rot, such that the fixed plate of other_rot is the fixed plate for the result.
        # If sense is 'on' then the fixed plate of this rotation will be the fixed plate for the result.
        if not ages:
            ages=list(set(self.summary().EndAge.tolist()+other_rotset.summary().EndAge.tolist()))
            ages.sort()
        if sense=='to':
            rotset1,rotset2=self.interpolate(ages),other_rotset.interpolate(ages)
        else:
            rotset1,rotset2=other_rotset.interpolate(ages),self.interpolate(ages)
        added=[]
        for rot1,rot2 in zip(rotset1.rotations,rotset2.rotations):
            added.append(rot1.addrot(rot2))
        return FiniteRotationSet(added,rotset1.MovingPlate,rotset2.FixedPlate) 
        
    def newagemodel(self, timescale='CK95'):
        return FiniteRotationSet([rotation.newagemodel(timescale) for rotation in self.rotations[1:]],self.FixedPlate,self.MovingPlate)  
   
    def summary(self):
        return pd.DataFrame([[item.MovingPlate,item.FixedPlate,item.StartAge,item.EndAge,item.RotPars[0],item.RotPars[1],item.RotPars[2]] for item in self.rotations],
                                columns=['MovingPlate','FixedPlate','StartAge','EndAge','RotLat','RotLong','RotAng'])
                

class StageRotationSet(object):
    """
    a collection of stage EulerRotations with common moving and fixed plates. Methods to easily list set parameters such as ages, chrons,
    generate interpolated rotations and stage rotations
    """
    def __init__(self, rotations,movingplate,fixedplate):
        #where rotations is a list of EulerRotation objects
        #note unlike finite rotationset does not add in a zero rotation.
        self.MovingPlate=movingplate
        self.FixedPlate=fixedplate
        self.rotations=rotations
        self.N=len(rotations)
        
    def invert(self):
        """
        returns an inverted rotationset by sequentially calling the invert method on each rotation
        AFAICT inversion is only meaningful in time for stage poles, so currently that is what it does
        """
        inverted=FiniteRotationSet([rotation.invert('time') for rotation in self.rotations],self.MovingPlate,self.FixedPlate)
        #don't want zeroed rotation for stage poles
        inverted.rotations=inverted.rotations[1:]
        return inverted
                          
    def finiterots(self):
        """
        Adds the stage rotations sequentially to produce net/finite rotations.
        Rotations must be in order.
        """
        finiterots=[self.rotations[0]]
        for rotation in self.rotations[1:]:
            finiterots.append(finiterots[-1].addrot(rotation))
        return FiniteRotationSet(finiterots, self.MovingPlate,self.FixedPlate)
   
    def summary(self):
        return pd.DataFrame([[item.MovingPlate,item.FixedPlate,item.StartAge,item.EndAge,item.RotPars[0],item.RotPars[1],item.RotPars[2]] for item in self.rotations],
                                columns=['MovingPlate','FixedPlate','StartAge','EndAge','RotLat','RotLong','RotAng'])
                                  
class EulerRotation(object):
    """
    Euler rotation for objects on Earth's surface. At the moment it does not make a distinction between finite rotations
    (where one of the age points is 0) and a stage rotation: this may change if different methods are required.
    Attributes:
        MovingPlate: tectonic plate code for moving plate
        FixedPlate: tectonic plate code for fixed plate
        StartAge: starting age of rotation (will normally be 0 for loaded finite rotations)
        StartChron: magnetic chron for starting age ('None' if does not correspond to any chron due to e.g., interpolation)
        EndAge: end age of rotation
        EndChron: magnetic chron for end age ('None' if does not correspond to any chron due to e.g., interpolation)
        RotPars: pandas Series with Latitude, Longitude and Rotation Angle
        Covariances: pandas Series with kappahat and covariance values a,b,c,d,e,f
        HellingerInfo: pandas Series with Points and Segments used for Hellinger fit, 2/3 Plate system, and corresponding degrees of freedom (DOF)
        Source: Published source, or 'Calculated' if created by addition/interpolation etc.
        PlotColor: assigned colour (defaults to black)
    """
    def __init__(self, rotation, plotcolor='black'):
        """Return object
        input is a pandas DataFrame row with neccessary parameters. Should probably try to make it fail gracefully if no covariances?
        the getattr is a way of giving a default value if none is detected. So could potentially do something like that an give default covariances...
        """ 
        self.MovingPlate = rotation.MovingPlate
        self.FixedPlate=rotation.FixedPlate
        self.Timescale=getattr(rotation,'Timescale','Undefined') #potentially useful attribute to know...
        self.StartAge=getattr(rotation,'StartAge',0)
        self.EndAge=rotation.EndAge
        self.StartChron=getattr(rotation,'StartChron','None')
        self.EndChron=rotation.EndChron
        self.RotPars=pd.Series([rotation.RotLat,rotation.RotLong,rotation.RotAng],index=['RotLat','RotLong','RotAng'])
        self.Covariances=pd.Series([rotation.Kappahat,rotation.a,rotation.b,rotation.c,rotation.d,rotation.e,rotation.f],
                                    index=['Kappahat','a','b','c','d','e','f'])
        self.HellingerInfo=pd.Series([rotation.Points,rotation.Segs,rotation.Plates,rotation.DOF], index=['Points','Segs','Plates','DOF'])
        self.Source=rotation.Source
        self.PlotColor=plotcolor

    def make_covmat(self):
        return np.matrix([[self.Covariances.a,self.Covariances.b,self.Covariances.c],
                            [self.Covariances.b,self.Covariances.d,self.Covariances.e],
                            [self.Covariances.c,self.Covariances.e,self.Covariances.f]])
                            
    def make_rotmat(self):
        lat,lon,ang=self.RotPars.RotLat*np.pi/180,self.RotPars.RotLong*np.pi/180,self.RotPars.RotAng*np.pi/180,
        num=1-np.cos(ang)
        px=np.cos(lat)*np.cos(lon)
        py=np.cos(lat)*np.sin(lon)
        pz=np.sin(lat)
        return np.matrix([[(px*px*num )+ np.cos(ang),(px*py*num )-(pz*np.sin(ang)),(px*pz*num)+(py*np.sin(ang))],
                    [(py*px*num)+(pz*np.sin(ang)),(py*py*num)+np.cos(ang),(py*pz*num)-(px*np.sin(ang))],
                    [(px*pz*num)-(py*np.sin(ang)),(pz*py*num)+(px*np.sin(ang)),(pz*pz*num)+np.cos(ang)]])

    def invert(self, inverting='plate'):
        """
        Calculate the inverse of a rotation and its convariance matrix
        By default will also invert the plate codes, set second argument to 'time'
        to invert starting and ending ages for rotation instead.
        """
        # calculate inverse rotation and its covariance matrix - adapted from invrot.f (original November 3 1988 J-Y R)
        # a=b**t; cova=b*covb*(b**t) where b, covb =rotation and covariance matrices of original rotation, and a is the transpose of b
        invrot=EulerRotation(self.details())
        invrot.RotPars.RotAng=-invrot.RotPars.RotAng
        rot_matrix=invrot.make_rotmat()
        cov_matrix=invrot.make_covmat()
        t_rot_matrix=np.transpose(rot_matrix)
        inv_cov=rot_matrix*cov_matrix*t_rot_matrix
        invrot.Covariances=pd.Series([invrot.Covariances.Kappahat,inv_cov[0,0],inv_cov[0,1],inv_cov[0,2],inv_cov[1,1],inv_cov[1,2],inv_cov[2,2]],
                                    index=['Kappahat','a','b','c','d','e','f'])
        if inverting=='time': 
            invrot.StartAge=self.EndAge
            invrot.EndAge=self.StartAge
            invrot.StartChron=self.EndChron
            invrot.EndChron=self.StartChron
        else: 
            invrot.FixedPlate=self.MovingPlate
            invrot.MovingPlate=self.FixedPlate
        return invrot
        
    def half_rotation(self):
        """
        Return a copy of the rotation with the rotation angle halved
        (Useful for ridge reconstructions, although be careful with sense of rotation)
        """        
        half_rot=EulerRotation(self.details())
        half_rot.RotPars.RotAng=half_rot.RotPars.RotAng/2.
        return half_rot
        
    def addrot(self,other_rot,sense='to'):
        #adds together this rotation to another one. 
        # Default sense is that you are adding 'to' other_rot, such that the fixed plate of other_rot is the fixed plate for the result.
        # If sense is 'on' then the fixed plate of this rotation will be the fixed plate for the result.
        # Rotation vector addition adapted from fortran code provided by David Rowley, based on vector addition as described in Cox and Hart  p. 232 and 233
        # Treatment of rotation covariance per Chang 1990 - covariance of combined rotations A and B C=AB, with covariance cova and covb: covc ~= B(t)*covA*B+covB
        # a linear approximation that assumes rotations *associated with the errors* are small enough to be treated as infintesimal so can be added as vectors.
        # Condition for this that det(cova) and det (covb) <<1, which seems to be more than met with typical rotation pole covariances.
        # Adapted from addrot.f (most notes below copied from there) ======  November 3 1988 J-Y R ======= MODIFIED OCTOBER 14, 1991
      
        if sense=='to':
            rot1=EulerRotation(self.details())
            rot2=other_rot
        else:
            rot1=other_rot
            rot2=EulerRotation(self.details())
         
        #1. Check if rotation angles are 0    
        if rot1.RotPars[2]==0 and rot2.RotPars[2]==0: newRotPars=[0,0,0]
        #2. Check if poles are identical (sum to zero)
        elif -rot1.RotPars[2]==rot2.RotPars[2] and rot1.RotPars[0]==rot2.RotPars[0] and rot1.RotPars[1]==rot2.RotPars[1]: newRotPars=[rot1[0],rot1[1],0]
        elif -rot1.RotPars[2]==rot2.RotPars[2] and rot1.RotPars[0]==-rot2.RotPars[0] and rot1.RotPars[1]==rot2.RotPars[1]+180: newRotPars=[rot1[0],rot1[1],0]
        elif -rot1.RotPars[2]==rot2.RotPars[2] and rot1.RotPars[0]==-rot2.RotPars[0] and rot1.RotPars[1]==rot2.RotPars[1]-180: newRotPars=[rot1[0],rot1[1],0]
        else: newRotPars=rotmat_to_pole(rot2.make_rotmat()*rot1.make_rotmat())
        
        #now to deal with covariances. Currently, the kappa estimation when DOF<50 for one or both kappahat is not implemented.
        new_covmat=(np.transpose(rot2.make_rotmat())*(rot2.make_covmat()/rot2.Covariances.Kappahat)
                                *rot2.make_rotmat())+(rot1.make_covmat()/rot1.Covariances.Kappahat)
        newKappa=1
        newDOF=10000
        
        return EulerRotation(pd.Series([rot1.MovingPlate,rot2.FixedPlate,self.Timescale,rot1.StartChron,rot2.EndChron,rot1.StartAge,rot2.EndAge,
                    newRotPars[0],newRotPars[1],newRotPars[2],
                    newKappa,new_covmat[0,0],new_covmat[0,1],new_covmat[0,2],new_covmat[1,1],new_covmat[1,2],new_covmat[2,2],
                    rot1.HellingerInfo.Points+rot2.HellingerInfo.Points,'NA','NA',newDOF,'Calculated'],
                    index=['MovingPlate','FixedPlate','Timescale','StartChron','EndChron','StartAge','EndAge','RotLat','RotLong',
                    'RotAng','Kappahat','a','b','c','d','e','f','Points','Segs','Plates','DOF','Source']))             
    
    def rotate(self,rotation):
        """
        Rotates rotation pole by EulerRotation rotation (i.e. into a different reference frame)
        Most relevent for stage rotations: makes the FixedPlate/MovingPlate a bit tricky potentially.
        """ 
        result=rotkit_f.rotatepts(np.array([[self.RotPars.RotLat,self.RotPars.RotLong]]),
                np.array([rotation.RotPars.tolist()+rotation.Covariances.tolist()+[rotation.EndAge]]),1)
        return EulerRotation(pd.Series([self.MovingPlate,self.FixedPlate,self.Timescale,self.StartChron,self.EndChron,self.StartAge,self.EndAge,result[:,0][0],result[:,1][0],self.RotPars[2],
                        self.Covariances.Kappahat,self.Covariances.a,self.Covariances.b,self.Covariances.c,self.Covariances.d,self.Covariances.e,self.Covariances.f,
                        self.HellingerInfo.Points,self.HellingerInfo.Segs,self.HellingerInfo.Plates,self.HellingerInfo.DOF,self.Source],
                            index=['MovingPlate','FixedPlate','Timescale','StartChron','EndChron','StartAge','EndAge','RotLat','RotLong','RotAng',
                            'Kappahat','a','b','c','d','e','f','Points','Segs','Plates','DOF','Source']))
    
    def details(self):
        return pd.Series([self.MovingPlate,self.FixedPlate,self.Timescale,self.StartChron,self.EndChron,self.StartAge,self.EndAge,self.RotPars[0],self.RotPars[1],self.RotPars[2],
                        self.Covariances.Kappahat,self.Covariances.a,self.Covariances.b,self.Covariances.c,self.Covariances.d,self.Covariances.e,self.Covariances.f,
                        self.HellingerInfo.Points,self.HellingerInfo.Segs,self.HellingerInfo.Plates,self.HellingerInfo.DOF,self.Source],
                            index=['MovingPlate','FixedPlate','Timescale','StartChron','EndChron','StartAge','EndAge','RotLat','RotLong','RotAng',
                            'Kappahat','a','b','c','d','e','f','Points','Segs','Plates','DOF','Source'])

    def newagemodel(self, timescale='CK95'):
        """
        allows polarity timescale switching. Defaults to CK95
        """
        if self.EndChron != 'None': NewEndAge=get_chron_age(self.EndChron, timescale)
        else: NewEndAge=self.EndAge
        #hopefully this will work for stage rotations where StartChron is specified
        if self.StartChron != 'None': NewStartAge=get_chron_age(self.StartChron, timescale)
        else: NewStartAge=self.StartAge
        
        return EulerRotation(pd.Series([self.MovingPlate,self.FixedPlate,timescale,self.StartChron,self.EndChron,NewStartAge,NewEndAge,self.RotPars[0],self.RotPars[1],self.RotPars[2],
                        self.Covariances.Kappahat,self.Covariances.a,self.Covariances.b,self.Covariances.c,self.Covariances.d,self.Covariances.e,self.Covariances.f,
                        self.HellingerInfo.Points,self.HellingerInfo.Segs,self.HellingerInfo.Plates,self.HellingerInfo.DOF,self.Source],
                            index=['MovingPlate','FixedPlate','Timescale','StartChron','EndChron','StartAge','EndAge','RotLat','RotLong','RotAng',
                            'Kappahat','a','b','c','d','e','f','Points','Segs','Plates','DOF','Source']))
        
                            
    def summary(self):
        return pd.Series([self.MovingPlate,self.FixedPlate,self.StartAge,self.EndAge,self.RotPars[0],self.RotPars[1],self.RotPars[2]],
                            index=['MovingPlate','FixedPlate','StartAge','EndAge','RotLat','RotLong','RotAng'])   

class Point(object):
    """ Baseclass for a point that can be acted on by rotations
    
        Attributes:
        PointName: string describing feature
        PlateCode: tectonic plate code on which point(s) are located
        LocPars: pandas Series with Latitude and Longitude, plus rotation error ellipse parameters
        FeatureAge: age of feature
        ReconstructionAge: age that current set of points has been reconstructed at
        PlotColor: assigned colour
        PlotLevel: plot level (defaults to 5)
    """    
    def __init__(self, PointPars,PlotColor='grey',PlotLevel=5,PlotSymbolSize=2):
        """Return object
        PointPars should be a Series/DataFrame row with Name,PlateCode,Lat,Lon,FeatureAge,ReconstructionAge;
        optionally MaxError,MinError,MaxBearing, otherwise 0 by default
        """ 
        self.PointName = PointPars.Name
        self.PlateCode=PointPars.PlateCode
        #not sure how much, but information about reference frame could be useful
        self.ReferencePlate=PointPars.PlateCode
        self.FeatureAge=PointPars.FeatureAge
        self.ReconstructionAge=PointPars.ReconstructionAge
        self.PlotColor=PlotColor
        self.PlotLevel=PlotLevel
        self.PlotSymbolSize=PlotSymbolSize
        #if input does not contain rotation error parameters, creates empty columns 
        if not 'MaxError' in PointPars:
            self.LocPars=pd.Series([PointPars.Lat,PointPars.Lon,0.,0.,0.],index=['PointLat','PointLong','MaxError','MinError','MaxBearing'])
        else:
            self.LocPars=pd.Series([PointPars.Lat,PointPars.Lon,PointPars.MaxError,PointPars.MinError,PointPars.MaxBearing]
                                        ,index=['PointLat','PointLong','MaxError','MinError','MaxBearing'])

    
    def mapplot(self,ellipseflag=0,pointcol=""):
        """
        plots point on preexisting axes (currently set up for Cartopy Geoaxes) 
        If ellipseflag set to 1, will plot the associated error ellipse
        """
        if pointcol=="": pointcol=self.PlotColor
        plt.plot(self.LocPars.PointLong,self.LocPars.PointLat, marker='o',
                ms=self.PlotSymbolSize, color=pointcol, zorder=self.PlotLevel,
                transform=ccrs.Geodetic())
        ax=plt.gca()
        if ellipseflag==1: 
            ax.add_patch(ellipse_pars(self.LocPars.PointLong, self.LocPars.PointLat, self.LocPars.MaxError,self.LocPars.MinError,self.LocPars.MaxBearing,
                        facecolor='none', edgecolor=pointcol, zorder=self.PlotLevel-1))

    def rotate(self,rotation):
        """Rotates pointset by EulerRotation rotation
        """ 
        #intially tried to make a copy and modify it, but it seems you can call the object type to initialise a new version
        result=rotkit_f.rotatepts(np.array([[self.LocPars.PointLat,self.LocPars.PointLong]]),
                np.array([rotation.RotPars.tolist()+rotation.Covariances.tolist()+[rotation.EndAge]]),1)
                
        rotated=Point(pd.Series([self.PointName,self.PlateCode,self.FeatureAge,rotation.EndAge,result[:,0][0],result[:,1][0],result[:,3][0],result[:,4][0],result[:,5][0]],
                                index=['Name','PlateCode','FeatureAge','ReconstructionAge','Lat','Lon','MaxError','MinError','MaxBearing']))
        rotated.ReferencePlate=rotation.FixedPlate
        return rotated
        
    def reconstruct(self,age,refplate,rotmodel):
        """
        Finds the EulerRotation for specified age, then rotates it
        """
        #first check if the reference plate is this plate.
        if refplate==self.PlateCode:
            return self
        #another place where adding the zero rotation could trip you up.
        else: return self.rotate(rotmodel.get_rots(self.PlateCode,refplate,[age]).rotations[-1])
        
    def flowline(self,ages,refplate,rotmodel,SetName='Flowline',PlotColor='grey',PlotLevel=5):
        """
        Finds coordinates of line that charts motion of point relative to the reference plate for the specified age points
        Current output: list of Point objects. Eventually will be a flowline object
        """
        return Flowline(pd.DataFrame([self.reconstruct(age,refplate,rotmodel).summary() for age in ages]),refplate,SetName,PlotColor,PlotLevel)
  
    def motion_vector(self,refplate,rotmodel,startage=0.78,endage=0,fixed=0):
        """
        calculates the magnitude and direction of the plate motion vector for this locality relative to the specified
        reference plate. By default, it will calculate the most recent interval (last 0.78 Ma of motion).
        outputs a bearing and a rate in mm/yr or km/Myr; also N-S and E-W components of velocity and total displacement.
        set fixed=1 to show vector for refplate relative to point 
        """
        #velocity of point on Earth's surface=cross product of Euler vector (omega) and position vector (r)
        #N-S component of v vNS = a*|rotrate|*cos(polelat)*sin(sitelong-polelong)
        #E-W component of v vEW =  a*|rotrate|* [cos(sitelat)*sin(polelat)-sin(sitelat)*cos(polelat)*cos(sitelong-polelong)]
        #where a=Earth radius
        #rate of motion = sqrt (vNS^2+vEW^2)
        #azimuth = 90-atan[vNS/vEW]

        #get Lat and Long of starting position of point in relevant reference frame
        #if start age or endage is 0, causes problems... (because if not there is an additional 0 Ma rotation)
        if startage==0 or endage==0: rotindex=0
        else: rotindex=1
        
        reconstruction_rots=rotmodel.get_rots(self.PlateCode,refplate,[startage,endage])
        latlong=self.rotate(reconstruction_rots.rotations[rotindex]).LocPars   
        pointlat=latlong.PointLat*np.pi/180
        pointlong=latlong.PointLong*np.pi/180
        
        #get the stage rotation over the relevant interval, and invert if going forward in time (normally will be)
        stagerot=reconstruction_rots.stagerots().rotations[rotindex]
        if startage>endage:
            stagerot=stagerot.invert('time')
        if fixed==1:
            stagerot=stagerot.invert()        
        rotlat=stagerot.RotPars.RotLat*np.pi/180
        rotlong=stagerot.RotPars.RotLong*np.pi/180
        rotrate=stagerot.RotPars.RotAng*np.pi/180
        vNS=6371.*rotrate*np.cos(rotlat)*np.sin(pointlong-rotlong)
        vEW=6371.*rotrate*(np.cos(pointlat)*np.sin(rotlat)-np.sin(pointlat)*np.cos(rotlat)*np.cos(pointlong-rotlong))
        azimuth=90.-(180/np.pi*np.arctan2(vNS,vEW))
        if azimuth<0.: azimuth=azimuth+360.
        return pd.Series([stagerot.StartAge,stagerot.EndAge,np.sqrt(vNS**2+vEW**2)/abs(stagerot.StartAge-stagerot.EndAge),
                            azimuth,np.sqrt(vNS**2+vEW**2),vNS,vEW],
                            index=['StartAge','EndAge','Rate','Bearing','Distance','NS_component','EW_component'])
                            
    def predict_DI(self,ages,rotmodel,abs_ref_frame=3):
        """
        Predicts the Declination and Inclination that should be observed from samples of age ages
        at this site - a range of ages can be used because of the possibility of remagnetisation.
        absolute_ref_frame should be an absolute or hotspot frame of reference (e.g. Pacific=3)
        """
        VGPs=rotmodel.synthetic_APWP(self.PlateCode,abs_ref_frame,ages)
        reconstruction_rots=rotmodel.get_rots(self.PlateCode,abs_ref_frame,ages)
        rotated=[self.rotate(rotation) for rotation in reconstruction_rots.rotations[1:]]
        paleoI=[np.arctan(2*np.tan(point.LocPars.PointLat*np.pi/180))*180/np.pi for point in rotated]
        pp=np.array(([np.sin((90-vgp_lat)*np.pi/180) for vgp_lat in VGPs.Lat]))
        dphi=np.array(([np.sin((vgp_lon-self.LocPars.PointLong)*np.pi/180) for vgp_lon in VGPs.Lon]))
        pm=np.array(([np.sin((90-point.LocPars.PointLat)*np.pi/180) for point in rotated]))
        paleoD=np.arcsin(pp*dphi/pm)*180/np.pi
        return pd.DataFrame(np.column_stack((ages,[point.LocPars.PointLat for point in rotated],[point.LocPars.PointLong for point in rotated],paleoD,paleoI)), 
                                        columns=['Age','Lat','Lon','PredDec','PredInc'])
        
    def summary(self):
        return pd.Series([self.PointName,self.PlateCode,self.FeatureAge,self.ReconstructionAge]+self.LocPars.tolist(),
                                index=['Name','PlateCode','FeatureAge','ReconstructionAge','Lat','Lon','MaxError','MinError','MaxBearing'])       
        
class PointSet(object):
    """ Baseclass for a set of Points assigned to the same plate that can be operated on. 
        Methods match those for Point: basically loop through them
        Mapplot plots as individual points.
    
        Attributes:
        SetName: string describing feature
        PlateCode: tectonic plate code on which points are located
        
        Other Point attributes (e.g. FeatureAge, colors may vary by point.
    """
    def __init__(self,PointList,SetName='PointSet',PlotColor='grey',PlotLevel=5):
        """Return object
        PointList should be a DataFrame with columns name,PlateCode,Lat,Lon,FeatureAge,ReconstructionAge;
        optionally MaxError,MinError,MaxBearing
        """
        self.points=[Point(point,PlotColor,PlotLevel) for i,point in PointList.iterrows()]
        self.SetName = SetName
        self.PlateCode=PointList.iloc[0].PlateCode
        #not sure how much, but information about reference frame could be useful
        self.ReferencePlate=self.PlateCode
        self.FeatureAge=PointList.iloc[0].FeatureAge
        self.ReconstructionAge=PointList.iloc[0].ReconstructionAge

    def mapplot(self,ellipseflag=0):
        for point in self.points:
            point.mapplot(ellipseflag)
    
    def rotate(self,rotation):
        """Rotates pointset by EulerRotation rotation
        """ 
        rotated=copy.deepcopy(self)
        #this is a lot less fiddly than converting rotated point list into format where can create PointSet de novo
        rotated.points=[point.rotate(rotation) for point in self.points]
        rotated.ReferencePlate=rotation.FixedPlate
        rotated.ReconstructionAge=rotation.EndAge
        return rotated
    
    def reconstruct(self,age,refplate,rotmodel):
        """
        Finds the EulerRotation for specified age, then rotates it
        Todo: drop points where FeatureAge<ReconstructionAge?
        """
        #first check if the reference plate is this plate.
        if refplate==self.PlateCode:
            return self
        #another place where adding the zero rotation could trip you up.
        else: return self.rotate(rotmodel.get_rots(self.PlateCode,refplate,[age]).rotations[-1])
  
    def motion_vectors(self,refplate,rotmodel,age_range=[1,0]):
        """
        calculates the magnitudes and directions of the plate motion vector for each locality relative to the specified
        reference plate. By default, it will calculate the contemporary vector (last 1 Ma of motion).
        outputs bearings and rates in mm/yr or km/Myr; also N-S and E-W components of velocity and total displacement.
        """
        return pd.DataFrame([point.motion_vector(refplate,rotmodel,age_range) for point in self.points])

    def summary(self):
        return pd.DataFrame([[point.PointName,point.PlateCode,point.FeatureAge,point.ReconstructionAge]+point.LocPars.tolist() for point in self.points],
                            columns=['Name','PlateCode','FeatureAge','ReconstructionAge','Lat','Lon','MaxError','MinError','MaxBearing'])
  
               
class Boundary(PointSet):        
    """ line that can be acted on by rotations
    
        SetName: string describing feature
        PlateCode: tectonic plate code on which points are located
        ReconstructionAge: age that current set of points has been reconstructed at
        points: List of Points
    """
    def __init__(self,PointList,SetName='PointSet',PlotColor='grey',PlotLevel=5):
        """Return object
        PointList should be a DataFrame with columns name,PlateCode,Lat,Lon,FeatureAge,ReconstructionAge;
        optionally MaxError,MinError,MaxBearing
        """
        PointSet.__init__(self,PointList,SetName,PlotColor,PlotLevel)
        self.PlotColor=PlotColor
        self.PlotLevel=PlotLevel

        
    def mapplot(self,thickness=2):          
        plt.plot(self.summary().Lon.values,self.summary().Lat.values, 
                linewidth=thickness, color=self.PlotColor, zorder=self.PlotLevel,
                transform=ccrs.Geodetic())

class Platelet(PointSet):
    """ closed Polygon that can be acted on by rotations but no defined subsegments as PlatePolygon
    
        Attributes (inherited from PointSet):
        name: string describing feature
        PlateCode: tectonic plate code on which point(s) are located
        FeatureAge: age of feature
        ReconstructionAge: age that current set of points has been reconstructed at
        LatLons: pandas DataFrame with Lat and Long points 
        PlotColor: assigned colour (defaults to grey)
        PlotLevel: plot level (defaults to 5)
    """
    def __init__(self,PointList,SetName='PointSet',PlotColor='grey',PlotLevel=5):
        """Return object
        PointList should be a DataFrame with columns name,PlateCode,Lat,Lon,FeatureAge,ReconstructionAge;
        optionally MaxError,MinError,MaxBearing
        """
        PointSet.__init__(self,PointList,SetName,PlotColor,PlotLevel)
        self.PlotColor=PlotColor
        self.PlotLevel=PlotLevel

    def mapplot(self,m,thickness=2,transparency=0.5):
        """plot polygon on specified basemap"""
        self.pltx,self.plty=self.summary().Lon.values,self.summary().Lat.values
        self.polygon=Polygon(zip(self.pltx,self.plty),
                             facecolor=self.PlotColor, alpha=transparency, zorder=self.PlotLevel,
                             transform=ccrs.Geodetic())
        plt.gca().add_patch(self.polygon)
        plt.plot(self.pltx,self.plty, color=self.PlotColor, linewidth=thickness,zorder=self.PlotLevel,transform=ccrs.Geodetic())

        
class Flowline(PointSet):
    """
    A set of points tracking a feature over a range of reconstruction ages
    Attributes:      
    name: string describing feature
    
    MovingPlate: tectonic plate code for point being modelled
    FixedPlate: reference frame for point rotation
    PlotLevel: plot level (defaults to 5)
    
    Methods:
    mapplot: 
    mapplot_age:
    summary:
    """
    def __init__(self,PointList,RefPlate,SetName='PointSet',PlotColor='grey',PlotLevel=5):
        PointSet.__init__(self,PointList,SetName,PlotColor,PlotLevel)
        self.PlotColor=PlotColor
        self.PlotLevel=PlotLevel
        self.MovingPlate=PointList.iloc[0].PlateCode
        self.FixedPlate=RefPlate
        
    def mapplot(self,thickness=2):          
        plt.plot(self.summary().Lon.values,self.summary().Lat.values, 
                linewidth=thickness, color=self.PlotColor, zorder=self.PlotLevel,
                transform=ccrs.Geodetic())
                            
    def mapplot_age(self,colourmap='plasma_r',plotbar='N',ellipseflag=1):
        age_cmap=plt.get_cmap(colourmap)
        maxage=np.round(max(self.summary().ReconstructionAge),-1)
        minage=np.round(min(self.summary().ReconstructionAge),-1)
        cNorm  = colors.Normalize(vmin=minage, vmax=maxage)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=age_cmap)
        i=0
        x0,y0,age0=0.,0.,0.
        for point in self.points:
            point.mapplot(ellipseflag,pointcol=scalarMap.to_rgba(point.ReconstructionAge))
            if i>0:
                plt.plot([x0,point.LocPars.PointLong], [y0,point.LocPars.PointLat],
                        linewidth=2,
                        color=scalarMap.to_rgba(point.ReconstructionAge+((point.ReconstructionAge-age0)/2)),
                        transform=ccrs.Geodetic(),
                        zorder=self.PlotLevel-1)
            x0,y0,age0=point.LocPars.PointLong,point.LocPars.PointLat,point.ReconstructionAge
            i=i+1 
            
        #add colorbar for age
        if plotbar=='Y':   
            ax1 = plt.gca()
            #sneaky way to extend axes with right dimensions
            divider=make_axes_locatable(ax1)
            ax2 = divider.append_axes("bottom", size="5%", pad=0.1)
            cb1=colorbar.ColorbarBase(ax2, cmap=age_cmap,norm=cNorm,orientation='horizontal')
            cb1.set_label('Age (Ma)')                


class APWP(PointSet):
    """
    A flowline specifically for APWPs, which unlike the more standard flowline are a reconstruction of time 0
    of the position of the VGP at a particular age. Currently set up to 
    Attributes:      
    name: string describing feature
    PlateCode: tectonic plate code for point being modelled
    FixedPlate: reference frame for point rotation
    PlotLevel: plot level (defaults to 5)
    
    Methods:
    plot
    """
    def __init__(self,PointList,RefPlate,SetName='PointSet',PlotColor='grey',PlotLevel=5):
        PointSet.__init__(self,PointList,SetName,PlotColor,PlotLevel)
        self.PlotColor=PlotColor
        self.PlotLevel=PlotLevel
        self.MovingPlate=PointList.iloc[0].PlateCode
        self.FixedPlate=RefPlate
        
    def mapplot(self,thickness=2):          
        plt.plot(self.summary().Lon.values,self.summary().Lat.values, 
                linewidth=thickness, color=self.PlotColor, zorder=self.PlotLevel,
                transform=ccrs.Geodetic())
                            
    def mapplot_age(self,colourmap='plasma_r',plotbar='N',ellipseflag=1):
        age_cmap=plt.get_cmap(colourmap)
        maxage=np.round(max(self.summary().FeatureAge),-1)
        minage=np.round(min(self.summary().FeatureAge),-1)
        cNorm  = colors.Normalize(vmin=minage, vmax=maxage)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=age_cmap)
        i=0
        x0,y0,age0=0.,0.,0.
        for point in self.points:
            point.mapplot(ellipseflag,pointcol=scalarMap.to_rgba(point.FeatureAge))
            if i>0:
                plt.plot([x0,point.LocPars.PointLong], [y0,point.LocPars.PointLat],
                        linewidth=2,
                        color=scalarMap.to_rgba(point.FeatureAge+((point.FeatureAge-age0)/2)),
                        transform=ccrs.Geodetic(),
                        zorder=self.PlotLevel-1)
            x0,y0,age0=point.LocPars.PointLong,point.LocPars.PointLat,point.FeatureAge
            i=i+1                         
                                                
                                                            
                                        
class AMS_Locality(Point):
    """ A sampling site with associated AMS data 
    New Attributes:
    AMS_info: a pandas series that contains AMS ellipsoid information (currently just lineation orientation and error)
    Age_err: Age Error associated with site
    """
    def __init__(self, AMSpars,PlotColor='grey',PlotLevel=5,PlotSymbolSize=2):
        Point.__init__(self,AMSpars,PlotColor,PlotLevel,PlotSymbolSize)
        #need to add whole ellispoid info here eventually...
        self.AMS_info=pd.Series([AMSpars.AMS_max,AMSpars.max_err],index=['k_max','k_max_err'])
        self.Age_err=AMSpars.Age_err
        
    def mapplot(self,m,plottrend=0):
        """
        plottrend=1 will add the trend of the sigma1 direction to the symbol.
        """
        self.pltx,self.plty=self.LatLons.Lon,self.LatLons.Lat
        plt.plot(self.pltx,self.plty, marker='o',ms=self.PlotSymbolSize, color=self.PlotColor, zorder=self.PlotLevel,transform=ccrs.Geodetic())
        if plottrend==1:
            trend=self.sigma_1()
            plt.quiver(self.pltx,self.plty,np.sin(trend*np.pi/180.),np.cos(trend*np.pi/180.), pivot='tip',color=self.PlotColor,zorder=self.PlotLevel,transform=ccrs.Geodetic())
            plt.quiver(self.pltx,self.plty,-np.sin(trend*np.pi/180.),-np.cos(trend*np.pi/180.), pivot='tip',color=self.PlotColor,zorder=self.PlotLevel,transform=ccrs.Geodetic())
   
    def sigma_1(self, quadrant='E'):
        s1_dir=self.AMS_info.k_max+90.
        a,b=np.cos(s1_dir*np.pi/180.), np.sin(s1_dir*np.pi/180.)
        if quadrant=='N':
            if a<0: s1_dir=s1_dir-180.
            elif b<0: s1_dir=s1_dir-360.
        elif quadrant=='S':
            if a>0: s1_dir=s1_dir+180.
        elif quadrant=='W':
            if b>0: s1_dir=s1_dir+180.
        elif quadrant=='E':
            if b<0: s1_dir=s1_dir-180. 
        if s1_dir>360.: s1_dir=s1_dir-360   
        return s1_dir
        
    def rotate(self,rotation):
        """Rotates pointset by EulerRotation rotation
        Also rotates the AMS_max direction: doesn't presently incorporate rotation errors
        """ 
        #to rotate the fabric eigenvector declinations,need to have some kind of reference.
        ref_lat,ref_long=sphere_point_along_bearing(self.LocPars.PointLat,self.LocPars.PointLong,self.AMS_info.k_max,0.1)
        result=rotkit_f.rotatepts(np.column_stack(([self.LocPars.PointLat,ref_lat],[self.LocPars.PointLong,ref_long])),
                np.array([rotation.RotPars.tolist()+rotation.Covariances.tolist()+[rotation.EndAge]]),1)
        #by calculating the bearing to the rotated reference point can rotate the eigenvector dec               
        rotated=AMS_Locality(pd.Series([self.PointName,self.PlateCode,self.FeatureAge,self.Age_err,rotation.EndAge,
                                    result[:,0][0],result[:,1][0],result[:,3][0],result[:,4][0],result[:,5][0],
                                    sphere_bearing(result[:,0][0],result[:,1][0],result[:,0][1],result[:,1][1]),self.AMS_info.k_max_err],
                                index=['Name','PlateCode','FeatureAge','Age_err','ReconstructionAge','Lat','Lon','MaxError','MinError','MaxBearing','AMS_max','max_err']),
                                self.PlotColor,self.PlotLevel,self.PlotSymbolSize)
        rotated.ReferencePlate=rotation.FixedPlate
        return rotated

    def strain_history(self,rotmodel,ages,moving_plate,changeref=-1):
        """
        Reconstructs plate vector directions over ages according to the rotation model
        Assumes this corresponds to maximum strain axis for region.
        change_ref option allows you to look at plate motion vector relative to the fixed plate
        inboard of the deforming region, if have specified small block codes; set it to the relevent platecode.
        Note: assumes that you want forward motion, so will reverse age order.
        """
        ages.reverse()
        if changeref>0:
            ActualPlate=self.PlateCode
            self.PlateCode=changeref
        result=pd.DataFrame([self.motion_vector(moving_plate,rotmodel,age1,age2,1) for age1,age2 in zip (ages[:-1],ages[1:])])
        if changeref>0:
            self.PlateCode=ActualPlate
        return result

class PMag_Locality(Point):
    """ A sampling site with associated AMS data 
    New Attributes (in addition to those inherited from Point):
    Age_err: Age error associated with locality
    DI_info: a pandas series that contains remanence vector information (GeoD,GeoI,TiltD,TiltI,k,alpha95)
    Tilt_info: a pandas series that contains structural information (DipAzimuth,Dip and optionally FoldMinAge,FoldMaxAge)
    (Use of Azimuth/Dip rather than strike/dip because that is input to dotilt(), and don't need to worry about different conventions.
    """
    def __init__(self, PMag_data,PlotColor='grey',PlotLevel=5,PlotSymbolSize=2):
        """Return object
        
        """ 
        Point.__init__(self,PMag_data,PlotColor,PlotLevel,PlotSymbolSize)
        self.Age_err=PMag_data.Age_err
        self.DI_info=pd.Series([PMag_data.GeoD,PMag_data.GeoI,PMag_data.TiltD,PMag_data.TiltI,PMag_data.k,PMag_data.alpha95],
                                index=['GeoD','GeoI','TiltD','TiltI','k','alpha95'])
        self.Tilt_info=pd.Series([PMag_data.DipAzimuth,PMag_data.Dip,PMag_data.FoldMaxAge,PMag_data.FoldMinAge],
                                    index=['DipAzimuth','Dip','FoldMaxAge','FoldMinAge']) 
                                    
    def vgp(self, tiltcorr=1.):
        """
        Calculates VGP for locality. By default assumes that the tilt-corrected remanence is the
        correct one.
        """
        if tiltcorr==1:
            dec=self.DI_info.TiltD*np.pi/180
            inc=self.DI_info.TiltI*np.pi/180
        elif tiltcorr==0:
            dec=self.DI_info.GeoD*np.pi/180
            inc=self.DI_info.GeoI*np.pi/180           
        else: 
            direction=self.untilt(tiltcorr)
            dec=direction.Dec*np.pi/180
            inc=direction.Inc*np.pi/180
        colat=np.arctan2(inc,2)    
        vgplat=np.arcsin((np.sin(self.LocPars.PointLat*np.pi/180)*np.cos(colat))+(np.cos(self.LocPars.PointLat*np.pi/180)*np.sin(colat)*np.cos(dec)))
        beta=np.arcsin((np.sin(colat)*np.sin(dec))/np.cos(vgplat))
        if np.cos(colat)>=np.sin(self.LocPars.PointLat)*np.sin(vgplat):
            vgplong=self.LocPars.PointLong*np.pi/180+beta
        else:
            vgplong=self.LocPars.PointLong*np.pi/180+np.pi-beta
        dp=self.DI_info.alpha95*np.pi/180*((1+3*np.cos(colat)**2)/2)
        dm=self.DI_info.alpha95*np.pi/180*(np.sin(colat)/np.cos(inc))    
        #eventually am going to want to create a VGP object for this
        return pd.Series([vgplat*180/np.pi,vgplong*180/np.pi,dp*180/np.pi,dm*180/np.pi],index=['PoleLat','PoleLong','dp','dm'])

    def untilt(self,untilting=1):
        """
        performs a partial tilt correction according to untilting (0=geographic, 1=stratigraphic)
        adapted from the dotilt routine in PMagPy (Lisa Tauxe) 
        """
        rad=np.pi/180. # converts from degrees to radians
        X=dir2cart([self.DI_info.GeoD,self.DI_info.GeoI,1.]) # get cartesian coordinates of dec,inc
    # get some sines and cosines of new coordinate system
        partial_dip=self.Tilt_info.Dip*untilting
        sa,ca= -np.sin(self.Tilt_info.DipAzimuth*rad),np.cos(self.Tilt_info.DipAzimuth*rad)
        cdp,sdp= np.cos(partial_dip*rad),np.sin(partial_dip*rad)
    # do the rotation
        xc=X[0]*(sa*sa+ca*ca*cdp)+X[1]*(ca*sa*(1.-cdp))+X[2]*sdp*ca
        yc=X[0]*ca*sa*(1.-cdp)+X[1]*(ca*ca+sa*sa*cdp)-X[2]*sa*sdp
        zc=X[0]*ca*sdp-X[1]*sdp*sa-X[2]*cdp
    # convert back to direction:
        Dir=cart2dir([xc,yc,-zc])
        return pd.Series([untilting,Dir[0],Dir[1]],index=['Untilting','Dec','Inc'])

    def untilt_history(self):
        """
        returns the D,I at 1 Ma steps between formation age and present, with higher resolution steps during the
        interval defined by FoldMinAge and FoldMaxAge
        """
        #Assumes the folding is linear over the specified folding interval, which assumes that the folding age is well
        # defined. Will be a problem if not, so will need more options.
        # Probably also at some stage want to make it able to return single age or specific age range..? 
        tiltages=np.concatenate(((np.arange(0,self.Tilt_info.FoldMinAge)),
                            np.arange(self.Tilt_info.FoldMinAge,self.Tilt_info.FoldMaxAge,0.1),
                            np.arange(self.Tilt_info.FoldMaxAge,self.FeatureAge+self.Age_err),
                            np.array([self.FeatureAge+self.Age_err])))
        #set untilting to a linear function of position within folding interval
        tiltcorr=(tiltages-self.Tilt_info.FoldMinAge)/(self.Tilt_info.FoldMaxAge-self.Tilt_info.FoldMinAge)
        #set untilting to 0 for age points less than Fold Minimum Age
        tiltcorr[[i for i, age in enumerate(tiltages) if age <= self.Tilt_info.FoldMinAge]]=0.
        #set untilting to 1 for age points greater than Fold Maximum Age
        tiltcorr[[i for i, age in enumerate(tiltages) if age>=self.Tilt_info.FoldMaxAge]]=1
        tiltDI=[]
        for age, untilt in zip (tiltages,tiltcorr):
            if untilt==0: tiltDI.append ([age,untilt,self.DI_info.GeoD,self.DI_info.GeoI])
            elif untilt==1: tiltDI.append ([age,untilt,self.DI_info.TiltD,self.DI_info.TiltI])
            else: 
                corrected=self.untilt(untilt)
                tiltDI.append ([age,untilt,corrected.Dec,corrected.Inc])
        return pd.DataFrame(tiltDI, columns=['Age','Untilting','Dec','Inc'])
        
    def compare_to_predicted(self,rotmodel,abs_ref_frame=3):
        """
        between the maximum locality age and the present, compare the D,I values predicted from the specified rotation model to
        the geographic, synfolding and tilt-corrected D,I (depending on age)
        """
        actual_DI=self.untilt_history()
        predicted_DI=self.predict_DI(actual_DI.Age.values.tolist(),rotmodel,abs_ref_frame)
        return pd.DataFrame(np.column_stack((actual_DI.Age,actual_DI.Dec,actual_DI.Inc,predicted_DI.PredDec,predicted_DI.PredInc,
                                    angle([[row.Dec,row.Inc] for i,row in actual_DI.iterrows()],[[row.PredDec,row.PredInc] for i,row in predicted_DI.iterrows()]))),
                                        columns=['Age','Dec','Inc','PredDec','PredInc','Seperation_Angle'])





