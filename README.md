## What is RotKit?
RotKit is a python-based set of tools and routines for manipulating Euler rotations and using them for reconstructing plate motions and related tectonic parameters.

## Current Status
- As of March 2019, many of the core routines are working, stability of code (or said code’s behaviour) is not guaranteed! 
- Also note that some functions are written in Fortran and have been linked to python using f2py. For these functions to work on a new machine, RotKit.f and RotKit-nonpython.f need to be compiled. The instructions listed in RotKit_compile.txt are currently very specific to my own set-up (compiled using gfortran on a Mac running OSX 10.10) but may provide some guidance.

# EulerRots.py
Basic summary: an object oriented framework where various rotation objects (Rotations bundled into RotationSets bundled into a RotationModel) can be manipulated and used to reconstruct various point, line and polygon objects on the Earth's surface.

## Individual functions

#### rotmat_to_pole(rot_matrix)
Given a 3x3 rotation matrix rot_matrix, converts it into a Euler rotation (lat,long,angle)

### Spherical Trigonometry Functions

#### sphere_ang_dist(lat1,long1,lat2,long2,degrees=True):
Calculates the length of the great circle path between points (lat1, long1) and (lat2, long2) using haversine formula. By default, assumes input coordinates are in degrees and returns distance in degrees. If degrees=False, then assumes input is in radians and returns distance in radians.

Source: [http://www.movable-type.co.uk/scripts/latlong.html](http://www.movable-type.co.uk/scripts/latlong.html)

#### sphere_bearing(lat1,long1,lat2,long2,degrees=True):
Returns bearing of point (lat2, long2) wrt point (lat1,long1). By default, assumes input coordinates are in degrees and returns distance in degrees. If degrees=False, then assumes input is in radians and returns distance in radians.
    
Source: [http://www.movable-type.co.uk/scripts/latlong.html](http://www.movable-type.co.uk/scripts/latlong.html)

#### sphere_dist_bearing(lat1,long1,lat2,long2,degrees=True)
Returns a list containing angular distance and bearing of point (lat2, long2) wrt point (lat1,long1). By default, assumes input coordinates are in degrees and returns distance in degrees. If degrees=False, then assumes input is in radians and returns distance in radians.

#### sphere_point_along_bearing(lat1,long1,bearing,d,degrees=True):
Returns lat and long of a point angular distance d at azimuth bearing from initial point (lat1,long1).  By default, assumes input coordinates are in degrees and returns distance in degrees. If degrees=False, then assumes input is in radians and returns distance in radians.   

Source: [http://www.movable-type.co.uk/scripts/latlong.html](http://www.movable-type.co.uk/scripts/latlong.html)

#### dir2cart(d)
Converts list or array of vector directions [D,I] to array of cartesian coordinates [x,y,z].

Source: [PMagPy (Lisa Tauxe)](https://github.com/PmagPy/PmagPy)

#### cart2dir(cart):
Converts list of cartesian [x,y,z] values to spherical coordinates [D,I]. 

Source: [PMagPy (Lisa Tauxe)](https://github.com/PmagPy/PmagPy)

#### angle(D1,D2):
Given list of two directions D1,D2 where D1 is for example, [[Dec1,Inc1],[Dec2,Inc2],etc.], returns angle between D1 and D2

Source: [PMagPy (Lisa Tauxe)](https://github.com/PmagPy/PmagPy)

### Look-up Functions

#### get_chron_age(chron,timescale='CK95')
Given chron - a chron name as a string with appended 'y', 'm' or 'o' - will attempt to look up and return the age of the end/young end, middle or beginning/old end of the designated
chron. 

Currently two timescales are available: 'CK95' (Cande & Kent 1995, default) and'GTS12' (Geological Timescale with some astronomically tuned chron boundaries). The relevant .txt files can be found in ./DataFiles - different timescales could also be added here.

#### find_plate_from_name(text)
Returns (in the form of a pandas DataFrame) the numerical and text codes and descriptions of plates with the input string fragment 'text' within their descriptions.

Warning: can return multiple matches! (e.g., Pacific will return any matches with 'Pacific' in the plate description, which is not just the Pacific Plate).

Data source is ./DataFiles/PlateCodes.txt 

#### find_plate_from_number(code)
Given a numerical plate code, will return a pandas DataFrame row with the matching plate code, text code, and description.

Data source is ./DataFiles/PlateCodes.txt 

### Rotation File Functions

## load_rotsets(rotfile)
Loads rotation file rotfile, and returns a list of FiniteRotationSets for each unique set of rotations (which can then be used to build an EulerRotationModel). Assumes that finite rotations for a given plate pair are listed together, in chronological order.

Required header/columns for rotfile: 'MovingPlate','FixedPlate','Chron','EndAge','RotLat','RotLong','RotAng','Kappahat','a','b','c','d','e','f','Points','Segs','Plates','DOF','Source'

- kappahat, a,b,c,d,e,f are parameters from the covariance matrix.
- Points, Segs, Plates and DOF/Degrees of Freedom are parameters for the original Hellinger fit for the specified rotation.

__TODO__ These parameters are not always available: currently they are filled in with dummy values in the input file, which is both inelegant and opaque. Either need some way of dealing with null values, or consider just doing away with them entirely. 
    
## load stage_rots(rotfile)
Takes a rotation file containing stage rotations for a particular plate pair and returns a StageRotationSet object. Note: unlike load_rotsets() expects there to be only plate pair in the file.

__TODO__  check for (and deal with) stage rotations for multiple plate pairs (though probably not commonly required).

Required header/columns for rotfile: 'MovingPlate','FixedPlate','Chron','EndAge','RotLat','RotLong','RotAng','Kappahat','a','b','c','d','e','f','Points','Segs','Plates','DOF','Source'

- kappahat, a,b,c,d,e,f are parameters from the covariance matrix.
- Points, Segs, Plates and DOF/Degrees of Freedom are parameters for the original Hellinger fit for the specified rotation.

#### load_points(pointfile)
Takes a point file and returns a list of Point objects, which can then be used to build PointSet, Path, Platelet objects
    
colums in tab-delimited pointfile: 

'Name','PlateCode','Lat','Lon','FeatureAge','ReconstructionAge'; optionally 'MaxError','MinError','MaxBearing', otherwise 0 by default.
\
# Euler rotation objects

## EulerRotationModel
A class for organising and manipulating multiple sets of finite rotations between different plate pairs. Create by command EulerRotationModel(FiniteRotationSets)

Attributes:
    
- .rotationsets: list of FiniteRotationSets for a particular plate pair 
- .FixedPlate (currently set to 'None' and not currently used for anything)
    
### General/summary functions

#### .summary(self)
Returns a summary DataFrame for FiniteRotationSets within the current rotation model:	plate pair codes and names, number of finite rotations, age range. 

#### .plates(self):
Returns a list of unique plate IDs within the current rotation model.

### Functions for editing rotation model

#### .insert_rots(self,newrots ask=True)
Given a list of FiniteRotationSets (calculated, from load_rotsets(), from another EulerRotationModel, etc.), adds them to the rotation model. 

If ask=True, when a FiniteRotationSet for the same plate pair exists, will ask if you want to replace (default) or skip and keep the original set. If ask=False, replacement is automatic with a notification that it has occurred

#### .remove_rots(self,plate1,plate2=-1)
Remove FiniteRotationSets for the pair plate1, plate2. If plate2 is not defined, it will removed any FiniteRotationSet that involves plate1.

#### .newagemodel(self, timescale='CK95')
WARNING: experimental. **Should** overwrite old rotationsets with new ones where rotation ages timed to reversal ages have been recalibrated to ages in the specified timescale. Current options: 'CK95', 'GTS12'. 

__TODO__ Have some way of defining the current agemodel, possibly by defining it on loading and assigning anomaly ages within rotationsets at that stage. 

### Functions for searching within the rotation model

#### .find_pairs(self,plate)
Searches for any FiniteRotationSets in which the specified plate is one of the plate pair; returns the pairs of platecodes for those FiniteRotationSets. 

(used internally by .find_circuits(), which is probably more immediately useful).

#### .find_circuit(self,startplate,endplate,preferred=-1)
Searches for ways to link startplate to endplate using rotations within the currently defined rotation model. The returned circuit is a list of platecodes that starts with	startplate and ends with endplate. If multiple possible circuits exist, they will all be returned;	can filter by setting a preferred platecode, although this may not always restrict the list to just one! 

NB: This works by building out a list of all possible circuits from the startplate, and then selecting the ones that end with the endplate. 

#### .get_rots(self,movingplate,fixedplate,ages=[],preferred=-1)
Returns a FiniteRotationSet for the movingplate-fixedplate pair. If a list of ages is supplied, 	then the rotationset will consist of the interpolated finite rotations for those ages; otherwise it will be all	available ages in the rotationset. If more than one circuit is possible, then the circuit can be forced through preferred.

If multiple movingplate-fixedplate rotationsets already exist within the rotation model (hopefully unlikely), or no viable plate circuit can be found, an empty list is returned. If multiple viable plate circuits are found, it will use the first one in the list returned from find_circuit().

__TODO__ Interactively select from multiple plate circuits (or existing rotationsets?)

### Functions that generate other useful things using rotation model

#### .synthetic_APWP(self,moving_plate,absolute_ref_frame,ages,SetName='APWP',PlotColor='orange',PlotLevel=5)
Returns a Flowline object that predicts the Apparent Polar Wander path that should have been generated by motion of moving_plate in absolute_ref_frame (i.e reconstructed position of geographic North Pole in the moving_plate reference frame) for specified list of age points. Note that currently, there is no restriction on what reference frame is used, but this will only be meaningful if it is an absolute frame (e.g, hotspot frames 001/Atlantic or 003/Pacific)

__TODO__ restrict possible absolute_reference_frames to ones that produce meaningful results?

#### .hotspot_track(self,point,absolute_ref_frame,ages,SetName='Hotspot Track',PlotColor='orange',PlotLevel=5)
Returns a Flowline object that predicts the track made by a hotspot currently located at the coordinates of the input Point object at the supplied list of ages if it is fixed to absolute_ref_frame. Note that currently, there is no restriction on what reference frame is used, but this will only be meaningful if it is an absolute frame (e.g, hotspot frames 001/Atlantic or 003/Pacific)

__TODO__ restrict possible absolute_reference_frames to ones that produce meaningful results?

## FiniteRotationSet

## StageRotationSet

## EulerRotation

# Objects acted upon by Euler rotations

## Point


## PointSet

## Boundary 
Derived from PointSet but explicitly defined as a line. 

## Platelet
Derived from PointSet but explicitly defined as a closed polygon.

## Flowline

## APWP

Also currently AMS_Locality and PMag_Locality classes, but these may be too specialised; likely to be put into a separate file eventually.