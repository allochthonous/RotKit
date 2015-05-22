import numpy as np
import pandas as pd
import sys,os
import mpl_toolkits.basemap.pyproj as pyproj
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from triangle import triangulate, plot as tplot
#imports fortran complied point rotation and ellipse routines - MUST BE COMPILED using fortran2py on your host system.
sys.path.append('/Users/crowan/PlateRotPy/')
import rotatetest as rot
#Imports a table for looking up plate codes - potentially useful when building new rotations?
os.chdir('/Users/crowan/Dropbox/Research/Platerots/')
platecodes=pd.read_table('PlateCodes.txt', header=None, names=['NumCode','LetterCode','Description'],index_col='NumCode')

def dir2cart(d):
   # converts list or array of vector directions, in degrees, to array of cartesian coordinates, in x,y,z
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
    cart= np.array([ints*np.cos(decs)*np.cos(incs),ints*np.sin(decs)*np.cos(incs),ints*np.sin(incs)]).transpose()
    return cart

def cart2dir(cart):
    """
    converts a direction to cartesian coordinates
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

def find_plate_from_name(text):
    #given a string fragment (e.g., Pacific will return any matches from the plate description. Warning: can return multiple matches!
    return platecodes[platecodes.Description.str.contains(text)]

def find_plate_from_number(code):
    return platecodes.loc[code]

def load_rots(file):
    #load rotation poles into a pandas table. Eventually may need a file type flag if loading from different sources. And maybe a seperate one for stage poles?
    rot_data=pd.read_table(file, header=None, names=['RotName','Method','EndAge','RotLat','RotLong','RotAng','Kappahat','a','b','c','d','e','f','Points','Segs','Plates','DOF','Source'])
    start_ages=pd.Series(np.repeat(0,len(rot_data)), name='StartAge')
    #currently I'm using named plates - if I move to plate codes (as I should), could do a lookup here to add in plate names as well.
    plates=pd.DataFrame(rot_data.RotName.str.split('-').tolist(), columns=['MovingPlate','FixedPlate'])
    #combine covariances into single entity because they are always used together - may not be the most elegant way to do it...
    covs=[]
    for a,b,c,d,e,f in zip(rot_data.a,rot_data.b, rot_data.c,rot_data.d,rot_data.e, rot_data.f): covs.append([a,b,c,d,e,f])
    covs=pd.Series(covs,name='Covariances')
    rot_data=pd.concat([rot_data.RotLat,rot_data.RotLong,rot_data.RotAng,rot_data.Kappahat,covs,start_ages,rot_data.EndAge,plates,rot_data.Points,rot_data.DOF,rot_data.Source], axis=1)
    return rot_data
    
def load_points(file):
    #again - at some stage could well need to deal with multiple formats: this is my preferred one.
    point_data=pd.read_table(file, header=None, names=['PointLat','PointLong','Plate','ReconstructionAge', 'FeatureAge', 'ErrEllipseMax','ErrEllipseMin','ErrEllipseMaxBearing'])
    #if no error ellipse data, will fill it with 0s rather than NaNs. May need to tie this to specific columns?
    point_data=point_data.fillna(0)
    return point_data

def make_points(lats,lons,ages,platecode='None'):
    #turn lists of points and ages into pandas dataframe:  by default assumes they are current position
    rotages=np.repeat(0,len(lats))
    plates=np.repeat(platecode,len(lats))
    errmax=np.repeat(0,len(lats))
    errmin,erraz=errmax,errmax
    return pd.DataFrame(np.column_stack((lats,lons,plates,rotages,ages,errmax,errmin,erraz)),
                            columns=['PointLat','PointLong','Plate','ReconstructionAge', 'FeatureAge', 'ErrEllipseMax','ErrEllipseMin','ErrEllipseMaxBearing'])

def rotate_points(points,rots, order=1):
    if order<1 or order>2:
        #normally it wouldn't matter but fortran will complain
        print 'Order needs to be 1 or 2.\n'
        rotated=[]
    else:
        #rotates the points by the specified Euler poles. Both input as pandas dataframes.
        #currently sends to the wrapped fortran functions - so must put path to specifically compliled fortran2py .so file 
        #order=1; rotate all the points by each pole in sequence. order=2: rotate each point in sequence by every pole. 
        inpoints=np.column_stack((points.PointLat,points.PointLong))
        # need to extract covariances back into separate rows to feed to fortran
        covs=[[],[],[],[],[],[]]
        for index,row in rots.iterrows():
            for cov in np.arange(0,6):
                covs[cov].append(row.Covariances[cov])
        inpoles=np.column_stack((rots['RotLat'],rots['RotLong'],rots['RotAng'],rots['Kappahat'],covs[0],covs[1],covs[2],covs[3],covs[4],covs[5],rots['EndAge']))     
        result=rot.rotatepts(inpoints,inpoles,order)
        if order==1:
            plates=np.tile(points.Plate,len(rots))
            featages=np.tile(points.FeatureAge,len(rots))
            rotages=np.repeat(rots.EndAge.values,len(points)) #strange how tile and repeat require different forms...
        else:
            #if order=2 implied
            plates=np.repeat(points.Plate.values,len(rots))
            featages=np.repeat(points.FeatureAge.values,len(rots))
            rotages=np.tile(rots.EndAge,len(points))   
        #rotated=np.column_stack((result[:,0],result[:,1],plates,rotages,featages,result[:,3],result[:,4],result[:,5]))
        rotated=pd.DataFrame(np.column_stack((result[:,0],result[:,1],plates,rotages,featages,result[:,3],result[:,4],result[:,5])),
                            columns=['PointLat','PointLong','Plate','ReconstructionAge', 'FeatureAge', 'ErrEllipseMax','ErrEllipseMin','ErrEllipseMaxBearing'])
    return rotated   

def fetch_rotation(panda_row):
    # create [lat,lon,angle] of rotation from pole dataframe row or series
    rotation=[panda_row.iloc[0]['RotLat'],panda_row.iloc[0]['RotLong'],panda_row.iloc[0]['RotAng']]
    return rotation

def make_rotmat(rotation):
    # works on array [lat,long,angle]
    lat,lon,ang=rotation[0]*np.pi/180,rotation[1]*np.pi/180,rotation[2]*np.pi/180,
    num=1-np.cos(ang)
    px=np.cos(lat)*np.cos(lon)
    py=np.cos(lat)*np.sin(lon)
    pz=np.sin(lat)
    b=np.matrix([[(px*px*num )+ np.cos(ang),(px*py*num )-(pz*np.sin(ang)),(px*pz*num)+(py*np.sin(ang))],
                 [(py*px*num)+(pz*np.sin(ang)),(py*py*num)+np.cos(ang),(py*pz*num)-(px*np.sin(ang))],
                 [(px*pz*num)-(py*np.sin(ang)),(pz*py*num)+(px*np.sin(ang)),(pz*pz*num)+np.cos(ang)]])
    return b

def rotmat_to_pole(rot_matrix):
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

def make_covmat(covs):
    a=np.matrix([[covs[0],covs[1],covs[2]],[covs[1],covs[3],covs[4]],[covs[2],covs[4],covs[5]]])
    return a 

def add_rots(rot1,rot2):
    #where rot1,rot2 are arrays of [lat,long,angle]
    #adapted from fortran code provided by David Rowley, based on vector addition as described in Cox and Hart  p. 232 and 233
    #1. Check if rotation angles are 0    
    if rot1[2]==0 and rot2[2]==0: newpole=[0,0,0]
    #2. Check if poles are identical (sum to zero)
    elif -rot1[2]==rot2[2] and rot1[0]==rot2[0] and rot1[1]==rot2[1]: newpole=[rot1[0],rot1[1],0]
    elif -rot1[2]==rot2[2] and rot1[0]==-rot2[0] and rot1[1]==rot2[1]+180: newpole=[rot1[0],rot1[1],0]
    elif -rot1[2]==rot2[2] and rot1[0]==-rot2[0] and rot1[1]==rot2[1]-180: newpole=[rot1[0],rot1[1],0]
    else:
        rot_matrix=make_rotmat(rot2)*make_rotmat(rot1)
        newpole=rotmat_to_pole(rot_matrix)
    return newpole

def invert_covrots(rotpars, inverting='plate'):
    # calculate inverse rotation and its covariance matrix - adapted from invrot.f (original November 3 1988 J-Y R)
    # a=b**t; cova=b*covb*(b**t) where b, covb =rotation and covariance matrices of original rotation, and a is the transpose of b
    # Input: pandas dataframe (will cycle through multiple rows and invert them)
    # Output: pandas dataframe
    # By default, switches the fixed and moving plate. If inverting set to 'time' will switch
    inverted=pd.DataFrame()
    for index, row in rotpars.iterrows():
        rotation=[row.RotLat,row.RotLong,row.RotAng]
        rot_matrix=make_rotmat(rotation)
        cov_matrix=make_covmat(row.Covariances)
        t_rot_matrix=np.transpose(rot_matrix)
        inv_cov=rot_matrix*cov_matrix*t_rot_matrix
        covariances=[inv_cov[0,0],inv_cov[0,1],inv_cov[0,2],inv_cov[1,1],inv_cov[1,2],inv_cov[2,2]]
        rotation[2]=-rotation[2]
        if inverting=='time': 
            startage,endage=row.EndAge,row.StartAge
            fixedplate,movingplate=row.FixedPlate,row.MovingPlate
        else: 
            fixedplate,movingplate=row.MovingPlate, row.FixedPlate
            startage,endage=row.StartAge,row.EndAge    
        newrow=pd.DataFrame([[rotation[0],rotation[1],rotation[2],rotpars.iloc[0].Kappahat,covariances,startage,endage,fixedplate,movingplate,rotpars.iloc[0].Points,rotpars.iloc[0].DOF,rotpars.iloc[0].Source]],
                                columns=['RotLat','RotLong','RotAng','Kappahat','Covariances','StartAge','EndAge','MovingPlate','FixedPlate','Points','DOF','Source'])
        inverted=pd.concat([inverted,newrow], ignore_index=True)
    return inverted

def add_covrots(rotationstoadd):
    # Add poles and their covariances. 
    # If more than two sequentially adds them until final product is moving plate of 1st rotation wrt fixed plate of final rotation; endage of first to startage of last
    # input: pandas dataframe: output Dataframe of all sums; final row is final rotation 
    # Chang 1990 - covariance of combined rotations A and B C=AB, with covariance cova and covb: covc ~= B(t)*covA*B+covB
    # a linear approximation that assumes rotations *associated with the errors* are small enough to be treated as infintesimal so can be added as vectors.
    # Condition for this that det(cova) and det (covb) <<1, which seems to be more than met with typical rotation pole covariances.
    # Adapted from addrot.f (most notes below copied from there) ======  November 3 1988 J-Y R ======= MODIFIED OCTOBER 14, 1991
    #     chat=bhat.ahat 
    #         =b*phi(tbhat)*a*phi(tahat)
    #         =b*a*(a**t)*phi(tbhat)*a*phi(tahat)
    #         =b*a*phi((a**t)*tbhat)*phi(tahat)
    #         ~b*a*phi((a**t)*tbhat +tahat)  up to linear approximations
    #
    #     tchat~(a**t)*tbhat +tahat
    #     cov(tchat)~(a**t)*cov(tbhat)*a + cov(tahat)
    #
    # IF THE DEGREES OF FREEDOM IN KAPPAHAT FOR A AND FOR B ARE BOTH
    # GREATER THAN 50, THEY ARE EACH ASSUMED TO THE THE TRUE KAPPAS.
    # OTHERWISE A POOLED KAPPAHAT IS USED.  IN THIS CASE, AN ERROR
    # MESSAGE IS PRODUCED IF AN F-TEST WOULD CONTRADICT THE EQUALITY
    # OF THE KAPPAS.
    # 
    # IF KAPPA IS KNOWN FOR ONE ROTATION, BUT ESTIMATED FOR THE OTHER,
    # AND IF THE ESTIMATED KAPPA HAS FEWER THAN 50 DEGREES OF FREEDOM,
    # THE ESTIMATED KAPPA IS ASSUMED TO BE CORRECT, AND AN ERROR MESSAGE
    # IS PRINTED.
    # 
    added=pd.DataFrame()
    rot1= rotationstoadd.ix[0]
    endage=rot1.EndAge
    movingplate=rot1.MovingPlate
    for index,rot2 in rotationstoadd[1::].iterrows():
        rotation1=[rot1.RotLat,rot1.RotLong,rot1.RotAng]
        rotation2=[rot2.RotLat,rot2.RotLong,rot2.RotAng]
        newrot=add_rots(rotation1,rotation2)
        rot_matrix_b=make_rotmat(rotation2)
        cov_matrix_a=make_covmat(rot1.Covariances)
        cov_matrix_b=make_covmat(rot2.Covariances)
        
        dof_a, dof_b = rot1.DOF, rot2.DOF
        #calculate combined kappa
        ndata=rot1.Points+rot2.Points
        
        if (dof_a>=50 and dof_b>=50):
            cov_matrix_a=cov_matrix_a/rot1.Kappahat
            cov_matrix_b=cov_matrix_b/rot2.Kappahat
            dof_c = 10000.
            Kappahat_c = 1.
            plev = -1.
            fkap1, fkap2 = 0.,0.
        else:
            # I am cheating here. if one or both kappahats < 50 then should be calculating a pooled kappahat.Should be implemented at some point 
            cov_matrix_a=cov_matrix_a/rot1.Kappahat
            cov_matrix_b=cov_matrix_b/rot2.Kappahat
            dof_c = 10000.
            Kappahat_c = 1.
            plev = -1.
            fkap1, fkap2= 0.,0.
    
        cov_matrix_c=(np.transpose(rot_matrix_b)*cov_matrix_b*rot_matrix_b)+cov_matrix_a
        covariances=[cov_matrix_c[0,0],cov_matrix_c[0,1],cov_matrix_c[0,2],cov_matrix_c[1,1],cov_matrix_c[1,2],cov_matrix_c[2,2]]
        #Correctly gets plates and ages right if adding two poles with either 1) same start age and end age or 2) same moving plate and fixed plate
        startage=rot2.StartAge
        fixedplate=rot2.FixedPlate
        newrow=pd.DataFrame([[newrot[0],newrot[1],newrot[2],Kappahat_c,covariances,startage,endage,movingplate,fixedplate,ndata,dof_c,'Addition']],
                                columns=['RotLat','RotLong','RotAng','Kappahat','Covariances','StartAge','EndAge','MovingPlate','FixedPlate','Points','DOF','Source'])
        added=pd.concat([added,newrow], ignore_index=True)
        rot1=newrow.ix[0] #so next pole is added to the previous sum
    return added


def get_stage_rot(rot1,rot2, return_pseud=0):
    # stage rotation between two poles where rot1, rot2 are arrays of [lat,long,angle]
    A1=make_rotmat(rot1)
    A2=make_rotmat(rot2)
    S=A2*np.transpose(A1)
    stagerot=rotmat_to_pole(S)
    if return_pseud==0:
        return stagerot
    else: 
        stagelat,stagelon,stageang=stagerot[0]*np.pi/180,stagerot[1]*np.pi/180,stagerot[2]*np.pi/180
        pseud=[stageang*np.cos(stagelat)*np.cos(stagelon),stageang*np.cos(stagelat)*np.sin(stagelon)]
        return stagerot,pseud      

def get_stage_covrots(rotations):
    #from a list of poles, get the stage poles between them
    #currently a bit inflexibile on how it gets startage and endage - may need to tweak a little bit for when not using reconstruction poles
    rots1=np.column_stack((rotations.RotLat[:-1],rotations.RotLong[:-1],rotations.RotAng[:-1]))
    rots2=np.column_stack((rotations.RotLat[1:],rotations.RotLong[1:],rotations.RotAng[1:]))
    covs1=rotations.Covariances[:-1]
    covs2=rotations.Covariances[1:]
    stagerots=np.array(())
    covariances=[]
    for rot1,rot2,cov1,cov2 in zip(rots1,rots2,covs1,covs2): 
        stagerots=np.append(stagerots,get_stage_rot(rot1,rot2))
        #calculation of covariance matrix for stage rotation between A1 and A2 = A1*(cov2-cov1)
        #from Appendix to Doubrovine & Tarduno 2008
        cov_matrix_s=make_rotmat(rot1)*(make_covmat(cov2)-make_covmat(cov1))
        covariances.append([cov_matrix_s[0,0],cov_matrix_s[0,1],cov_matrix_s[0,2],cov_matrix_s[1,1],cov_matrix_s[1,2],cov_matrix_s[2,2]])
    stagerots=np.reshape(stagerots,(len(rots1),3))
    covariances=pd.Series(covariances,name='Covariances')
    #not sure that this is totally valid - may be better to just put in smallest number of points?
    points=(rotations.Points.values[:1]+rotations.Points.values[:-1])/2
    DOF=np.repeat(10000,len(stagerots))
    sources=np.repeat('Interpolated Stage Pole',len(stagerots))
    output=pd.DataFrame(np.column_stack((stagerots[:,0],stagerots[:,1],stagerots[:,2],rotations.Kappahat[:-1],covariances,
                    rotations.EndAge[:-1],rotations.EndAge[1:],rotations.MovingPlate[:-1],rotations.FixedPlate[:-1], points,DOF,sources)),
                    columns=['RotLat','RotLong','RotAng','Kappahat','Covariances','StartAge','EndAge','MovingPlate','FixedPlate','Points','DOF','Source'])
    return output

def interpolate_rots(rot1,rot2,age, return_pseud=0):
    # interpolate between two rotations where rot1, rot2 are arrays of [lat,long,angle,age]
    # Adapted from code provided by Pavel Doubrovine - ibfrlib0.3f
    # make sure age is a float - if not, xi will be calculated as 0
    # 1. Check ages are ordered correctly
    if rot2[3]<rot1[3]:
        print 'age of second pole must be greater than age of first pole.'
        intpole=[]
    # 2. Check if time step is actually within ages of two poles
    elif age<rot1[3] or age>rot2[3]: 
        print 'requested timestep does not lie between pole ages.'
        intpole=[]
    # 3. Check if poles are the same     
    elif rot1[0:3]==rot2[0:3]: intpole=rot1[0:3].append(age)
    else:
        #xi is fraction of interval between two poles
        xi=(age-rot1[3])/(rot2[3]-rot1[3])
        A1=make_rotmat(rot1)
        #calculate the stage rotation and calculate rotation matrix for fractional rotation about it
        stage_rot, pseud=get_stage_rot(rot1,rot2,1)
        stage_rot[2]=stage_rot[2]*xi
        S_rotmat=make_rotmat(stage_rot)
        #calculate interpolated rotation
        A=S_rotmat*A1
        intpole=rotmat_to_pole(A)
    if return_pseud==0: return intpole
    else: return intpole,pseud
          
def interpolate_covrots(rotations,ages):
    # Produce interpolated poles and covariances, following the method of Doubrovine and Tarduno (2008)
    # Adapted from code provided by Pavel Doubrovine - ibfrlib0.3f 
    # Input: rotations is a pandas DataFrame, ages is a list of target interpolation ages
    # If there is a list of rotations, will search for the two bracketing rotations for each target age
    # Output: pandas dataframe
    ages.sort() #just in case
    for age in ages:
        print 'Age: '+`age`+'  '
        age=float(age) #also just in case
        if age<rotations.EndAge.min() or age>rotations.EndAge.max(): print `age`+' is outside age range of given poles.\n'
        #checks there is a point in interpolating. If a rotation exists with the given age, can just use that.
        elif len(rotations[rotations.EndAge==age])>0: 
            print `age`+' already exists.\n'
            interpolated=rotations[rotations.EndAge==age][:1].values
        else: 
            A1=rotations[rotations.EndAge<age][-1:]
            A2=rotations[rotations.EndAge>age][:1]
            A1rot=[A1.iloc[0].RotLat,A1.iloc[0].RotLong,A1.iloc[0].RotAng,A1.iloc[0].EndAge]
            A2rot=[A2.iloc[0].RotLat,A2.iloc[0].RotLong,A2.iloc[0].RotAng,A2.iloc[0].EndAge]
            newrot,pseud=interpolate_rots(A1rot,A2rot,age,1)
            #if Euler pole location and covariance matrices the same, keep them the same.
            if A1rot[0:2]==A2rot[0:2] and A1.Covariances.values==A2.Covariances.values: newcov=A1.Covariances.values
            #otherwise its complicated calculation time!
            else:
                xi=(age-A1rot[3])/(A2.rot[3]-A1rot[3])
                cov_matrix_A1=make_covmat(A1.Covariances)
                cov_matrix_A2=make_covmat(A2.Covariances)             
                pseud1=pseud*xi
         
        print interpolated
        print '\n'

def get_plate_velocity(point,rot,duration=1):
    #calculates the magnitude and direction of the plate motion vector from a point [latitude,longitude] and a Euler rotation [latitude, longitude,rate]
    #if duration given will convert to a mm/yr or km/Myr rate
    #velocity of point on Earth's surface=cross product of Euler vector (omega) and position vector (r)
    #N-S component of v vNS = a*|rotrate|*cos(polelat)*sin(sitelong-polelong)
    #E-W component of v vEW =  a*|rotrate|* [cos(sitelat)*sin(polelat)-sin(sitelat)*cos(polelat)*cos(sitelong-polelong)]
    #where a=Earth radius
    #rate of motion = sqrt (vNS^2+vEW^2)
    #azimuth = 90-atan[vNS/vEW]
    Earthrad=6371
    pointlat=point[0]*np.pi/180
    pointlong=point[1]*np.pi/180
    rotlat=rot[0]*np.pi/180
    rotlong=rot[1]*np.pi/180
    rotrate=abs(rot[2]*np.pi/180)
    vNS=Earthrad*rotrate*np.cos(rotlat)*np.sin(pointlong-rotlong)
    vEW=Earthrad*rotrate*(np.cos(pointlat)*np.sin(rotlat)-np.sin(pointlat)*np.cos(rotlat)*np.cos(pointlong-rotlong))
    rate=np.sqrt(vNS**2+vEW**2)/duration
    az=90-(180/np.pi*np.arctan(vNS/vEW))
    return [rate,az,vNS,vEW]

def get_velsquared(unit,m,rotation):
    #input: mesh nodes in lat,long array + map projection
    tlon,tlat=m(unit[:,0],unit[:,1])
    thing=Polygon(np.column_stack((tlon,tlat)))
    area=thing.area/1e6
    x,y=thing.centroid.xy
    x,y=m(x,y,inverse=True)
    rate=get_plate_velocity(([y[0],x[0]]),rotation)[0]
    velsquared=area*rate**2
    print rate, area
    return area,velsquared

def RMS_calc(lons,lats,rotation):
    #Define an equal area projection around the plate
    m = pyproj.Proj("+proj=aea +lat_1="+`min(lats)`+" +lat_2="+`max(lats)`+" +lat_0="+`(min(lats)+max(lats))/2`+" +lon_0="+`(min(lons)+max(lons))/2`)
    #Create irregular triangular mesh for the plate polygon
    data={}
    data['vertices']=np.column_stack((lons,lats))
    segs=[[0,len(lons)-1]]
    for i in np.arange(0,len(lons)-1): segs.append([i,i+1])
    data['segments']=np.array(segs)
    t = triangulate(data, 'pa50q30') #starting off with too small an area can cause crashes
    #can refine further using r switch - note I do so right now but the difference in the overall result is small
    t2= triangulate(t, 'ra10q30')
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    tplot.plot(ax1, **t)
    plt.title('1st Triangulation - min area 50, min angle 30')
    tplot.plot(ax2, **t2)
    plt.title('Refined mesh - min area 10, min angle 30')
    plt.tight_layout()
    units=t2['vertices'][t2['triangles']]
    vel_lat,vel_lon,vel_mag,vel_az= [],[],[],[]
    totvelsquared=0.
    totarea=0.
    for unit in units:
        tlon,tlat=m(unit[:,0],unit[:,1])
        thing=Polygon(np.column_stack((tlon,tlat)))
        area=thing.area/1e6
        x,y=thing.centroid.xy
        x,y=m(x,y,inverse=True)
        platevel=get_plate_velocity(([y[0],x[0]]),rotation)
        totvelsquared=totvelsquared+(area*platevel[0]**2)
        totarea=totarea+area
        vel_lon.append(x)
        vel_lat.append(y)
        vel_mag.append(platevel[0])
        vel_az.append(platevel[1])
    rms=np.sqrt(totvelsquared/totarea)
    velgrid=pd.DataFrame(np.column_stack((vel_lon,vel_lat,vel_mag,vel_az)), columns=['Lons','Lats','Rate','Azimuth'])
    return rms, totarea, velgrid

def interpolate_plate_boundary(lons,lats,steplength=100):
    # make higher resolution plate boundaries by interpolating along great circle lines. Can be a segment or close polygon.
    # something that is done automatically in map projections but an issue for long line segments when creating plate meshes, for example.
    # by default, looks to interpolate to a roughly 100 km step length.    
    g = pyproj.Geod(ellps='WGS84')
    temp=np.column_stack((lons[:-1],lats[:-1],lons[1:],lats[1:]))
    #calculates distance and azimuths between adjacent points
    f_az,b_az,dist=g.inv(temp[:,0],temp[:,1],temp[:,2],temp[:,3])
    #for each point, interpolate along its bearing to the next at the specified spacing
    int_lon,int_lat=[],[]
    for slon,slat,az,d in zip(temp[:,0],temp[:,1],f_az,dist):
        int_lon.append(slon)
        int_lat.append(slat)
        step=round(d/(steplength*1000))
        for i in np.arange(1,step):
            nlon,nlat,naz=g.fwd(slon,slat,az,i*d/step)
            int_lon.append(nlon)
            int_lat.append(nlat)
    int_lon.append(lons[-1])        
    int_lat.append(lats[-1])
    int_lat=np.array(int_lat)
    int_lon=np.array(int_lon)
    return int_lon, int_lat
 
def reduce_plate_boundary(lons,lats,threshold_length=100, threshold_angle=5):
    #remove excess nodes which do not actually record a sharp change in direction.
    #becuase the triangulation seems to run into memory problems if go above ~400 nodes (or, potentially, if nodes too closely spaced)
    g = pyproj.Geod(ellps='WGS84')
    temp=np.column_stack((lons[:-1],lats[:-1],lons[1:],lats[1:])) 
    f_az,b_az,dist=g.inv(temp[:,0],temp[:,1],temp[:,2],temp[:,3])
    dist=dist/1000.
    red_lon=[lons[0]]
    red_lat=[lats[0]]
    summed_length=0.
    for lon,lat,length,angle in zip(temp[:,0],temp[:,1],dist,abs(f_az[:-1]-f_az[1:])):
        summed_length=summed_length+length
        #print summed_length, angle
        if summed_length>threshold_length and angle>threshold_angle:
            red_lon.append(lon)
            red_lat.append(lat)
            summed_length=0
    return red_lon, red_lat        
  
   
