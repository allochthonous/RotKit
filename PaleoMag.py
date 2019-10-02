from pmag import fshdev, dir2cart, cart2dir, dotilt, fisher_mean, angle, vfunc, dokent, Tmatrix, dimap,circ
import numpy as np
import pandas as pd
from EulerRots import sphere_ang_dist, sphere_bearing, Point, PointSet
import f_untilt as f_untilt

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.lines as mlines


def load_poles(polefile, use_A95=False):
    """
    Takes a pole file and returns a list of Point objects based on the pole parameters, including dm/dp
    
    colums in tab-delimited file: Name,PlateCode,PoleLat,PoleLon,FeatureAge,ReconstructionAge;
    dm,dp for error ellipse, or set use_a95=True to look for column a95.
    """
    pole_data=pd.read_csv(polefile, sep='\t')
    if use_A95==True: pole_data['MaxError'], pole_data['MinError']=pole_data['a95'], pole_data['a95'] 
    else: pole_data['MaxError'], pole_data['MinError']=pole_data['dm'], pole_data['dp']
    pole_data['MaxBearing']=0.
    return [PMagMeanPole(row) for i,row in pole_data.iterrows()]

class PMagMeanPole(Point):
    """
    Represents a paleomagnetic mean pole; behaves as point does, with some additional attributes and 
    functions
    """
    def __init__(self, PointPars,PlotColor='grey',PlotLevel=5,PlotSymbolSize=12):
        Point.__init__(self,PointPars,PlotColor,PlotLevel,PlotSymbolSize)
        self.PMagPars=pd.Series([PointPars.a95,PointPars.k,PointPars.n], 
                                index=['a95','k','n'])
            
class PMagMeanPoleSet(PointSet):
    """
    A collection of paleomagnetic site mean poles
    """
    def __init__(self,PointList,age,SetName='PMagPoleSet',PlotColor=None,PlotLevel=None):
        # assumes all the poles have a common age/fall within common age window 
        PointSet.__init__(self,PointList,SetName,PlotColor,PlotLevel)
        self.FeatureAge=age
        self.PlateCode=PointList[0].PlateCode
    
    def fisher_mean(self, raw=False):
        """
        calculates the Fisher mean of the poles in this set and by default returns it as a Point object; 
        if instead raw is set to True it returns a Pandas Series  
        """
        result=fisher_mean(np.column_stack((self.summary().Lon, self.summary().Lat)))
        if raw==True:
            return pd.Series([result['dec'],result['inc'],result['alpha95'],result['k'],result['n']],
                            index=['Dec','Inc','a95','k','n'])
        else:
            return PMagMeanPole(pd.Series([self.SetName+' Mean Pole', self.PlateCode, result['inc'],result['dec'],
                                           self.ReconstructionAge, self.FeatureAge, result['alpha95'],result['alpha95'],0.,
                                           result['alpha95'],result['k'],result['n']],
                               index=['Name','PlateCode','Lat','Lon','FeatureAge','ReconstructionAge',
                                      'MaxError','MinError','MaxBearing','a95','k','n']), PlotLevel=10) 
    
    def bootstrap_mean(self, trials=1000, raw=False):
        """
        Performs parametric bootstrap to calculate mean direction
        """
        decs,incs=[],[]
        for i in range(trials):
            newmean=fisher_mean([pole.resample(1)[0] for pole in self.points])
            decs.append(newmean['dec'])
            incs.append(newmean['inc'])
        kmean=dokent(zip(decs,incs),1)
        if raw==True:
            return pd.Series([kmean['dec'],kmean['inc'],kmean['Zeta'],kmean['Eta'],
                              sphere_bearing(kmean['inc'],kmean['dec'],kmean['Zinc'],kmean['Zdec'])],
                              index=['Dec','Inc','MaxError','MinError','MaxBearing'])               
        else: 
        #because it is not a Fisher mean, no a95, k... might need to think about how to deal with that...
            return Point(pd.Series([self.SetName+' Mean Pole', self.PlateCode, kmean['inc'],kmean['dec'],
                                           self.ReconstructionAge, self.FeatureAge, kmean['Zeta'],kmean['Eta'],
                                            sphere_bearing(kmean['inc'],kmean['dec'],kmean['Zinc'],kmean['Zdec'])],
                               index=['Name','PlateCode','Lat','Lon','FeatureAge','ReconstructionAge',
                                      'MaxError','MinError','MaxBearing']), PlotLevel=10) 
    

   
class Synthetic_PMagSite(object):
    """
    A paleomagnetic site 'mean direction' that does not have an underlying set of
    directions associated with it
    """
    def __init__(self,direction,tiltinfo=[],name='SyntheticDir'):
        """Return object
        direction should be a DataFrame with columns age,D,I (in geographic coordinates),n,a95,and optionally k
        tiltinfo should be a DataFrame with columns BedStrike,BedDip,FoldMinAge,FoldMaxAge
        """
        self.Name=name
        self.Dec=direction.D
        self.Inc=direction.I
        self.a95=direction.a95
        self.n=int(direction.n)
        if 'k' in direction: self.k=direction.k
        else: self.k=((140/self.a95)**2)/self.n
        self.R=self.n-(self.n-1)/self.k
        if 'age' in direction: self.SiteAge=direction.age
        else: self.SiteAge=0
        self.BedStrike=tiltinfo.BedStrike
        self.BedDip=tiltinfo.BedDip
        self.FoldMinAge=tiltinfo.FoldMinAge
        self.FoldMaxAge=tiltinfo.FoldMaxAge
        
    def tilt_info(self):
        return pd.Series([self.BedStrike,self.BedDip,self.FoldMinAge,self.FoldMaxAge],index=['BedStrike','BedDip','FoldMinAge','FoldMaxAge'])
    
    def resample(self, rotate=True):
        """
        Parametric sample of fisher distribution for the direction.
        Setting rotate to 'False' will produce a fisher
        distribution with the site k, centred on 0,90 (useful for implementation 
        of MM90 fold test, unlikely to be needed otherwise)
        """
        tempD,tempI=[],[]
        rad=np.pi/180. # converts from degrees to radians
        #rotmat=np.matrix([dir2cart([self.Dec,90-self.Inc,1.]),dir2cart([self.Dec-90,0.,1.]),dir2cart([self.Dec,self.Inc,1.])])
        for i in range(self.n):
            dec,inc=fshdev(self.k)
            if rotate==False:
                tempD.append(dec)
                tempI.append(inc)  
            else: 
                #newdir=cart2dir(dir2cart([dec,inc,1])*rotmat)
                X=dir2cart([dec,inc,1.]) # get cartesian coordinates of dec,inc
                # get some sines and cosines of new coordinate system
                sa,ca= -np.sin((self.Dec)*rad),np.cos((self.Dec)*rad)
                cdp,sdp= np.cos((90-self.Inc)*rad),np.sin((90-self.Inc)*rad)
                # do the rotation
                xc=X[0]*(sa*sa+ca*ca*cdp)+X[1]*(ca*sa*(1.-cdp))+X[2]*sdp*ca
                yc=X[0]*ca*sa*(1.-cdp)+X[1]*(ca*ca+sa*sa*cdp)-X[2]*sa*sdp
                zc=X[0]*ca*sdp-X[1]*sdp*sa-X[2]*cdp
                # convert back to direction:
                newdir=cart2dir([xc,yc,-zc])
                tempD.append(newdir[0])
                tempI.append(newdir[1])
        return np.column_stack((tempD,tempI))
                                        
    def mean_resample(self):
        """
        Returns the mean direction of parametric sampling (i.e. should fall within the
        95% error ellipse of the original direction 95% of the time.
        """                
        return fisher_mean(self.resample())
        
    def fast_resample(self, rotate=True):
        """
        Parametric sample of fisher distribution for the direction, with 
        a call to a fortran function wrapped with f2py. Much faster
        when taking multiple draws for a bootstrap simulation. 
        Setting rotate to 'False' will produce a fisher
        distribution with the site k, centred on 0,90 (useful for implementation 
        of MM90 fold test, unlikely to be needed otherwise)
        """
        if rotate==True:
            return f_untilt.get_fish(self.Dec,self.Inc,self.k,self.n)
        else:
            return f_untilt.get_fish(0.,90.,self.k,self.n)
        
    def untilt(self,untilting=1):
        """
        performs a partial tilt correction according to untilting (0=geographic, 1=stratigraphic)
        adapted from the dotilt routine in PMagPy (Lisa Tauxe) 
        """
        rad=np.pi/180. # converts from degrees to radians
        X=dir2cart([self.Dec,self.Inc,1.]) # get cartesian coordinates of dec,inc
        # get some sines and cosines of new coordinate system
        partial_dip=self.BedDip*untilting
        sa,ca= -np.sin((self.BedStrike+90.)*rad),np.cos((self.BedStrike+90.)*rad)
        cdp,sdp= np.cos(partial_dip*rad),np.sin(partial_dip*rad)
        # do the rotation
        xc=X[0]*(sa*sa+ca*ca*cdp)+X[1]*(ca*sa*(1.-cdp))+X[2]*sdp*ca
        yc=X[0]*ca*sa*(1.-cdp)+X[1]*(ca*ca+sa*sa*cdp)-X[2]*sa*sdp
        zc=X[0]*ca*sdp-X[1]*sdp*sa-X[2]*cdp
        # convert back to direction:
        Dir=cart2dir([xc,yc,-zc])
        return [Dir[0],Dir[1]]
        
    def fast_untilt(self, untilting=1):
        """
        performs a partial tilt correction according to untilting (0=geographic, 1=stratigraphic) 
        with a call to a fortran function wrapped with f2py. Much faster
        when taking multiple draws for a bootstrap simulation. 
        """
        return f_untilt.untilt([self.Dec,self.Inc],self.BedStrike,self.BedDip*untilting).tolist()
        
    def di_map(self, untilting=0):
        """
        returns xy coordinates for a stereonet plot
        if untilting is specified then will perform the tilt correction beforehand.
        """
        if untilting==0:
            return dimap(self.Dec,self.Inc)
        else: 
            result=self.untilt(untilting)
            return dimap(result[0],result[1])
          
    def geographic(self):
        return [self.Dec,self.Inc]
        
    def stratigraphic(self):
        return self.untilt()
    
    def at_age(self,age):
        if age>self.SiteAge:
            return [nan,nan]
        else:
            if age<self.FoldMinAge:
                return[self.Dec,self.Inc]
            elif age>self.FoldMaxAge:
                return self.untilt()
            else:
                return self.untilt((age-self.FoldMinAge)/(self.FoldMaxAge-self.FoldMinAge))
                
    def sep_angle(self,otherdir):
        return angle((self.Dec,self.Inc),(otherdir.Dec,otherdir.Inc))[0]
        
    def plot_mean(self,untilting=0, color='k',marker='o',markersize=50,ellipsewidth=3):
        """
        Plots the mean direction at the specified untilting factor (default=geographic), by a 
        Assumes a pre-existing plot such as pretty.pmag_plots.subplot_net
        for no ellipse set ellipsewidth to 0.
        """
        # need this to test inclination
        plotdir=self.untilt(untilting)
        XY=self.di_map(untilting)
        if plotdir[1] < 0:
            plt.scatter(XY[0],XY[1],
                edgecolors=color ,facecolors='white',
                marker=marker,s=markersize,label=self.Name+'untilting '+`untilting`)
        else:
            plt.scatter(XY[0],XY[1],
                edgecolors=color,facecolors=color,
                marker=marker,s=markersize,label=self.Name+'untilting '+`untilting`)
        Xcirc,Ycirc=[],[]
        Da95,Ia95=circ(plotdir[0],plotdir[1],self.a95)
        for k in range(len(Da95)):
            XY=dimap(Da95[k],Ia95[k])
            Xcirc.append(XY[0])
            Ycirc.append(XY[1])
        plt.plot(Xcirc,Ycirc,c=color,linewidth=ellipsewidth)      
                
                    
    def common_mean_boot(self,otherdir,bootsteps=1000, fast=True):
        """
        Tests against another Pmag_Site for a common direction using the bootstrap method. For Synthetic_PMag_Sites where 
        no underlying directions are available, the input to the bootstrap simulations is a direct draw from the fisher 
        distribution using the mean direction parameters, rather than a pseudoselection from specified points. By default 
        uses the wrapped fortran untilt function rather than the pure python, since it is much faster.
        """
        # generate bootstrap distribution of mean directions in xyz
        if fast==True:
            BDI1=dir2cart([[fmean['dec'],fmean['inc']] for fmean in [fisher_mean(self.fast_resample()) for _ in range(bootsteps)]])
            BDI2=dir2cart([[fmean['dec'],fmean['inc']] for fmean in [fisher_mean(otherdir.fast_resample()) for _ in range(bootsteps)]])           
        else:            
            BDI1=dir2cart([[fmean['dec'],fmean['inc']] for fmean in [fisher_mean(self.resample()) for _ in range(bootsteps)]])
            BDI2=dir2cart([[fmean['dec'],fmean['inc']] for fmean in [fisher_mean(otherdir.resample()) for _ in range(bootsteps)]])
        #now want to check if is a pass or a fail - only pass if error bounds overlap in x,y,and z
        bresult=[]
        for b in range(3):
            if (np.percentile(BDI1[:,b],2.5)>np.percentile(BDI2[:,b],97.5) or np.percentile(BDI1[:,b],97.5)<np.percentile(BDI2[:,b],2.5)):
                bresult.append(0)
            else: bresult.append(1)
        if sum(bresult)==3: outcome='Pass'
        else: outcome='Fail'    
        return pd.Series([outcome,self.sep_angle(otherdir)],index=['Outcome','angle'])
        
    def common_mean_MM90(self,otherdir,NumSims=5000, fast=True):
        #largely based on iWatsonV routine of Swanson-Hyell in IPMag    
        cart_1=dir2cart([self.Dec,self.Inc,self.R])
        cart_2=dir2cart([otherdir.Dec,otherdir.Inc,otherdir.R])
        Sw=self.k*self.R+otherdir.k*otherdir.R # k1*r1+k2*r2
        xhat_1=self.k*cart_1[0]+otherdir.k*cart_2[0] # k1*x1+k2*x2
        xhat_2=self.k*cart_1[1]+otherdir.k*cart_2[1] # k1*y1+k2*y2
        xhat_3=self.k*cart_1[2]+otherdir.k*cart_2[2] # k1*z1+k2*z2
        Rw=np.sqrt(xhat_1**2+xhat_2**2+xhat_3**2)
        V=2*(Sw-Rw)
        # keep weighted sum for later when determining the "critical angle" 
        # let's save it as Sr (notation of McFadden and McElhinny, 1990)
        Sr=Sw 
    
        # do monte carlo simulation of datasets with same kappas as data, 
        # *but a common mean*, and take 95th percentile to get Vcrit
        if fast==True:
            Vcrit=np.percentile([vfunc(fisher_mean(self.fast_resample(rotate=False)),fisher_mean(otherdir.fast_resample(rotate=False))) for _ in range(NumSims)],95)
        else:
            Vcrit=np.percentile([vfunc(fisher_mean(self.resample(rotate=False)),fisher_mean(otherdir.resample(rotate=False))) for _ in range(NumSims)],95)

        # equation 18 of McFadden and McElhinny, 1990 calculates the critical
        # value of R (Rwc)
    
        Rwc=Sr-(Vcrit/2)

        # following equation 19 of McFadden and McElhinny (1990) the critical
        # angle is calculated. If the observed angle (also calculated below)
        # between the data set means exceeds the critical angle the hypothesis 
        # of a common mean direction may be rejected at the 95% confidence
        # level. The critical angle is simply a different way to present 
        # Watson's V parameter so it makes sense to use the Watson V parameter
        # in comparison with the critical value of V for considering the test
        # results. What calculating the critical angle allows for is the 
        # classification of McFadden and McElhinny (1990) to be made
        # for data sets that are consistent with sharing a common mean.

        critical_angle=np.degrees(np.arccos(((Rwc**2)-((self.k*self.R)**2)
                                                -((otherdir.k*otherdir.R)**2))/
                                                (2*self.k*self.R*otherdir.k*otherdir.R)))

        if V<Vcrit: 
            outcome='Pass'
            if critical_angle<5: MM90class='A'
            elif critical_angle<10: MM90class='B'
            elif critical_angle<20: MM90class='C'
            else: MM90class='INDETERMINATE'
        else:
            outcome='Fail'
            MM90class='FAIL'
        
        return pd.Series([outcome,V,Vcrit,self.sep_angle(otherdir),critical_angle,MM90class], index=['Outcome','VWatson','Vcrit','angle','critangle','MM90class']) 
        
    def summary(self,untilt='False'):
        if untilt==True:
            outdir,ref_frame=self.untilt(),'tilt'
        else: outdir,ref_frame=[self.Dec,self.Inc],'geo'
        return pd.Series([self.Name,self.SiteAge,ref_frame,outdir[0],outdir[1],self.n,self.k,self.a95],index=['name','age','frame','D','I','n','k','a95'])


class PMagSite(Synthetic_PMagSite):
    """
    A paleomagnetic site with an underlying distribution of sample mean directions. 
    """
    def __init__(self,directions,age=0.,tiltinfo=[],name='PMagSiteDir'):
        """Return object
        directions should be a DataFrame with columns of D,I (in geographic coordinates) - will calculate other parameters on initiation
        tiltinfo should be a DataFrame with columns BedStrike,BedDip,FoldMinAge,FoldMaxAge
        """
        meandir=fisher_mean(directions)
        Synthetic_PMagSite.__init__(self,
            pd.Series([meandir['dec'],meandir['inc'],len(directions),meandir['alpha95'],meandir['k']],
            index=['D','I','n','a95','k']),tiltinfo,name)
        self.dirs=directions 
        
    def resample(self, rotate=True, method='auto'):
        """
        If type is set to 'auto', then it will take a parametric sample from the fisher
        distribution for small n (<25) and a pseudo sample for large n (>25), per Tauxe.
        Can also directly specify which method (e.g. in circumstances where data are
        clearly not fisher distributed) by setting to 'param' or 'pseud'
        Setting rotate to 'False' will produce a fisher
        distribution with the site k, centred on 0,90 (useful for implementation 
        of MM90 fold test, unlikely to be needed otherwise) - this neccessitates 
        parametric sampling (which is fine because MM90 assumes Fisherian data).
        """
        if method=='auto':
            if self.n<25: method='param'
            else: method='pseud'
        if rotate==False: method='param'
        
        if method=='param':
            tempD,tempI=[],[]
            rotmat=np.matrix([dir2cart([self.Dec,90-self.Inc,1.]),dir2cart([self.Dec-90,0,1.]),dir2cart([self.Dec,self.Inc,1.])])
            for i in range(self.n):
                dec,inc=fshdev(self.k)
                if rotate==False:
                    tempD.append(dec)
                    tempI.append(inc)  
                else: 
                    #newdir=cart2dir(dir2cart([dec,inc,1])*rotmat)
                    X=dir2cart([dec,inc,1.]) # get cartesian coordinates of dec,inc
                    # get some sines and cosines of new coordinate system
                    sa,ca= -np.sin((self.Dec)*rad),np.cos((self.Dec)*rad)
                    cdp,sdp= np.cos((90-self.Inc)*rad),np.sin((90-self.Inc)*rad)
                    # do the rotation
                    xc=X[0]*(sa*sa+ca*ca*cdp)+X[1]*(ca*sa*(1.-cdp))+X[2]*sdp*ca
                    yc=X[0]*ca*sa*(1.-cdp)+X[1]*(ca*ca+sa*sa*cdp)-X[2]*sa*sdp
                    zc=X[0]*ca*sdp-X[1]*sdp*sa-X[2]*cdp
                    # convert back to direction:
                    newdir=cart2dir([xc,yc,-zc])
                    tempD.append(newdir[0])
                    tempI.append(newdir[1])
            return np.column_stack((tempD,tempI))
            
        else: return self.dirs[np.random.randint(len(self.dirs),size=len(self.dirs))]
        
    def fast_resample(self, rotate=True, method='auto'):
        """
        When taking a parametric sample of fisher distribution (for n<25 when method is 'auto',
        or when method specifically set to 'param', or rotate set to 'False'), calls fortran 
        function wrapped with f2py. Much faster when taking multiple draws for a bootstrap simulation. 
        Pseudosampling does not call a fortan function but is a fast numpy call that doesn't appear to
        be a bottleneck. 
        """
        if method=='auto':
            if self.n<25: method='param'
            else: method='pseud'
        if rotate==False: method='param'
        
        if method=='param':
            if rotate==True:
                return f_untilt.get_fish(self.Dec,self.Inc,self.k,self.n)
            else:
                return f_untilt.get_fish(0.,90.,self.k,self.n)
                
        else: return self.dirs[np.random.randint(len(self.dirs),size=len(self.dirs))] 
                                              
    def untilt_dirs(self,untilting=1):
        """
        performs a partial tilt correction according to untilting (0=geographic, 1=stratigraphic)
        adapted from the dotilt routine in PMagPy (Lisa Tauxe) 
        """
        rad=np.pi/180. # converts from degrees to radians
        newdirs=[]
        for dir in self.dirs:
            X=dir2cart([dir[0],dir[1],1.]) # get cartesian coordinates of dec,inc
            # get some sines and cosines of new coordinate system
            partial_dip=self.BedDip*untilting
            sa,ca= -np.sin((self.BedStrike+90.)*rad),np.cos((self.BedStrike+90.)*rad)
            cdp,sdp= np.cos(partial_dip*rad),np.sin(partial_dip*rad)
            # do the rotation
            xc=X[0]*(sa*sa+ca*ca*cdp)+X[1]*(ca*sa*(1.-cdp))+X[2]*sdp*ca
            yc=X[0]*ca*sa*(1.-cdp)+X[1]*(ca*ca+sa*sa*cdp)-X[2]*sa*sdp
            zc=X[0]*ca*sdp-X[1]*sdp*sa-X[2]*cdp
            # convert back to direction:
            Dir=cart2dir([xc,yc,-zc])
            newdirs.append([Dir[0],Dir[1]])
        return newdirs
                       
    def fast_untilt_dirs(self,untilting=1):
        """
        performs a partial tilt correction according to untilting (0=geographic, 1=stratigraphic)
        with a call to a fortran function wrapped with f2py. Considerably faster for bootstrap simulations
        """
        return f_untilt.untilt_dirs(self.dirs,self.BedStrike,self.BedDip*untilting).tolist()
                                                        
