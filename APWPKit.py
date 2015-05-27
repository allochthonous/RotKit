import numpy as np
import pandas as pd
import PlateRotKit as rotkit

def reconstruct_DI(forwardrots,point):
    #for a given set of rotations for specified point, find 'preicted D & I of ChRm in present day.
    backwardsrots=rotkit.invert_covrots(forwardrots)
    pastpoints=pd.concat([point,rotkit.rotate_points(point,backwardsrots)], ignore_index=True)
    npole=rotkit.make_points([90],[0],[0],[0])
    vgps=pd.concat([npole, rotkit.rotate_points(npole,forwardrots)], ignore_index=True)
    paleoI=np.arctan(2*np.tan(pastpoints.PointLat*np.pi/180))*180/np.pi
    pp=np.sin((90-vgps.PointLat)*np.pi/180)
    sitelong=point.iloc[0]['PointLong']
    dphi=np.sin((vgps.PointLong-sitelong)*np.pi/180)
    pm=np.sin((90-pastpoints.PointLat)*np.pi/180)
    paleoD=np.arcsin(pp*dphi/pm)*180/np.pi
    predDI=pd.DataFrame(np.column_stack((paleoD,paleoI)), columns=['PredDec','PredInc'])
    return pd.concat([pastpoints, predDI],axis=1),vgps