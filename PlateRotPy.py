import numpy as np
import re

def load_points(pointfile, point_file_format=[('Lat', np.float), ('Long', np.float)]):
    #if want more than Lat and Long, have to specify
    #typical point_file_format=[('Lat', np.float), ('Long', np.float), ('PCode', np.int32), ('ReconAge', np.float), ('FeatAge', np.float), ('Chron', np.str_,16)]
    data=np.loadtxt(pointfile, dtype=point_file_format)
    return data

def load_poles(polefile):
    pole_file_format=[('Rotation', np.str_,40),('Chron', np.str_,16),('Age', np.float),('Lat', np.float),('Long', np.float),('Rot', np.float), ('kappa', np.float),('a', np.float), ('b', np.float),('c', np.float), ('d', np.float),('e', np.float), ('f', np.float),('n', np.int32),('Segs', np.str_,4),('Plates', np.str_,2),('DOF', np.int32),('Ref', np.str_,100)]
    data=np.loadtxt(polefile, dtype=pole_file_format)
    return data

def plt_line(longs,lats,m,colour='black'):
	x,y= m(longs,lats)
	#note: zorder is a very useful thing. Higher number=higher plot level.
	m.plot(x,y,color=colour, zorder=3)

def plt_segs_file(file,m,colour='black',order=1):
	f=open(file, 'r')
	ridgeseg=[]
	flag=0
	for line in f:
		if re.match('^>', line):
			if flag==1: 
				seglat=[]
				seglong=[]
				for data in ridgeseg:
					data=data.rstrip('\n')
					seglong.append(float(data.split("\t")[0]))
					seglat.append(float(data.split("\t")[1]))
				x, y = m(seglong,seglat)
				#this converts lat and long into basemap coordinates
				m.plot(x,y, color=colour, zorder=order)				
			ridgeseg = []
			flag=1
		else:
			ridgeseg.append(line)
	f.close()
