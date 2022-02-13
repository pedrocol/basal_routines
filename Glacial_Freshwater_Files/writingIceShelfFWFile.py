import netCDF4
import sys
sys.path.append('/home/pedro/DEV/NEMO_Glacial_Freshwater_Files')
import numpy as np
from NEMO_ICB_File import *


#List of Ice shelves to consider#
listShelf=['lbc','fris','brl','jf','ar','ne','ais','w','sha','van','tot','mu','por','ade','mer','nin','coo','ren','dry','ris','sul','lan'
            ,'getz','cd','thw','pig','cosg','abb','ven','geo','wor']

listSectors=['westIndian','eastIndian','rossSea','amundsen','bellingshausen','weddell']

#listShelf=['pig']
#First coastal point, needed to use extractCoastNew. It is an ocean point with land at (x,y-1)
# 1/4 xinit=1 yinit=124
# 1 deg xinit=1 yinit=31
xinit=1
yinit=124


#Open needed Files
ncfile = netCDF4.Dataset('./topog.nc','r')
varBathy = np.array(ncfile.variables['depth'])[:,:]
ncfile.close()

varBathy[0,:] = 0 #First values masked
varBathy[1,:] = 0 #First values masked

ncfile = netCDF4.Dataset('./ocean_hgrid.nc','r')
#varLon = np.array(ncfile.variables['glamt'])[:,:]
#varLat = np.array(ncfile.variables['gphit'])[:,:]
#e1t = np.array(ncfile.variables['e1t'])[:,:]
#e2t = np.array(ncfile.variables['e2t'])[:,:]
varLon = np.array(ncfile.variables['x'])[:,:]
varLat = np.array(ncfile.variables['y'])[:,:]
e1t = np.array(ncfile.variables['dx'])[:,:]
e2t = np.array(ncfile.variables['dy'])[:,:]
ncfile.close()

varLon = varLon[0::2,0::2][:-1,:-1]
varLat = varLat[0::2,0::2][:-1,:-1]
#1/4 Deg
e1t = e1t[0::2,0::2][:-1,:-1]
e2t = e2t[0::2,0::2][:-1,:-1]
#1deg
#e1t = e1t[0::2,0::2][:-1,:]
#e2t = e2t[0::2,0::2][:-1,:]
###########
area=e1t*e1t

[ydim,xdim]=varBathy.shape

#######1Clean the Bathymetry, 0->land, 1000->ocean
var=cleanBathy(varBathy,1000)
########

######## Create rmpty fresh flux variable (units:kg/m2/s)
FreshwaterFlux=np.zeros([12,ydim,xdim])
varGL=np.zeros([12,ydim,xdim])
varFront=np.zeros([12,ydim,xdim])
########

######## Extract the list of consecutive coastal points
listP=extractCoast(var,xinit,yinit)
#List of X and Y values in the grid referencial
xP=list(zip(*listP))[0]
yP=list(zip(*listP))[1]
#####

####### Fill freshwater flux variable with data from the list of shelves and sectors   
createIceShelfFluxFile(FreshwaterFlux,listShelf,xP,yP,varLon,varLat,area)
#createIceShelfFluxFile(FreshwaterFlux,listSectors,xP,yP,varLon,varLat,area)

####### Fill GL and claving front depths with data from the list of shelves 
writeDepthFwFluxes(varGL,varFront,varBathy,listShelf,xP,yP,varLon,varLat)

#FreshwaterFlux writen in kg/m2/s as needed by NEMO

ncfile = netCDF4.Dataset('./basal_fw.nc','w')
ncfile.createDimension('time')
ncfile.createDimension('ny', varLon.shape[0])
ncfile.createDimension('nx', varLon.shape[1])
time = ncfile.createVariable('time','f4',('time'))
time.lon_name='time'
time.cartesian_axis = "T"
time.axis = "T"
time.calendar_type = "noleap"
time.modulo = " "
time.standard_name = "time"
time.climatology = "climatology_bounds"
time.units='months since 0001-01-01 00:00:00'
dMax = ncfile.createVariable('sodepmax_isf','f4',('time','ny','nx'))
dMin = ncfile.createVariable('sodepmin_isf','f4',('time','ny','nx'))
fw = ncfile.createVariable('meltingFlux','f4',('time','ny','nx'))
fw[:,:,:]=FreshwaterFlux[:,:,:]
dMax[:,:,:]=varGL[:,:,:]
dMin[:,:,:]=varFront[:,:,:]
ncfile.close()
