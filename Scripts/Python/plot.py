import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.path as mpath
from matplotlib import ticker
from netCDF4 import Dataset
#Cartopy
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature
import cartopy.feature as cfeature

def fmt(x, pos):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def plot_salt(variable="so",tt=None,ii=None,jj=None):
    #plt.figure()
    ncfile = "/home/552/pc5520/MOM6/control/AMERY/archive/output000/ocean_daily_z.nc"
    fh = Dataset(ncfile, mode='r')
    fh.set_auto_mask(False)
    #data = fh.variables[str(variable)][tt,:,jj,ii]
    data1 = fh.variables[str(variable)][tt,:,172,48]
    data2 = fh.variables[str(variable)][tt,:,36,16]
    vert = fh.variables["e"][tt,1:,jj,ii] 
    vert1 = fh.variables["e"][tt,1:,172,48]
    vert2 = fh.variables["e"][tt,1:,36,16]
    fh.close() #Close file
    #data = np.ma.masked_where(data==1e+20,data)
    data1 = np.ma.masked_where(data1==1e+20,data1)
    data2 = np.ma.masked_where(data2==1e+20,data2)
    plt.plot(-vert1,data1,color='b')
    plt.plot(-vert2,data2,color='g')
    plt.xlim(xmin=0,xmax=2500)
    plt.ylim(ymin=33.5,ymax=35)
    plt.show(block=False)
