import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.path as mpath
from netCDF4 import Dataset
#Cartopy
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature
import cartopy.feature as cfeature

def import_data1d(ncfile=None,variable=None):
    fh = Dataset(ncfile, mode='r')
    fh.set_auto_mask(False)
    #Get data
    data = fh.variables[str(variable)][:]
    fh.close() #Close file
    return data

def import_data2d(ncfile=None,variable=None):
    fh = Dataset(ncfile, mode='r')
    fh.set_auto_mask(False)
    #Get data
    data = fh.variables[str(variable)][0,:,:]
    fh.close() #Close file
    return data

def import_data3d(ncfile=None,variable=None):
    fh = Dataset(ncfile, mode='r')
    fh.set_auto_mask(False)
    #Get data
    data = fh.variables[str(variable)][:,:,:]
    fh.close() #Close file
    return data

def import_hgrid(ncfile=None):
    fh = Dataset(ncfile, mode='r')
    fh.set_auto_mask(False)
    #Get data
    datalon = fh.variables["geolon_t"][:,:]
    datalat = fh.variables["geolat_t"][:,:]
    fh.close() #Close file
    return datalon,datalat

def plot2d(ncfile=None,var=None,level="3446"):
    #Data
    plt.figure()
    fs = 10
    root = "/scratch/e14/pc5520/OUTPUT/acces-om2-01-GPC001/extract/y2150/diffs_vertave/"+var+"/"+level+"/"
    ncfile = root+"avet_avek_ncdiff_rregionocean_daily_3d_k"+level+"_"+var+"_y2150.nc"
    gridfile = "/home/552/pc5520/access-om2/control/01deg_jra55v13_ryf9091_rerun_for_easterlies/archive/rregionocean_ocean_grid.nc"
    data_plot = import_data2d(ncfile,var)
    #Axes
    projection = ccrs.PlateCarree(central_longitude=0)
    #projection = ccrs.Orthographic(180,-90)
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    #Grid
    lon,lat = import_hgrid(ncfile=gridfile)
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.475, 0.475], 0.475
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform = ax.transAxes)
    #Plot
    if var == "salt": vmin,vmax,contourstep,rounds=-0.3,0.3,0.01,4
    if var == "temp": vmin,vmax,contourstep,rounds=-0.75,0.75,0.025,4
    cmap=cm.seismic
    nb_contours = int(np.ceil((vmax-vmin)/contourstep))
    contours = [round(x , rounds) for x in np.linspace(vmin,vmax,nb_contours+1)]
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    CS = ax.contourf(lon,lat,data_plot,contours,\
                     cmap = cmap,norm=norm,vmin=vmin,vmax=vmax,transform=projection,extend='both')
    cbar = plt.colorbar(CS,ticks=np.linspace(vmin,vmax,7))
    cbar.ax.tick_params(labelsize=fs)
    #plt.ylim(ymax=-55,ymin=-78)
    plt.show(block=False)
