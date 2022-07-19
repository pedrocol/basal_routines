from netCDF4 import Dataset
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature
import cartopy.feature as cfeature
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os


def import_depth(ncfile=None,xx=None,tt=None):
    fh = Dataset(ncfile, mode='r')
    #Get data
    data = fh.variables[('e')][tt,:,:,xx]
    fh.close() #Close file
    return data

def import_depth2d(ncfile=None,tt=None):
    fh = Dataset(ncfile, mode='r')
    #Get data
    data = fh.variables[('e')][tt,:]
    fh.close() #Close file
    return data

def import_data(ncfile=None,xx=None,var=None,tt=None):
    fh = Dataset(ncfile, mode='r')
    #Get data
    data = fh.variables[(var)][tt,:,:,xx]
    fh.close() #Close file
    return data

def import_data2d(ncfile=None,var=None,tt=None):
    fh = Dataset(ncfile, mode='r')
    #Get data
    data = fh.variables[(var)][tt,:,:]
    fh.close() #Close file
    return data

def import_lat(ncfile=None,xx=None):
    fh = Dataset(ncfile, mode='r')
    #Get data
    data = fh.variables[('geolat')][:,xx]
    fh.close() #Close file
    return data

def import_lat_0(ncfile=None):
    fh = Dataset(ncfile, mode='r')
    #Get data
    data = fh.variables[('geolat')][:,0]
    fh.close() #Close file
    return data

def import_lonlat(ncfile=None):
    fh = Dataset(ncfile, mode='r')
    #Get data
    lat = fh.variables[('geolat')][:,:]
    lon = fh.variables[('geolon')][:,:]
    fh.close() #Close file
    return lon,lat

def import_data3d_TT(ncfile=None,variable=None,tt=None):
    fh = Dataset(ncfile, mode='r')
    fh.set_auto_mask(False)
    #Get data
    data = fh.variables[str(variable)][tt,:,:,:]
    fh.close() #Close file
    return data


def plot_isomip_section(var = "thetao",xx0=0,tt=0):
    plt.figure()
    root = "/home/552/pc5520/MOM6/control/ISOMIP/layer/archive/output011/"
    ncfile = root+"ocean_daily_z.nc"
    depthfile = root+"ocean_daily_z.nc"
    latfile = root+"ocean_static.nc"
    depth = import_depth(ncfile=depthfile,xx=xx0,tt=tt)
    lat = import_lat(ncfile=latfile,xx=xx0)
    dataplot = import_data(ncfile=ncfile,xx=xx0,var=var,tt=tt)
    if var == "thetao": vmax,vmin,contourstep=-1.9,-2.1,0.01
    if var == "vo": vmax,vmin,contourstep=0.1,-0.1,0.01
    nb_contours = np.ceil((vmax-vmin)/contourstep)
    contours = [round(x , 3) for x in np.linspace(vmin,vmax,int(nb_contours+1))]
    cmap=cm.jet
    #ax1.set_facecolor('xkcd:gray')
    xx = lat
    if var == "vo": dataplot = dataplot[:,:-1]
    yy = depth[:-1,:]
    z_levels = yy.shape[0]
    yi = np.linspace(0,5000,z_levels) #Dummy, just for have z_levels
    xi, yyi = np.meshgrid(xx, yi) #We repeat the horizontal part of the grid by the number of zlevels
    CS = plt.contourf(xi,yy,dataplot,contours,vmin=vmin,vmax=vmax,\
                     cmap = cmap,extend='both',corner_mask=True)
    cbar = plt.colorbar(CS)
    plt.ylim(ymin=-900,ymax=0)
    plt.show(block=False)

def plot_isomip_section_2dave(var = "vo"):
    root = "/home/552/pc5520/MOM6/control/ISOMIP/layer/archive/"
    ncfile = root+"avey_avet_ocean_daily_000_011.nc"
    depthfile = ncfile
    latfile = root+"output000/ocean_static.nc"
    tt = 0
    depth = import_depth2d(ncfile=depthfile,tt=tt)
    lat = import_lat_0(ncfile=latfile)
    dataplot = import_data2d(ncfile=ncfile,var=var,tt=tt)
    if var == "thetao": vmax,vmin,contourstep=-1.9,-2.1,0.01
    if var == "vo": vmax,vmin,contourstep=0.01,-0.01,0.001
    nb_contours = np.ceil((vmax-vmin)/contourstep)
    contours = [round(x , 4) for x in np.linspace(vmin,vmax,int(nb_contours+1))]
    cmap=cm.jet
    #ax1.set_facecolor('xkcd:gray')
    xx = lat
    if var == "vo": dataplot = dataplot[:,:-1]
    yy = depth[:-1,:]
    z_levels = yy.shape[0]
    yi = np.linspace(0,5000,z_levels) #Dummy, just for have z_levels
    xi, yyi = np.meshgrid(xx, yi) #We repeat the horizontal part of the grid by the number of zlevels
    CS = plt.contourf(xi,yy,dataplot,contours,vmin=vmin,vmax=vmax,\
                     cmap = cmap,extend='both',corner_mask=True)
    cbar = plt.colorbar(CS)
    plt.ylim(ymin=-900,ymax=0)
    plt.show(block=False)

def plot_data_along_line(tt=0):
    fig, ax = plt.subplots(figsize=(15,5))
    ax.set_facecolor('xkcd:gray')
    # Import data
    root = "/home/552/pc5520/MOM6/control/AMERY/archive/output004/"
    datafile   = root+"ocean_daily_model.nc"
    lonlatfile   = root+"ocean_static.nc"
    depth = import_data3d_TT(ncfile=datafile,variable="e",tt=tt)
    depth = depth[:-1,:,:]
    temp = import_data3d_TT(ncfile=datafile,variable="thetao",tt=tt)
    lon,lat = import_lonlat(ncfile=lonlatfile)
    # Define segment
    i_main,j_main = [6,20,32,41,46],[9,44,73,125,173] #side
    #i_main,j_main = [6,26,76,107],[9,50,132,173] #center
    iii,jjj = generate_line(i_main,j_main)
    iii = iii.astype(int)
    jjj = jjj.astype(int)
    temp_plot = temp[:,jjj,iii]
    depth_plot = depth[:,jjj,iii]
    lon,lat = lon[jjj,iii],lat[jjj,iii]
    temp_plot = np.ma.masked_where(temp_plot==1e+20,temp_plot)
    depth_plot = np.ma.masked_where(depth_plot==1e+20,depth_plot)
    #Plot setup
    vmax,vmin,contourstep=-1.6,-2,0.01
    nb_contours = np.ceil((vmax-vmin)/contourstep)
    contours = [round(x , 4) for x in np.linspace(vmin,vmax,int(nb_contours+1))]
    cmap=cm.Blues_r
    #Build grid to plot
    dataplot = temp_plot
    yy = depth_plot
    z_levels = yy.shape[0]
    yi = np.linspace(0,5000,z_levels)
    xi, yyi = np.meshgrid(lat, yi)
    CS = plt.contourf(xi,yy,dataplot,contours,vmin=vmin,vmax=vmax,\
                     cmap = cmap,extend='both',corner_mask=True)
    cbar = plt.colorbar(CS)
    plt.ylim(ymin=-3300,ymax=0)
    plt.show(block=False)



def generate_line(i_main=None,j_main=None):
    ii_f,jj_f = [],[]
    xshape = len(i_main)
    for seg in range(1,xshape):
        ii1,jj1 = i_main[seg-1],j_main[seg-1]
        ii2,jj2 = i_main[seg],j_main[seg]
        if (jj2-jj1) > (ii2-ii1):
           yy = np.linspace(jj1,jj2,jj2-jj1+1)
           b = ii1
           m = (ii1-ii2)/(jj1-jj2)
           xx = np.floor(m*(yy-jj1) + b)
        if (ii2-ii1) > (jj2-jj1):
           xx = np.linspace(ii1,ii2,ii2-ii1+1)
           b = jj1
           m = (jj1-jj2)/(ii1-ii2)
           yy = np.floor(m*(xx-ii1) + b)
        if ii2 == ii1:
           yy = np.linspace(jj1,jj2,jj2-jj1+1)
           xx = yy * 0 + ii1
        if jj2 == jj1:
           xx = np.linspace(ii1,ii2,ii2-ii1+1)
           yy = jj * 0 + jj1
        ii_f = np.append(ii_f,xx)
        jj_f = np.append(jj_f,yy)
    return ii_f,jj_f

