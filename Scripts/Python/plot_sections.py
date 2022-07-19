from netCDF4 import Dataset
#Cartopy
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


def import_data_cut(ncfile=None,variable=None):
    fh = Dataset(ncfile, mode='r')
    #Get data
    data = fh.variables[str(variable)][:,:,:,0]
    vert_lev = fh.variables[('st_ocean')][:]
    lat = fh.variables[('yt_ocean_sub01')][:]
    lon = fh.variables[('xt_ocean_sub01')][:]
    fh.close() #Close file
    return data,vert_lev,lat,lon

def loop_plot_cut_section(var=None):
    directory = "/scratch/e14/pc5520/Plots/"+str(var)+"/"
    if not os.path.exists(directory):
       os.makedirs(directory)
    string = "cat_rregionocean_daily_3d_i1123_"+str(var)
    for T in range(0,364):
        print(T)
        plot_cut_section(tt=T,var=var)
        print(directory+str(string)+"_"+"{0:03d}".format(T)+".png")
        plt.savefig(directory+str(string)+"_"+"{0:03d}".format(T)+".png")
        plt.close()

def import_data3d_TT(ncfile=None,variable=None,tt=None):
    fh = Dataset(ncfile, mode='r')
    fh.set_auto_mask(False)
    #Get data
    data = fh.variables[str(variable)][tt,:,:,:]
    fh.close() #Close file
    return data

def plot_cut_section(tt=None,var=None):
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(7,7))
    fs = 10
    if var == "salt": div = 7
    if var == "temp": div = 10
    #### acces-om2-01-GAM001 ####
    #Plot setup
    root = "/scratch/e14/pc5520/OUTPUT/acces-om2-01-GAM001/extract/y2150/sections/"+str(var)+"/1123/"
    ncfile = root+"cat_rregionocean_daily_3d_i1123_"+str(var)+"_2150.nc"
    data_plot,vert_lev,lat,lon = import_data_cut(ncfile=ncfile,variable=var)
    xx, yy = np.meshgrid(lat, vert_lev)    
    ax1.set_facecolor('xkcd:gray')
    data_t = data_plot[tt,:,:]
    if var == "temp": data_t = data_t - 273
    if var == "salt": vmax,vmin,contourstep=35,33.2,0.05
    if var == "temp": vmax,vmin,contourstep=2,-2.5,0.05
    nb_contours = np.ceil((vmax-vmin)/contourstep)
    contours = [round(x , 4) for x in np.linspace(vmin,vmax,int(nb_contours+1))]
    cmap=cm.jet
    CS = ax1.contourf(xx,yy,data_t,contours,vmin=vmin,vmax=vmax,\
                     cmap = cmap,extend='both',corner_mask=True)
    cbar = plt.colorbar(CS,ticks=np.linspace(vmin,vmax,div),ax=ax1)
    cbar.ax.tick_params(labelsize=fs)
    #### acces-om2-01-GAM001 ####
    #Plot setup
    root = "/scratch/e14/pc5520/OUTPUT/acces-om2-01-GPC001/extract/y2150/sections/"+str(var)+"/1123/"
    ncfile = root+"cat_rregionocean_daily_3d_i1123_"+str(var)+"_2150.nc"
    data_plot,vert_lev,lat,lon = import_data_cut(ncfile=ncfile,variable=var)
    xx, yy = np.meshgrid(lat, vert_lev)
    ax2.set_facecolor('xkcd:gray')
    data_t = data_plot[tt,:,:]
    if var == "temp": data_t = data_t - 273
    if var == "salt": vmax,vmin,contourstep=35,33.2,0.05
    if var == "temp": vmax,vmin,contourstep=2,-2.5,0.05
    nb_contours = np.ceil((vmax-vmin)/contourstep)
    contours = [round(x , 4) for x in np.linspace(vmin,vmax,int(nb_contours+1))]
    cmap=cm.jet
    CS = ax2.contourf(xx,yy,data_t,contours,vmin=vmin,vmax=vmax,
                     cmap = cmap,extend='both',corner_mask=True)
    cbar = plt.colorbar(CS,ticks=np.linspace(vmin,vmax,div),ax=ax2)
    cbar.ax.tick_params(labelsize=fs)
    #### Diff ####
    #Plot setup
    root = "/scratch/e14/pc5520/OUTPUT/acces-om2-01-GPC001/extract/y2150/sections/"+str(var)+"/1123/"
    ncfile = root+"diff_cat_rregionocean_daily_3d_i1123_"+str(var)+"_2150.nc"
    data_plot,vert_lev,lat,lon = import_data_cut(ncfile=ncfile,variable=var)
    xx, yy = np.meshgrid(lat, vert_lev)
    ax3.set_facecolor('xkcd:gray')
    data_t = data_plot[tt,:,:]
    if var == "salt": vmax,vmin,contourstep=0.3,-0.3,0.05
    if var == "temp": vmax,vmin,contourstep=0.75,-0.75,0.05
    nb_contours = np.ceil((vmax-vmin)/contourstep)
    contours = [round(x , 4) for x in np.linspace(vmin,vmax,int(nb_contours+1))]
    cmap=cm.seismic
    CS = ax3.contourf(xx,yy,data_t,contours,vmin=vmin,vmax=vmax,
                     cmap = cmap,extend='both',corner_mask=True)
    cbar = plt.colorbar(CS,ticks=np.linspace(vmin,vmax,7),ax=ax3)
    cbar.ax.tick_params(labelsize=fs)
    #Last details
    ax1.set_ylim([0, 3000])
    ax2.set_ylim([0, 3000])
    ax3.set_ylim([0, 3000])
    plt.autoscale(True,'x',True)
    plt.xlim(xmin=-79,xmax=-65)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    ax3.invert_yaxis()
    ax2.set_ylabel("Depth (m)",fontsize=fs)
    #plt.title("Longitue  = "+"{:.2f}".format(lon[0]),fontsize=fs)
    plt.xlabel("Latitude",fontsize=fs)
    plt.show(block=False)



