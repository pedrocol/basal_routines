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

def import_data2d_T(ncfile=None,variable=None):
    fh = Dataset(ncfile, mode='r')
    fh.set_auto_mask(False)
    #Get data
    data = fh.variables[str(variable)][:,:]
    fh.close() #Close file
    return data

def import_data2d_TT(ncfile=None,variable=None,tt=None):
    fh = Dataset(ncfile, mode='r')
    fh.set_auto_mask(False)
    #Get data
    data = fh.variables[str(variable)][tt,:,:]
    fh.close() #Close file
    return data

def import_data3d(ncfile=None,variable=None):
    fh = Dataset(ncfile, mode='r')
    fh.set_auto_mask(False)
    #Get data
    data = fh.variables[str(variable)][:,:,:]
    fh.close() #Close file
    return data

def import_data3d_TT(ncfile=None,variable=None,tt=None):
    fh = Dataset(ncfile, mode='r')
    fh.set_auto_mask(False)
    #Get data
    data = fh.variables[str(variable)][tt,:,:,:]
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

def plot2d(var=None,level="3446",diff="no"):
    #Data
    fs = 10
    root = "/scratch/e14/pc5520/OUTPUT/acces-om2-01-GPC002/extract/y2150/diffs/"+var+"/"+level+"/"
    ncfile = root+"cut_ncdiff_avet_cat_avek_k"+str(level)+"_"+var+".nc"
    gridfile = "/home/552/pc5520/access-om2/control/01deg_jra55v13_ryf9091_rerun_for_easterlies/archive/rregionocean_ocean_grid.nc"
    data_plot = import_data2d(ncfile,var)
    #Axes
    projection = ccrs.PlateCarree(central_longitude=0)
    #projection = ccrs.Orthographic(180,-90)
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    ax.set_facecolor('xkcd:gray')
    #Grid
    lon,lat = import_hgrid(ncfile=gridfile)
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.475, 0.475], 0.475
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform = ax.transAxes)
    #Plot
    if diff == "yes":
       if var == "salt": vmin,vmax,contourstep,rounds,div=-0.15,0.15,0.01,4,7
       if var == "temp": vmin,vmax,contourstep,rounds,div=-0.75,0.75,0.025,4,7
    else:
       if var == "salt": vmin,vmax,contourstep,rounds,div=34.3,34.8,0.05,4,10
       if var == "temp": vmin,vmax,contourstep,rounds,div=-2,2.8,0.4,4,12
    cmap=cm.seismic
    nb_contours = int(np.ceil((vmax-vmin)/contourstep))
    contours = [round(x , rounds) for x in np.linspace(vmin,vmax,nb_contours+1)]
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    CS = ax.contourf(lon,lat,data_plot,contours,\
                     cmap = cmap,norm=norm,vmin=vmin,vmax=vmax,transform=projection,extend='both')
    cbar = plt.colorbar(CS,ticks=np.linspace(vmin,vmax,div))
    cbar.ax.tick_params(labelsize=fs)
    plt.ylim(ymax=-55,ymin=-78)
    plt.show(block=False)

def plot_bathy():
    #Data
    fs = 10
    ncfile = "/home/552/pc5520/access-om2/control/01deg_jra55v13_ryf9091_rerun_for_easterlies/archive/rregionocean_ocean_grid.nc"
    gridfile = "/home/552/pc5520/access-om2/control/01deg_jra55v13_ryf9091_rerun_for_easterlies/archive/rregionocean_ocean_grid.nc"
    var = "ht"
    fh = Dataset(ncfile, mode='r')
    fh.set_auto_mask(False)
    data_plot = fh.variables[str(var)][:,:]
    fh.close() #Close file
    #Axes
    projection = ccrs.PlateCarree(central_longitude=0)
    #projection = ccrs.Orthographic(180,-90)
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    #ax.set_facecolor('xkcd:gray')
    #Grid
    lon,lat = import_hgrid(ncfile=gridfile)
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.475, 0.475], 0.475
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform = ax.transAxes)
    #Plot
    contours = [0]
    colors = ['k']
    CS = ax.contour(lon,lat,data_plot,contours,colors=colors,transform=projection)
    xx = -167.65
    plt.plot([xx,xx],[-63,-80],transform=projection)
    #plt.ylim(ymax=-55,ymin=-78)
    plt.show(block=False)

def plot_bathy_cart(varp='bathy',conf=None,tt=None):
    fs=12
    ncfile = "/home/552/pc5520/AMERY/AMERY_10/MOM_bathy.nc"
    #datafile = "/home/552/pc5520/MOM6/control/AMERY/archive/"+str(conf)+"/output000/ocean_hourly_model.nc"
    datafile = "/home/552/pc5520/MOM6-examples/AMERY/archive/output006_ocean.nc"
    datafile = "/home/552/pc5520/MOM6-examples/AMERY/archive/output006/19910601.ice_daily.nc"
    lon = import_data2d_T(ncfile=ncfile,variable='nav_lon')
    lat = import_data2d_T(ncfile=ncfile,variable='nav_lat')
    msk = import_data2d_T(ncfile=ncfile,variable='msk')
    bathy = import_data2d_T(ncfile=ncfile,variable='depth')
    isf = import_data2d_T(ncfile=ncfile,variable='isf_draft')
    if varp == 'isf':   var = import_data2d_T(ncfile=ncfile,variable='isf_draft')
    if varp == 'bathy': var = import_data2d_T(ncfile=ncfile,variable='depth')
    if varp in ['tfreeze',"melt","sithick"]: var = import_data2d_TT(ncfile=datafile,variable=varp,tt=tt)
    if varp == 'bathydiff': var = bathy - isf
    if varp in ['uo']: 
       #var = max_vo(ncfile=datafile,tt=tt,var='uo')
       #var = var[:,1:]
       var=import_data3d_TT(ncfile=datafile,variable=varp,tt=tt)
       var = var[41,:,1:]
       var = np.ma.masked_where(var==1e+20,var)
    if varp in ['vo']: 
       var=import_data3d_TT(ncfile=datafile,variable=varp,tt=tt)
       var = var[41,1:,:]
       var = np.ma.masked_where(var==1e+20,var)
    if varp in ["thetao",'e','h']:
       var=import_data3d_TT(ncfile=datafile,variable=varp,tt=tt)
       var = var[0,:,:]
       var = np.ma.masked_where(var==1e+20,var)
    mask = np.ma.masked_where(msk==0,msk)
    if varp == 'bathydiff':
       mask = np.ma.masked_where(isf==0,isf) * 0 + 1
       var = var * mask
    if varp not in ['uo','e','vo'] : var = var * mask
    projection = ccrs.PlateCarree(central_longitude=0)
    plt.figure(figsize=(8,10))
    ax = plt.axes(projection=projection)
    land = cfeature.NaturalEarthFeature('physical','land','50m',edgecolor='blue',facecolor=cfeature.COLORS['land'],linewidth=0.5)
    #ax.add_feature(land)
    ax.set_autoscale_on(False)
    ax.set_extent([lon[0,lon.shape[1]-1], lon[0,0], lat[0,lon.shape[1]-1], lat[lon.shape[0]-1,lon.shape[1]-1]],projection)
    ax.set_aspect('auto')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_facecolor('xkcd:gray')
    if varp == 'isf':    vmin,vmax,contourstep,rounds,div=0,2100,0.01,4,7
    if varp == 'bathy': vmin,vmax,contourstep,rounds,div=0,2700,0.01,4,7
    if varp == 'bathydiff': vmin,vmax,contourstep,rounds,div=0,1500,0.01,4,7
    if varp == 'uo': vmin,vmax,contourstep,rounds,div=-1e-3,1e-3,1e-4,6,7
    if varp == 'vo': vmin,vmax,contourstep,rounds,div=-0.1,0.1,0.01,6,6
    if varp == 'thetao': vmin,vmax,contourstep,rounds,div=-2.5,-1.5,0.1,4,5
    if varp == 'e': vmin,vmax,contourstep,rounds,div=-2100,0,0.1,4,7
    if varp == 'h': vmin,vmax,contourstep,rounds,div=0.001,50,0.1,4,7
    if varp == 'tfreeze': vmin,vmax,contourstep,rounds,div=-3.3,-1.5,0.1,4,7
    if varp == 'melt': vmin,vmax,contourstep,rounds,div=0,1.4,1,4,8
    if varp == "sithick": vmin,vmax,contourstep,rounds,div=0,0.4,0.01,4,5
    nb_contours = 50.
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    contourstep = (vmax-vmin)/nb_contours
    contours = np.arange(vmin,vmax+contourstep,contourstep)
    my_cmap=cm.jet
    CS = ax.contourf(lon,lat,var,contours,cmap = my_cmap,norm=norm,vmin=vmin,vmax=vmax,interpolation="none",transform=projection,extend='both')
    formats='%1.4g'
    #if varp == 'uo': formats=ticker.ScalarFormatter(useMathText=True);formats.set_powerlimits((0, 0))
    if varp == 'uo': formats=ticker.FuncFormatter(fmt)
    ticks = np.linspace(vmin,vmax,div)
    if varp == 'tfreeze': ticks = [-3.3,-3.0,-2.7,-2.4,-2.1,-1.8,-1.5]
    #if varp == 'uo': ticks = [1e-3,5e-4,0,-5e-4,-1e-3]
    cbar=plt.colorbar(CS,format=formats,ticks=ticks,shrink=0.7)
    cbar.ax.tick_params(labelsize= fs)
    #if varp in ['uo','theta']: 
    CS1 = plt.contour(lon,lat,isf,colors='k',transform=projection)
    #else:
    #CS1 = plt.contour(lon,lat,bathy,(250,500,1000),colors='k',transform=projection)
    ax.set_xticks([68,70,72,74,76,78])
    ax.set_yticks([-71,-69,-67])
    #Plot lines
    section = "none"
    if section == "center":
       plt.plot((lon[9,6],lon[50,24]),(lat[9,6],lat[50,24]),color='w',lw=4,transform=projection)
       plt.plot((lon[50,24],lon[132,71]),(lat[50,24],lat[132,71]),color='w',lw=4,transform=projection)
       plt.plot((lon[132,71],lon[173,107]),(lat[132,71],lat[173,107]),color='w',lw=4,transform=projection)
    if section == "side":
       plt.plot((lon[9,6],lon[44,20]),(lat[9,6],lat[44,20]),color='w',lw=4,transform=projection)
       plt.plot((lon[44,20],lon[73,32]),(lat[44,20],lat[73,32]),color='w',lw=4,transform=projection)
       plt.plot((lon[73,32],lon[125,41]),(lat[73,32],lat[125,41]),color='w',lw=4,transform=projection)
       plt.plot((lon[125,41],lon[173,46]),(lat[125,41],lat[173,46]),color='w',lw=4,transform=projection)
    #if varp == 'uo': plt.title("Zonal Velocity - Hour 0 (Day 1)")
    #plt.show(block=False)
    plt.savefig('/home/552/pc5520/Plots/plot.png')

def get_uo_pge():
    vfile="/home/552/pc5520/MOM6/control/AMERY/archive/AMERY-GPC002/ocean_hourly_z_000_007.nc"
    ttall = 119
    vo_all_max = np.zeros(ttall-1)
    vo_all_ave = np.zeros(ttall-1)
    for tt in range(0,ttall-1):
        vo=import_data3d_TT(ncfile=vfile,variable='uo',tt=tt)
        vo = np.ma.masked_where(vo==1e+20,vo)
        voabs = abs(vo)
        voave = voabs.mean()
        vomax = voabs.max()
        vo_all_ave[tt] = voave
        vo_all_max[tt] = vomax
        print(tt,vomax,voave)
    return vo_all_max,vo_all_ave

def max_vo(ncfile=None,tt=None,var=None):
    vfile=ncfile
    ttall = 240
    vo=import_data3d_TT(ncfile=vfile,variable=var,tt=tt)
    vo = np.ma.masked_where(vo==1e+20,vo)
    vomax = vo[0,:,:]*0
    for jj in range(0,vo.shape[1]-1):
        for ii in range(0,vo.shape[2]-1):
            if vo[:,jj,ii].max() > abs(vo[:,jj,ii].min()) : vomax[jj,ii] = vo[:,jj,ii].max()
            if vo[:,jj,ii].max() <= abs(vo[:,jj,ii].min()) : vomax[jj,ii] = vo[:,jj,ii].min()
    return vomax

def plot_series(array=None):
    plt.figure()
    #array = array*10**4
    plt.plot(array)
    plt.xlabel('Hours')
    plt.ylabel('Velocity (m/s)')
    plt.title("Max values")
    #plt.title("Mean absolute values x10-4")
    plt.autoscale(True,'x',True)
    plt.show(block=False)


