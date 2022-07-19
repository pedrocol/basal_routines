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
import os

def findij(lon_s=None,lat_s=None):
    #ncfile = '/g/data/ik11/inputs/access-om2/input_08022019/mom_01deg/ocean_hgrid.nc'
    ncfile = '/scratch/e14/pc5520/ISF/AMERY/ocean_hgrid_global_simple.nc'
    fh = Dataset(ncfile, mode='r')
    #Get data
    #lon = fh.variables[('x')][:2000,4000:]
    #lat = fh.variables[('y')][:2000,4000:]
    lon = fh.variables[('x')][:1000,2000:]
    lat = fh.variables[('y')][:1000,2000:]
    fh.close() #Close file
    npi,npj = lon.shape[0],lon.shape[1]
    dista = lon*0
    contx,conty=0,0
    while contx < npi:
       while conty < npj:
          dista[contx,conty] = dist(lon_s,lon[contx,conty],lat_s,lat[contx,conty])
          conty+=1
       contx+=1
       conty=0
       print(contx)
    minv = np.amin(dista)
    locmin = np.where(dista==minv)
    print(locmin[0],locmin[1],lon[locmin[0],locmin[1]],lat[locmin[0],locmin[1]])
    return dista
    #difflon = lon - lon_s
    #difflat = lat - lat_s
    #minv_lon = np.amin(np.abs(difflon))
    #minv_lat = np.amin(np.abs(difflat))
    #loclon = np.where(abs(difflon)==minv_lon)
    #loclat = np.where(abs(difflat)==minv_lat)
    #print(loclon[0],loclon[1], lon[loclon[0],loclon[1]],lon_s)
    #print(loclat[0],loclat[1], lat[loclat[0],loclat[1]],lat_s)

def dist(ddlona, ddlonb, ddlata, ddlatb):    
    dl_conv = np.pi/180  #for degree to radian conversion
    dl_r    = (6378.137+6356.7523)/2.0 #Earth radius km

    dl_latar   = ddlata*dl_conv
    dl_lonar   = ddlona*dl_conv
    dl_ux      = np.cos(dl_lonar)*np.cos(dl_latar)
    dl_uy      = np.sin(dl_lonar)*np.cos(dl_latar)
    dl_uz      = np.sin(dl_latar)
    dl_prevlat = ddlata
    dl_prevlon = ddlona

    dl_latbr = ddlatb*dl_conv
    dl_lonbr = ddlonb*dl_conv
    dl_vx    = np.cos(dl_lonbr)*np.cos(dl_latbr)
    dl_vy    = np.sin(dl_lonbr)*np.cos(dl_latbr)
    dl_vz    = np.sin(dl_latbr)

    dl_pds   = dl_ux*dl_vx + dl_uy*dl_vy + dl_uz*dl_vz

    if (dl_pds >= 1.):
       dista = 0.
    else:
       dista = dl_r*np.arccos(dl_pds)
    return dista

def grid_double_to_simple():
    cinfile = "/scratch/e14/pc5520/ISF/AMERY/ocean_hgrid_test_double.nc"
    #cinfile = "/scratch/e14/pc5520/ISF/AMERY/ocean_hgrid_original_10_2.nc"
    ncfile = Dataset(cinfile,'r+')
    #Get vars
    xi = ncfile.variables['x'][:,:]
    yi  = ncfile.variables['y'][:,:]
    ncfile.close()

    ncfile = Dataset("/scratch/e14/pc5520/ISF/AMERY/ocean_hgrid_global_simple.nc",'w')
    #Adapt var
    #vx  = ncfile.createVariable('x' , 'float',('nyp', 'nxp'))
    #vy  = ncfile.createVariable('y' , 'float',('nyp', 'nxp'))
    xx = xi[1::2,1::2]*0
    yy = yi[1::2,1::2]*0
    xx[:,:]  = xi[1::2,1::2]
    yy[:,:]  = yi[1::2,1::2]
    #xx = xx[:-1,:-1]
    #yy = yy[:-1,:-1]
    nyp,nxp = xx.shape[0],xx.shape[1]
    ncfile.createDimension('nxp', nxp)
    ncfile.createDimension('nyp', nyp)
    x = ncfile.createVariable('x', 'float', ('nyp', 'nxp'))
    y = ncfile.createVariable('y', 'float', ('nyp', 'nxp'))
    x[:,:] = xx[:,:]
    y[:,:] = yy[:,:]
    ncfile.close()

def modify_isf():
    cinfile = "/home/552/pc5520/AMERY/AMERY_10/MOM_bathy.nc"
    ncfile = Dataset(cinfile,'r+')
    #Get vars
    isf = ncfile.variables[('isf_draft')]
    bathy = ncfile.variables[('depth')][:,:]
    isf2 = isf[:,:]
    #We first remove the small cavities
    isf2[:,80:] = 0
    #We then leave some space between the isf and the bathy
    npy,npx = isf2.shape[0],isf2.shape[1]
    for ji in range(0,npx):
        for jj in range(0,npy):
            if isf2[jj,ji] > 0:
               if bathy[jj,ji] - isf2[jj,ji] < 200:
                  isf2[jj,ji] = max(50,bathy[jj,ji]-200)
    isf[:,:] = isf2[:,:]
    ncfile.close()

def modify_ts():
    #Vert coord data
    cinfile = "/home/552/pc5520/AMERY/AMERY_10/vcoord.nc"
    ncfile = Dataset(cinfile,'r+')
    vert = ncfile.variables[('st_edges_ocean')][:]
    ncfile.close()
    #T and S data
    cinfile = "/scratch/e14/pc5520/ISF/AMERY/linear_temp/ocean_daily_z_t0.nc"
    ncfile = Dataset(cinfile,'r+')
    #Get vars
    so = ncfile.variables[('so')]
    thetao = ncfile.variables[('thetao')]
    thetao[:,:,:,:] = thetao[:,:,:,:] * 0 - 1.9
    #npk = so[0,:,0,0].shape[0]-1
    #for kk in range(0,npk):
    #    so[0,kk,:,:] = 33.5 + 3*(vert[kk]/vert.max())

    npk = so[0,:,0,0].shape[0]-1
    npj = so[0,0,:,0].shape[0]-1
    npi = so[0,0,0,:].shape[0]-1
    for ii in range(0,npi):
        for jj in range(0,npj):
            for kk in range(0,npk):
                if np.ma.is_masked(so[0,kk,jj,ii]) == False: 
                   so[0,kk,jj,ii] = 33.5 + 3*(vert[kk]/vert.max()) 
                   #print(so[0,kk,jj,ii])
        print(ii)
    ncfile.close()

def extendinit():
    ncfile="/scratch/e14/pc5520/ISF/AMERY/scripts/init_TS.nc"
    ncfile = Dataset(ncfile,'r+')
    temp = ncfile.variables['temp']
    salt = ncfile.variables['salt']
    temp2 = ncfile.variables['temp'][:,:,:]
    salt2 = ncfile.variables['salt'][:,:,:]
    npz,npy,npx = temp2.shape[0],temp2.shape[1],temp2.shape[2]
    vertical = "no"
    horizontal = "yes"
    if vertical == "yes":
       for ji in range(0,npx):
           for jj in range(0,npy):
               for jk in range(1,npz):
                   if  np.ma.is_masked(temp2[jk,jj,ji]):
                       temp2[jk,jj,ji] = temp2[jk-1,jj,ji]
                       salt2[jk,jj,ji] = salt2[jk-1,jj,ji]
    if horizontal == "yes":
       for ji in range(0,npx):
           for jj in range(npy-2,-1,-1):
               for jk in range(0,npz):
                   if  np.ma.is_masked(temp2[jk,jj,ji]):
                       temp2[jk,jj,ji] = temp2[jk,jj+1,ji]
                       salt2[jk,jj,ji] = salt2[jk,jj+1,ji]

    temp[:,:,:] = temp2[:,:,:]
    salt[:,:,:] = salt2[:,:,:]
    ncfile.close()


