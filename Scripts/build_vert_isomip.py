import numpy as np
from netCDF4 import Dataset as nc

def main():
    ncfilename = "vert_isomip.nc"
    ncfile = nc(ncfilename,'w')
    #nk = 25
    #depth_max = 710
    nk = 30
    depth_max = 900
    depth_values = np.linspace(0,depth_max,nk+1)
    depth = depth_values[0:-1]
    depth_edges = np.zeros(nk+1)
    depth_edges[0:-1] = depth_values[0:-1] + (depth_values[1:] - depth_values[0:-1])/2
    depth_edges[-1] = depth_values[-1]

    ncfile.createDimension('st_edges_ocean', len(depth_edges))
    ncfile.createDimension('st_ocean', len(depth))
    st_edges_ocean = ncfile.createVariable('st_edges_ocean' ,'float','st_edges_ocean')
    st_edges_ocean[:] = depth_edges
    st_edges_ocean.long_name = "tcell zstar depth"
    st_edges_ocean.units = "meters"
    st_edges_ocean.cartesian_axis = "Z"
    st_edges_ocean.positive = "down"
    st_ocean        = ncfile.createVariable('st_ocean' ,'float','st_ocean')
    st_ocean[:] = depth
    st_ocean.long_name = "tcell zstar depth"
    st_ocean.units = "meters"
    st_ocean.cartesian_axis = "Z"
    st_ocean.positive = "down"
    ncfile.close()
    print(depth_edges)

