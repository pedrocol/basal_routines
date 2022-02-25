import numpy as np
import argparse
import pyproj
from netCDF4 import Dataset as nc

def get_args():
    ''' 
    read arguments
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-f"     , metavar='file_name'        , help="names of input files" , type=str  , required=True)
    parser.add_argument("-epsg"  , metavar='epsg code'        , help="epsg code"            , type=str  , required=True)
    return parser.parse_args()

def main():
    #args = get_args()

    print('===== open netcdf')
    #cinfile=args.f
    cinfile='/scratch/e14/pc5520/BedMachineAntarctica_2020-07-15_v02_test.nc'
    ncfile = nc(cinfile,'r+')

    print('===== READ Coordinates' )
    x = ncfile.variables['x'][:]
    y = ncfile.variables['y'][:]

    print('===== Get grid' )
    xin ,yin = np.meshgrid(x,y)

    print('===== Get lat/lon' )
    p = pyproj.Proj('epsg:3031')
    lon,lat = p(xin,yin,inverse=True)

    print('===== Print corner data' )
    print(lat[0,0], lon[0,0])
    print(lat[-1,-1], lon[-1,-1])
    print(lat[ 0,-1], lon[ 0,-1])
    print(lat[-1, 0], lon[-1, 0])

    print('===== Compute isf draft' )
    hsurf = ncfile.variables['surface'][:,:]
    hice  = ncfile.variables['thickness'][:,:]
    bed   = ncfile.variables['bed'][:,:]
    msk   = ncfile.variables['mask'][:,:]

    # floating ice is 3
    # ocean is 0
    msk[msk==0]  = -1 # bckup ocean cell
    msk[msk==3]  = -1 # bck up floating ice
    msk[msk>=0]  = 0  # mask everything else
    msk[msk==-1] = 1  # defined floating ice cell and ocean cell as ocean cell
    isfd = (hsurf - hice) * msk + bed * (1 - msk)
    test = (hsurf - hice) - bed

    print('===== Write data')
    vlat  = ncfile.createVariable('lat' , 'float',('y', 'x'))
    vlon  = ncfile.createVariable('lon' , 'float',('y', 'x'))
    visfd = ncfile.createVariable('isfd', 'float',('y', 'x'))
    vtest = ncfile.createVariable('test', 'float',('y', 'x'))
    vmsk  = ncfile.createVariable('msk_oce', 'float',('y', 'x'))
    vlat[:,:]  = lat
    vlon[:,:]  = lon
    visfd[:,:] = isfd
    vmsk[:,:] = msk
    vtest[:,:] = test

    print('===== Close file')
    ncfile.close()
     
#if __name__ == '__main__':
#    main()
