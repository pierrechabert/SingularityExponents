###### TO RUN IN THE TERMINAL
###### $ f2py -c -m loop_SE_fortran loop_SE_fortran.f90
###### $ ipython
###### In [1]: %run make_SE_fortran_EastMed_2RUN.py
###### pierre.chabertpc@gmail.com
###### https://oceancolor.gsfc.nasa.gov/cgi/browse.pl

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
import warnings
import time
from numba import jit
import loop_SE_fortran
import numpy.ma as ma
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

class SingularityExponents(object):
    """
    """
    def __init__(self, arr):
        """
        arr : values on a 2d grid representing an image or a physical dataset
        """
        self.arr = arr
        assert len(arr.shape) == 2



    @jit
    def _compute_SE(self, arr):
        """
        Compute the Singularity Exponents of the image - array
        """
        sstf = np.asfarray(arr)
        lst5 = np.empty(np.shape(sstf), order='F') * np.nan
        sstf[np.isnan(sstf)] = 0
        loop_SE_fortran.sh_treat(sstf, lst5)
        return(lst5)



    @jit
    def _export_SE(self, arr_exp, date='Not Communicated', filename='SE.nc',
                   lonmin=None, lonmax=None, latmin=None, latmax=None):
        dataset = Dataset(filename , 'w')
        dataset.description = 'Singularity Exponents computation on the %s L2 Aqua MODIS SST data (lonmin = %s, lonmax = %s, latmin = %s, latmax = %s' % (date, lonmin, lonmax, latmin, latmax)

        dataset.history = 'Created ' + time.ctime(time.time())
        dataset.source = 'IMEDEA (CSIC-UIB)'
        lat = dataset.createDimension('lat', lst5[0,:].size)
        lon = dataset.createDimension('lon', lst5[:,0].size)

        latitudes = dataset.createVariable('latitude', np.float32, ('lat',))
        longitudes = dataset.createVariable('longitude', np.float32, ('lon',))
        sedata = dataset.createVariable('se', np.float32, ('lat','lon'))

        latitudes[:] = np.arange(0,
                                 lst5[0,:].size,
                                 lst5[0,:].size / lst5[0,:].size)
        longitudes[:] = np.arange(0,
                                  lst5[:,0].size,
                                  lst5[:,0].size / lst5[:,0].size)
        sedata[:] = lst5

        latitudes.units = 'degree_north'
        longitudes.units = 'degree_east'
        se.units = 'dimensionless'

        dataset.close()




    def _ref_image(self, lon1, lat1, sst1, dte, tmin, tmax, M):
        '''
        creating reference figure, with raw sst values (in black and white)
        '''
        plt.figure()
        M.pcolormesh(lon1, lat1, sst1, cmap='Greys',vmin=tmin, vmax=tmax)
        plt.axis('off')
        plt.savefig('%s_plt_SST_ref.png'%dte,bbox_inches='tight', dpi=300)
        




    def _se_image(self, lon1, lat1, lst5, dte, M):
        '''
        creating computed SE figure
        '''
        plt.figure()
        M.pcolormesh(lon1, lat1, lst5, cmap='Greys',vmin=-.4, vmax=.1)
        plt.colorbar()
        M.fillcontinents()
        plt.axis('off')
        plt.savefig('%s_plt_SE.png'%dte,bbox_inches='tight', dpi=300)




    def _open_se_netcdf(self, filename):
        '''
        open se netCDF file previously created
        '''
        with Dataset(filename) as nc:
            lon = nc.variables['longitude'][:]
            lat = nc.variables['latitude'][:]
            se = nc.variables['se'][:]

        return(lon, lat, se)


    def _detect_edge(self, se_arr):
        '''
        detect minimum se in great number to enlight edges,
        return se_enlightenned, play with semin/ semax to enlight frontal areas
        '''
        se0 = se_arr
        nb = 0
        semin, semax = -.4, -.2
        for i in range(1, se0.shape[0] - 1):
            for j in range(1, se0.shape[1] - 1):
                if se0[i, j] > semin and se0[i, j] < semax:
                    se0[i, j] = -10000
                if np.isnan(se0[i, j]) == True:
                    se0[i, j] = 1000
                if ma.is_masked(se0[i, j]) == True:
                    se0[i, j] = 1000

        se_en = ma.masked_where(se0 > -1000, se0)

        return(se_en)

    def _plot_edges_detected(self, lon1, lat1, lst5, dte, se_en, M):
        '''
        create an image of enlighted edges
        '''
        plt.figure()
        M.pcolormesh(lon1, lat1, lst5, cmap='Greys',vmin=-.4, vmax=.1)
        plt.colorbar()
        M.fillcontinents()
        plt.axis('off')
        M.pcolormesh(lon1, lat1, se_en,
                     vmin=-10000, vmax=0, cmap='PuRd_r')

        plt.savefig('%s_plt_SE_en.png'%dte,bbox_inches='tight', dpi=300)
        



if __name__ == '__main__':

    t0 = time.time()

    warnings.filterwarnings("ignore", category=RuntimeWarning)
    # ------------------------------------------------------------------------------------

    #mur_file = 'A2017201014500.L2_LAC_SST4.nc'
    mur_file = "A2017206020500.L2_LAC_SST4.nc"

    dte = mur_file[1:8]

    # if night data: change sst -> sst4 in nc.variables
    night_sst4 = 'no'
    night_sst4 = 'yes'

    # ------------------------------------------------------------------------------------

    # parameters (coordinates limits, values of sst max/min) of the plot
    lonmin, lonmax, latmin, latmax = -3, 0, 35, 38

    M = Basemap(projection = 'cyl',llcrnrlon = lonmin,
                                   urcrnrlon = lonmax,
                                   llcrnrlat = latmin,
                                   urcrnrlat = latmax,
                                   lat_ts = 37.,
                                   resolution = 'h')


    # netCDF file processing

    with Dataset(mur_file) as nc:

        nav = nc.groups['navigation_data']
        lon = nav.variables['longitude'][:]
        lat = nav.variables['latitude'][:]
        geo = nc.groups['geophysical_data']

        if night_sst4 == 'no':
            sst = geo.variables['sst'][:]
        if night_sst4 == 'yes':
            sst = geo.variables['sst4'][:]

    x, y = lon, lat

    # --------------------------------------------------------------------------
    # data processing

    # take the data in the limits defined previously
    lon0 = lon[(lon > lonmin) & (lon < lonmax) & (lat > latmin) & (lat < latmax)]
    lat0 = lat[(lon > lonmin) & (lon < lonmax) & (lat > latmin) & (lat < latmax)]
    sst0 = sst[(lon > lonmin) & (lon < lonmax) & (lat > latmin) & (lat < latmax)]



    # creating a grid with all the values (interpolation made)
    m = 500
    lon1,lat1 = np.meshgrid(np.linspace(np.min(lon0),np.max(lon0),m),
                            np.linspace(np.min(lat0),np.max(lat0),m))
    sst1 = griddata((lon0, lat0), sst0,(lon1, lat1))

    tmin, tmax = 22, 28

    #--------------------------------------------------------------------------
    # REAL PROCESSING (calculating the Singularity Exponents)
    # calling the function sh_treat (fortran)


    # call SE module
    se = SingularityExponents(sst1)


    # creating reference figure
    se._ref_image(lon1, lat1, sst1, dte, tmin, tmax, M)

    # compute SE from SST
    lst5 = se._compute_SE(sst1)

    # create SE image
    se._se_image(lon1, lat1, lst5, dte, M)

    # export SE in a netCDF file
    fname = 'SE_%s.nc'%dte
    se._export_SE(lst5, date=dte, filename=fname,
                  lonmin=lonmin , lonmax=lonmax ,
                  latmin=latmin, latmax=latmax)

    # open previous se netCDF file
    lon_exp, lat_exp, se_exp = se._open_se_netcdf(filename=fname)

    # coordinates have to be recomputed
    lon_exp0 = np.arange(lonmin, lonmax, abs(lonmax - lonmin)/se_exp[0,:].size)
    lat_exp0 = np.arange(latmin, latmax, abs(latmax - latmin)/se_exp[:,0].size)

    plt.figure()
    M.pcolormesh(lon_exp0, lat_exp0, se_exp,
                 cmap='Greys', vmin=-.4, vmax=.1)
    M.fillcontinents()
    plt.savefig('%s_plt_SE_nc.png'%dte,bbox_inches='tight', dpi=300)
    plt.close()

    # detect edges
    se_en = se._detect_edge(lst5)

    # plot edges detected
    se._plot_edges_detected(lon1, lat1, lst5, dte, se_en, M)
    #--------------------------------------------------------------------------
    #print the processing duration

    print('Time after processing:', round(time.time()-t0,2))

