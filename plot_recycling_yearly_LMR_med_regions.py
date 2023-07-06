############################
###Import Python packages###
############################
from netCDF4 import Dataset
from numpy import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import sys
from scipy.interpolate import RegularGridInterpolator
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

#########################
###Plotting parameters###
#########################
savepath = 'C:/Users/jthe/Documents/Python/project_ob_1/figs_mediterranean/'

mpl.rcParams['figure.figsize'] = [20.0, 8.0]
mpl.rcParams['figure.dpi'] = 80
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['font.size'] = 12
mpl.rcParams['legend.fontsize'] = 12
mpl.rcParams['figure.titlesize'] = 12

########################
###Coordinate to grid###
########################
def coordinate_to_grid(lat,lon):
    grid_lat = 180 + (-2*lat)
    grid_lon = 360 + (2*lon)
    return grid_lat, grid_lon

#################
###Import data###
#################
def obtain_rec(month):
    f_lsm = Dataset('C:/Users/jthe/Documents/Python/project_ob_1/netCDF/shape_files/lsm_nan.nc')
    lsm = f_lsm.variables['lsm'][:,:]
    lsm[lsm==0.0] = np.nan

    f=Dataset('C:/Users/jthe/Documents/Python/project_ob_1/netCDF/metrics/group_recycling/9_grids/moisture_recycling_in_group_'+str(month).zfill(2)+'.nc')
    rec = f.variables['recycling'][:,:]
    rec = np.multiply(rec,lsm)
    plotrec=np.zeros(rec.shape)
    plotrec[:,:360]=rec[:,360:]
    plotrec[:,360:]=rec[:,:360]
    return plotrec

rec_months = np.zeros((12,360,720))
for i in np.arange(12):
    month=i+1
    rec_months[i,:,:] = obtain_rec(month)

ya_rec = np.nanmean(rec_months, axis=0)
print(ya_rec.shape)

##############
###Plotting###
##############
lats=np.arange(90,-90,-0.5)
lons=np.arange(-180,180,0.5)
#australia
extent_au = [106, 154, -48, -15]
central_lon_au, central_lat_au = 130, -31.5
#california
extent_ca = [-131, -105, 20, 52]
central_lon_ca, central_lat_ca = -118, 36
#Chile
central_lon_ch, central_lat_ch = -70, -32
extent_ch = [-80, -60, -54, -10]
#Med basin
central_lon_mb, central_lat_mb = 13, 36
extent_mb = [-20, 45, 23, 48]
#S. Africa
central_lon_sa, central_lat_sa = 24, -33
extent_sa = [10, 38, -40, -26]

##Shapefile
fname = r'C:/Users/jthe/Documents/Python/shape_file/shapefile_mediterranean.shp'

def plt_local_recycling(ya_rec,extent,central_lon,central_lat,region):
    cm=mpl.cm.get_cmap('YlGnBu')
    cols=[]
    ticks = [0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06]
    for i in np.arange(0,1,1.0/(len(ticks)+1)):
        cols.append(cm(i))
    cols.append(cm(0.95))
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.Orthographic(central_lon, central_lat))
    ax.set_extent(extent)
    ax.coastlines(resolution='50m')
    cs = ax.contourf(lons,lats,ya_rec,ticks,colors=cols,transform=ccrs.PlateCarree(),extend='max')
    cax,kw = mpl.colorbar.make_axes(ax,location='right',pad=0.05,shrink=0.6)
    out=fig.colorbar(cs,cax=cax,extend='both',**kw)#,ticks = ticks_label)
    out.set_label('Local moisture recycling [-]',size=20)
    out.ax.tick_params(labelsize=18)
    ax.set_title(region, size=20);


    shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                    ccrs.PlateCarree(), edgecolor='black')

    #ax.add_feature(shape_feature,edgecolor='black',facecolor='black',alpha=0.1)
    ax.add_feature(shape_feature,edgecolor='black',facecolor='none', hatch='//')#,alpha=0.1) #black area where locations are.
    #ax.add_feature(shape_feature,edgecolor='black', facecolor='none') #only a black border
    #plt.show()
    plt.savefig('C:/Users/jthe/Documents/Python/project_ob_1/figs_mediterranean/shape_LMR_'+region+'.png')
    plt.savefig('C:/Users/jthe/Documents/Python/project_ob_1/figs_mediterranean/shape_LMR_'+region+'.pdf')
    plt.close()

plt_local_recycling(ya_rec,extent_au,central_lon_au,central_lat_au,'Australia')
plt_local_recycling(ya_rec,extent_ca,central_lon_ca,central_lat_ca,'California')
plt_local_recycling(ya_rec,extent_ch,central_lon_ch,central_lat_ch,'Chile')
plt_local_recycling(ya_rec,extent_mb,central_lon_mb,central_lat_mb,'Mediterranean Basin')
plt_local_recycling(ya_rec,extent_sa,central_lon_sa,central_lat_sa,'South Africa')
