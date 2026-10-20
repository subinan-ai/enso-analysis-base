from netCDF4 import Dataset
import numpy as np
import os, sys
import pathlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Data Open (SST anomaly)
i_path = './data/'
i_name = 'ERA.sst_anom.195001_201912_1x1.nc'
f = Dataset(i_path+i_name,'r')
sst = f.variables['sst'][:,:,:]
f.close()

sst = np.array(sst)
sst = sst.reshape(-1,12,180,360) 
sst = sst[:-1,:,:,:].reshape(-1,12,180,360) # SST anomalies (excluding last year), shape: (69, 12, 73, 144)
sst = np.where(sst==-9.99e+08,np.nan,sst)


# Data Open (DJF Niño 3.4 index)
i_path = './data/'
i_name2 = 'DJF_NINO34_idx_1950_2019_1x1.nc'
f2 = Dataset(i_path+i_name2,'r')
nino = f2.variables['nino34'][:]
f2.close()

nino = np.array(nino)
nino = nino.reshape(-1) #69*12
nino = np.where(nino==-9.99e+08,np.nan,nino)
print("nino shape: ", nino.shape)
# Filter SST data where Niño 3.4 index > 0.5
sst_filtered = np.where(nino[:, None, None, None] < -0.5, sst, np.nan) # Replace missing values with NaN

# Get the corresponding years where Niño 3.4 index > 0.5
valid_indices = np.where(nino < -0.5)[0]  
years_djf = 1950 + np.arange(len(nino))    # 1951–2019

print("DJF Nino3.4 > 0.5 Years, Index")
for idx in valid_indices:
    print(f"Year: {years_djf[idx]}, Nino3.4: {nino[idx]:.3f}")

# Compute the composite mean
comp = np.nanmean(sst_filtered, axis=0)  
mcomp = np.nanmean(comp, axis=0) 
print("comp[90,199]: ",mcomp[90,199])
# Data Save

o_path = './data/'
o_name = 'comp_-0.5_DJF_nino_SST.195001_201912'

mcomp = np.array(mcomp)
mcomp.astype('float32').tofile(o_path+o_name+'.gdat')

ctl = open(o_path+o_name+'.ctl','w')
ctl.write('dset ^'+o_name+'.gdat\n')
ctl.write('undef  -9.99e+08\n')
ctl.write('xdef   360  linear   0.  1\n')
ctl.write('ydef    180  linear -90.  1\n')
ctl.write('zdef     1  linear 1 1\n')
ctl.write('tdef     1  linear jan1950 1mo\n')
ctl.write('vars   1\n')
ctl.write('sst    1   1  variable\n')
ctl.write('ENDVARS\n')
ctl.close()

os.system('cdo -f nc import_binary '+o_path+o_name+'.ctl '+o_path+o_name+'.nc')
os.system('rm -f '+o_path+o_name+'.ctl '+o_path+o_name+'.gdat')


i_path = './data/'
i_name = 'comp_-0.5_DJF_nino_SST.195001_201912.nc'
f = Dataset(i_path+i_name,'r')
f.set_auto_mask(False)
comp_sst = f.variables['sst'][:,:,:]
f.close()
print("comp_sst shape: ", comp_sst.shape) #(1,1,180,360)
comp_sst=np.squeeze(comp_sst)
print("comp_sst shape: ", comp_sst.shape)
comp_sst = np.where(comp_sst==-9.99e+08,np.nan,comp_sst) 

cmap = plt.colormaps['bwr']
clevs = np.arange(-1.0,1.0,0.1)

x, y  = np.meshgrid(np.arange(0,360.0,1.0), np.arange(-89.5, 90.0, 1.0))
fig = plt.figure(figsize = (11,8.5))

ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree(central_longitude=180.0))
grd = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                           linewidth=0.5, linestyle="dashed", color="gray", alpha=0.7, zorder=7)
lon_ticks = np.arange(-180, 181, 60)
lat_ticks = np.arange(-90, 91, 30)
grd.ylocator = mticker.FixedLocator(lat_ticks)
ax.set_xticks(lon_ticks,crs=ccrs.PlateCarree())
ax.set_yticks(lat_ticks,crs=ccrs.PlateCarree())
lon_fmt = LongitudeFormatter(number_format=".0f", degree_symbol=r"$^{\circ}$")
lat_fmt = LatitudeFormatter(number_format=".0f", degree_symbol=r"$^{\circ}$")
ax.xaxis.set_major_formatter(lon_fmt)
ax.yaxis.set_major_formatter(lat_fmt)
plt.setp(ax.get_xticklabels(), fontsize=8)
plt.setp(ax.get_yticklabels(), fontsize=8)

ax.coastlines(linewidth=0.7,color='dimgray',zorder=7)

cax = plt.contourf(x, y, comp_sst, levels=clevs, cmap=cmap, transform=ccrs.PlateCarree(),extend='both', zorder=4) 
  
plt.title('(a) DJF Nino3.4 < -0.5',loc='left',fontsize=10)

cax = plt.axes([0.1, 0.2, 0.8, 0.013])
cbar = plt.colorbar(cax=cax, orientation='horizontal')
cbar.ax.tick_params(labelsize=5,direction='in',length=6,width=0.4,color='black',zorder=6)

plt.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,hspace=0.3,wspace=0.3)
plt.savefig(i_path+i_name+'.png', dpi=300)
plt.close()
