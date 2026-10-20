from netCDF4 import Dataset
import numpy as np
import os, sys
import pathlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#=============DJF 강수와 DJF Nino3.4 index correlation=================
#1도x1도
#mtpr_anom
i_path='./data/'
i_name1='ERA.sst.195001_201912_1x1.nc'
i_name2='ERA.mtpr.195001_201912_1x1.nc'

#강수 아노말리
f=Dataset(i_path+i_name2,'r')
dat=f.variables['mtpr'][:]
f.close()
print("mtpr shape: ", dat.shape)

dat=np.array(dat)
dat=dat.reshape(-1,12, 180,360)
dat=np.where(dat==-9.99e+08, np.nan, dat)
print(dat.shape)
clim=np.mean(dat, axis=0)
anom=dat-clim
std = np.std(anom, axis=0)   # standard deviation 

var1 = anom /std              # nomalize
var1 = var1.reshape(-1,12,180, 360)

# Data Save
o_path = './data/'
o_name = 'ERA.mtpr_anom.195001_201912'

var = np.array(var1)
var.astype('float32').tofile(o_path+o_name+'.gdat')

ctl = open(o_path+o_name+'.ctl','w')
ctl.write('dset ^'+o_name+'.gdat\n')
ctl.write('undef  -9.99e+08\n')
ctl.write('xdef   360  linear   0.  1\n')
ctl.write('ydef    180  linear -90.  1\n')
ctl.write('zdef     1  linear 1 1\n')
ctl.write('tdef   840  linear jan1950 1mo\n')
ctl.write('vars   1\n')
ctl.write('mtpr    1   1  variable\n')
ctl.write('ENDVARS\n')
ctl.close()

# gdat, ctl file -> nc file
os.system('cdo -f nc import_binary '+o_path+o_name+'.ctl '+o_path+o_name+'.nc')
os.system('rm -f '+o_path+o_name+'.ctl '+o_path+o_name+'.gdat')

#---sst_anom--------------------------
i_path='./data/'
i_name1='ERA.sst.195001_201912_1x1.nc'
f3=Dataset(i_path+i_name1, 'r')
dat3=f3.variables['sst'][:]
f3.close()
print("sst shape: ",dat3.shape)

dat3=np.array(dat3)
dat3=dat3.reshape(-1,12, 180, 360)
dat3=np.where(dat3==-9.99e+08, np.nan, dat3)
print(dat3.shape)
clim2=np.mean(dat3, axis=0)
anom2=dat3-clim2
anom2_flat=anom2.reshape(-1, 180, 360)

#print(f"전체 Anomaly 값의 최소값: {np.nanmin(anom2_flat):.4f} °C")
#print(f"전체 Anomaly 값의 최대값: {np.nanmax(anom2_flat):.4f} °C")
#print(f"전체 Anomaly 값의 평균: {np.nanmean(anom2_flat):.4f} °C (거의 0이어야 함)")
#print("-" * 40)

std = np.std(anom2, axis=0)   # standard deviation 

var = anom2 /std              # nomalize
var = var.reshape(-1,12,180, 360)

# Data Save
o_path = './data/'
o_name = 'ERA.sst_anom.195001_201912_1x1'

var = np.array(var)
var.astype('float32').tofile(o_path+o_name+'.gdat')

ctl = open(o_path+o_name+'.ctl','w')
ctl.write('dset ^'+o_name+'.gdat\n')
ctl.write('undef  -9.99e+08\n')
ctl.write('xdef   360  linear   0.  1\n')
ctl.write('ydef    180  linear -90.  1\n')
ctl.write('zdef     1  linear 1 1\n')
ctl.write('tdef   840  linear jan1950 1mo\n')
ctl.write('vars   1\n')
ctl.write('sst    1   1  variable\n')
ctl.write('ENDVARS\n')
ctl.close()

# gdat, ctl file -> nc file
os.system('cdo -f nc import_binary '+o_path+o_name+'.ctl '+o_path+o_name+'.nc')
os.system('rm -f '+o_path+o_name+'.ctl '+o_path+o_name+'.gdat')

#-----DJF nino3.4 index---------------------------------
i_path='./data/'
i_name1='ERA.sst.195001_201912_1x1.nc'
i_name3='ERA.sst_anom.195001_201912_1x1.nc'
f2=Dataset(i_path+i_name3, 'r')
print(f2.variables.keys())
dat2=f2.variables['sst'][:]
print("sst_anom shape: ", dat2.shape)
f2.close()
#SST anomaly를 위도 코사인 가중으로 면적가중 평균 -> (69,)
dat2=np.array(dat2)
dat2=dat2.reshape(-1, 180, 360)
dat2=np.where(dat2==-9.99e+08, np.nan, dat2)
'''
lat = np.arange(-90, 91, 1)      # (37,)
lon = np.arange(0, 360, 1)       # (72,)
lat_idx = np.where((lat >= -5) & (lat <= 5))[0]         #Nino3.4 지역
lon_idx = np.where((lon >= 190) & (lon <= 240))[0]
print(lat_idx)  # 예: [85 86 ... 95]
print(lon_idx)  # 예: [190 191 ... 240]
pick_dat = dat[:, lat_idx.min():lat_idx.max()+1,
lon_idx.min():lon_idx.max()+1]
pick_lat = lat_3d[:, lat_idx.min():lat_idx.max()+1,
lon_idx.min():lon_idx.max()+1]
'''
lat=np.arange(-89.5, 90, 1)
lat=np.cos(np.radians(lat))
lat=np.full_like(dat2, lat[np.newaxis,...,np.newaxis])
lat=np.where(np.isnan(dat2)==True, np.nan, lat)
pick_dat=dat2[:, 85:96, 190:241]
pick_lat=lat[:, 85:96, 190:241]
idx = np.nansum(pick_dat*pick_lat, axis=(1,2)) / np.nansum(pick_lat, axis=(1,2))

idx = idx.reshape(-1)      # time(year*month)
idx = idx.reshape(-1,12)   # time(year, month)

nino34 = np.append(idx[:-1],idx[1:],axis=1)  # Concatenate overlapping years to create a (69, 24) array

djf_nino34 = np.mean(nino34[:,11:14],axis=1) # Compute the DJF (December-January-February) mean for each year 
djf_nino34 = djf_nino34.reshape(-1)

# Data Save
o_path = './data/'
o_name = 'DJF_NINO34_idx_1950_2019_1x1'

# save
djf_nino34 = np.array(djf_nino34)
djf_nino34.astype('float32').tofile(o_path+o_name+'.gdat')

ctl = open(o_path+o_name+'.ctl','w')
ctl.write('dset ^'+o_name+'.gdat\n')
ctl.write('undef  -9.99e+08\n')
ctl.write('xdef     1  linear   0.  1\n')
ctl.write('ydef     1  linear -90.  1\n')
ctl.write('zdef     1  linear 1 1\n')
ctl.write('tdef   '+str(len(djf_nino34))+'  linear jan1950 1mo\n')
ctl.write('vars   1\n')
ctl.write('nino34    1   1  variable\n')
ctl.write('ENDVARS\n')
ctl.close()

# gdat, ctl file -> nc file
os.system('cdo -f nc import_binary '+o_path+o_name+'.ctl '+o_path+o_name+'.nc')
os.system('rm -f '+o_path+o_name+'.ctl '+o_path+o_name+'.gdat')

#-------------correlation 
i_path = './data/'
name = 'ERA.mtpr_anom.195001_201912.nc'
mtpr_anom = Dataset(i_path+name,'r')

mtpr_anom = mtpr_anom.variables['mtpr'][:,:,:]
print("mtpr anom: " ,mtpr_anom.shape)

mtpr_anom = np.array(mtpr_anom)
mtpr_anom = mtpr_anom.reshape(-1,12,180, 360) 
mtpr_anom = np.where(mtpr_anom==-9.99e+08,np.nan,mtpr_anom)
dec=mtpr_anom[:-1,11,:,:]
jan=mtpr_anom[1:,0,:,:]
feb=mtpr_anom[1:,1,:,:]
djf_mtpr = np.nanmean(np.stack([dec,jan,feb], axis=1), axis=1) 
'''
dec/jan/feb stack -> (69, 3, 180, 360)
nanmean axis=1 -> (69, 180, 360)
'''

name2 = 'DJF_NINO34_idx_1950_2019.nc'
nino = Dataset(i_path+name2,'r')
nino34 = nino.variables['nino34'][:]
nino34 = np.array(nino34)
nino34 = nino34.reshape(-1) # 69
nino34 = np.where(nino34==-9.99e+08,np.nan,nino34)
x=nino34
y=djf_mtpr
x=x.reshape(69)
y=y.reshape(-1,180, 360)

xy=np.mean(x[:,None,None]*y, axis=0)
xx=np.mean(x**2, axis=0)
yy=np.mean(y**2, axis=0)
cor=xy/np.sqrt(xx*yy)

# save
o_path = './data/'
o_name = 'cor_DJF_nino_mtpr.195001_201912'


cor = np.array(cor)
cor.astype('float32').tofile(o_path+o_name+'.gdat')

ctl = open(o_path+o_name+'.ctl','w')
ctl.write('dset ^'+o_name+'.gdat\n')
ctl.write('undef  -9.99e+08\n')
ctl.write('xdef   360  linear   0.  1\n')
ctl.write('ydef    180  linear -90.  1\n')
ctl.write('zdef     1  linear 1 1\n')
ctl.write('tdef     1  linear jan1950 1mo\n')
ctl.write('vars   1\n')
ctl.write('mtpr    1   1  variable\n')
ctl.write('ENDVARS\n')
ctl.close()

os.system('cdo -f nc import_binary '+o_path+o_name+'.ctl '+o_path+o_name+'.nc')
os.system('rm -f '+o_path+o_name+'.ctl '+o_path+o_name+'.gdat')

#======correlation map
p = './data/'
n = 'cor_DJF_nino_mtpr.195001_201912.nc'
p2 = './data/'
f = Dataset(p2+n,'r')
f.set_auto_mask(False)

mtpr = f.variables['mtpr'][:]
print(mtpr.shape)
mtpr = mtpr.reshape(-1,180, 360) 
print("cor_mtpr shape: ", mtpr.shape)
mtpr = np.where(mtpr==-9.99e+08,np.nan,mtpr)    
mtpr=np.nanmean(mtpr, axis=0)    #전체 기간 평균 correlation map

# Plot
cmap = plt.colormaps['bwr']
clevs = np.arange(-1.0,1.0,0.3)

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

cax = plt.contourf(x, y, mtpr, levels=clevs, cmap=cmap, transform=ccrs.PlateCarree(),extend='both', zorder=4) 
  
plt.title('(a) COR DJF u_Nino3.4 & DJF PR',loc='left',fontsize=10)

cax = plt.axes([0.1, 0.2, 0.8, 0.013])
cbar = plt.colorbar(cax=cax, orientation='horizontal')
cbar.ax.tick_params(labelsize=5,direction='in',length=6,width=0.4,color='black',zorder=6)

plt.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,hspace=0.3,wspace=0.3)
plt.savefig(p+n+'.png', dpi=300)
plt.close()
