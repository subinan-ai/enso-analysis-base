from netCDF4 import Dataset
import numpy as np
import os, sys
import pathlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#===============DJF Nino3.4 index(69,)와 DJF u/v(일별)간의 각각의 correlation map============
#5도x5도
i_path='./data/'
i_name1='ERA.sst.195001_201912_5x5.nc'

f3=Dataset(i_path+i_name1, 'r')
dat3=f3.variables['sst'][:]
f3.close()
print("sst shape: ",dat3.shape)

dat3=np.array(dat3)
dat3=dat3.reshape(-1,12, 36, 72)
dat3=np.where(dat3==-9.99e+08, np.nan, dat3)
print(dat3.shape)
clim2=np.mean(dat3, axis=0)
anom2=dat3-clim2

# Data Save
o_path = './data/'
o_name = 'ERA.sst_anom.195001_201912'

var = np.array(anom2)
var.astype('float32').tofile(o_path+o_name+'.gdat')

ctl = open(o_path+o_name+'.ctl','w')
ctl.write('dset ^'+o_name+'.gdat\n')
ctl.write('undef  -9.99e+08\n')
ctl.write('xdef   72  linear   0.  5\n')
ctl.write('ydef    36  linear -87.5  5\n')
ctl.write('zdef     1  linear 1 1\n')
ctl.write('tdef   840  linear jan1950 1mo\n')
ctl.write('vars   1\n')
ctl.write('sst    1   1  variable\n')
ctl.write('ENDVARS\n')
ctl.close()

# gdat, ctl file -> nc file
os.system('cdo -f nc import_binary '+o_path+o_name+'.ctl '+o_path+o_name+'.nc')
os.system('rm -f '+o_path+o_name+'.ctl '+o_path+o_name+'.gdat')

#-----sst DJF nino3.4 index---------------------------------
i_path='./data/'
i_name3='ERA.sst_anom.195001_201912.nc'
f2=Dataset(i_path+i_name3, 'r')
print(f2.variables.keys())
dat2=f2.variables['sst'][:]
print("sst_anom shape: ", dat2.shape)
f2.close()
#SST anomaly를 위도 코사인 가중으로 면적가중 평균 -> (69,)
dat2=np.array(dat2)
dat2=dat2.reshape(-1, 36, 72)
dat2=np.where(dat2==-9.99e+08, np.nan, dat2)
lat=np.arange(-87.5, 90, 5)
lat=np.cos(np.radians(lat))
lat=np.full_like(dat2, lat[np.newaxis,...,np.newaxis])
lat=np.where(np.isnan(dat2)==True, np.nan, lat)
pick_dat=dat2[:, 18:20, 38:49]
pick_lat=lat[:, 18:20, 38:49]
idx = np.nansum(pick_dat*pick_lat, axis=(1,2)) / np.nansum(pick_lat, axis=(1,2))

idx = idx.reshape(-1)      # time(year*month)
idx = idx.reshape(-1,12)   # time(year, month)

nino34 = np.append(idx[:-1],idx[1:],axis=1)  # Concatenate overlapping years to create a (69, 24) array

djf_nino34 = np.mean(nino34[:,11:14],axis=1) # Compute the DJF (December-January-February) mean for each year 
djf_nino34 = djf_nino34.reshape(-1)

# Data Save
o_path = './data/'
o_name = 'DJF_NINO34_idx_1950_2019'

# save
djf_nino34 = np.array(djf_nino34)
djf_nino34.astype('float32').tofile(o_path+o_name+'.gdat')

ctl = open(o_path+o_name+'.ctl','w')
ctl.write('dset ^'+o_name+'.gdat\n')
ctl.write('undef  -9.99e+08\n')
ctl.write('xdef     1  linear   0.  5\n')
ctl.write('ydef     1  linear -87.5  5\n')
ctl.write('zdef     1  linear 1 1\n')
ctl.write('tdef   '+str(len(djf_nino34))+'  linear jan1950 1mo\n')
ctl.write('vars   1\n')
ctl.write('nino34    1   1  variable\n')
ctl.write('ENDVARS\n')
ctl.close()

# gdat, ctl file -> nc file
os.system('cdo -f nc import_binary '+o_path+o_name+'.ctl '+o_path+o_name+'.nc')
os.system('rm -f '+o_path+o_name+'.ctl '+o_path+o_name+'.gdat')

#===========u/v anom===============
#u/v_anom
i_path='./data/'
i_name1='v10m_1993-2018_2.5deg_noleap_5x5.nc'

f=Dataset(i_path+i_name1, 'r')
dat=f.variables['v10m'][:]
f.close()
print("v shape: ",dat.shape) #(9490, 1, 36, 72)
dat=np.array(dat)
dat=dat.reshape(-1,365, 36, 72) #(26, 365, 36, 72)
dat=np.where(dat==-9.99e+08, np.nan, dat)
print(dat.shape)
clim=np.mean(dat, axis=0)
print("clim shape: ", clim.shape) #(365, 36, 72)
anom=dat-clim

std = np.std(anom, axis=0)   # standard deviation 
var1 = anom /std              # nomalize
#var1 = var1.reshape(-1,365,36, 72)
#print("var1 shape: ", var1.shape) #(26, 365, 36, 72)
# Data Save
o_path = './data/'
o_name = 'v10m_anom.199301_201812'

var = np.array(var1)
var.astype('float32').tofile(o_path+o_name+'.gdat')

ctl = open(o_path+o_name+'.ctl','w')
ctl.write('dset ^'+o_name+'.gdat\n')
ctl.write('undef  -9.99e+08\n')
ctl.write('xdef   72  linear   0.  5\n')
ctl.write('ydef    36  linear -87.5  5\n')
ctl.write('zdef     1  linear 1 1\n')
ctl.write('tdef   9490  linear jan1993 1dy\n')
ctl.write('vars   1\n')
ctl.write('v10m    1   1  variable\n')
ctl.write('ENDVARS\n')
ctl.close()

# gdat, ctl file -> nc file
os.system('cdo -f nc import_binary '+o_path+o_name+'.ctl '+o_path+o_name+'.nc')
os.system('rm -f '+o_path+o_name+'.ctl '+o_path+o_name+'.gdat')


#-----DJF u/v index---------------------------------
i_path='./data/'
i_name2='v10m_anom.199301_201812.nc'
f2=Dataset(i_path+i_name2, 'r')
print(f2.variables.keys())
dat2=f2.variables['v10m'][:]  
print("v10m_anom shape: ", dat2.shape) #(9490, 1, 36, 72)
f2.close()

dat2 = np.array(dat2)
dat2=dat2[:,0,:,:] #level 차원 제거-> (9490, 36, 72)
nyear = dat2.shape[0] // 365
dat2 = dat2.reshape(nyear, 365, 36, 72)
dat2 = np.where(dat2==-9.99e+08, np.nan, dat2)

djf_u_2yr = np.append(dat2[:-1], dat2[1:], axis=1)

djf_days = djf_u_2yr[:, 334:424, :, :]
print(djf_days.shape)
djf_flat = np.nanmean(djf_days, axis=1)     # (nyear-1, ncell)
print("djf_flat shape: ", djf_flat)

o_path = './data/'
o_name = 'DJF_NINO34_v10m_idx_1993_2018'

djf_flat.astype('float32').tofile(o_path+o_name+'.gdat')

ctl = open(o_path+o_name+'.ctl','w')
ctl.write('dset ^'+o_name+'.gdat\n')
ctl.write('undef  -9.99e+08\n')
ctl.write('xdef     72  linear   0.  5\n')
ctl.write('ydef     36  linear -87.5  5\n')
ctl.write('zdef     1  linear 1 1\n')
ctl.write('tdef   '+str(len(djf_flat))+'  linear jan1994 1yr\n')
ctl.write('vars   1\n')
ctl.write('v_nino34    1   1  variable\n')
ctl.write('ENDVARS\n')
ctl.close()

os.system('cdo -f nc import_binary '+o_path+o_name+'.ctl '+o_path+o_name+'.nc')
os.system('rm -f '+o_path+o_name+'.ctl '+o_path+o_name+'.gdat')

#-------------correlation 
i_path = './data/'
o_path = './data/'
o_name = 'cor_DJF_NINO34_v10m_idx'
name2 = 'DJF_NINO34_idx_1950_2019.nc'
nino = Dataset(i_path + name2, 'r')
nino34 = nino.variables['nino34'][:]   
nino.close()

nino34 = np.array(nino34).reshape(-1)
nino34 = np.where(nino34 == -9.99e+08, np.nan, nino34) 
ny_djf = nino34.shape[0]  # 69

name = 'DJF_NINO34_v10m_idx_1993_2018.nc'
ds_u = Dataset(i_path + name, 'r')
print("variables in u file:", ds_u.variables.keys())

u = ds_u.variables['v_nino34'][:]   
print(u.shape) #(25, 1, 36, 72)
lat = ds_u.variables['lat'][:]
lon = ds_u.variables['lon'][:]
ds_u.close()

u = np.array(u)
u = u[:, 0, :, :]
u = np.where(u == -9.99e+08, np.nan, u)

# 기간 맞추기
base_year = 1951
start_year = 1994
end_year   = 2018
ny_season  =  25

start_idx = start_year - base_year       # 1994-1951 = 43
end_idx   = start_idx + ny_season        # 43+25 = 68 → index 43~67 = 1994~2018

x = nino34[start_idx:end_idx]            # (25,)

print("Nino3.4 shape:", x.shape) #(25,)
print("DJF u10m shape:", u.shape)  # (25, 36, 72)
                  # (36,72)

X = x[:, np.newaxis, np.newaxis]         # (25, 1, 1)
XY = X * u                       # (25, 36, 72)
xy = np.nanmean(XY, axis=0)              # (36, 72)
xx = np.nanmean(X**2)            
yy = np.nanmean(u**2, axis=0)   
cor_u = xy / np.sqrt(xx * yy)            # (36, 72)

o_name_u = 'cor_DJF_nino_v10m_1993-2018_2.5deg_noleap_5x5'

cor_u32 = np.array(cor_u, dtype='float32')
cor_u32.tofile(o_path + o_name_u + '.gdat')

ctl = open(o_path+o_name_u+'.ctl','w')
ctl.write('dset ^'+o_name_u+'.gdat\n')
ctl.write('undef  -9.99e+08\n')
ctl.write('xdef     72  linear   0.  5\n')
ctl.write('ydef     36  linear -87.5  5\n')
ctl.write('zdef     1  linear 1 1\n')
ctl.write('tdef   1  linear jan1950 1mo\n')
ctl.write('vars   1\n')
ctl.write('v    1   1  variable\n')
ctl.write('ENDVARS\n')
ctl.close()

os.system('cdo -f nc import_binary '+o_path+o_name_u+'.ctl '+o_path+o_name_u+'.nc')
os.system('rm -f '+o_path+o_name_u+'.ctl '+o_path+o_name_u+'.gdat')
#########correlation map

p = './data/'
n = 'cor_DJF_nino_v10m_1993-2018_2.5deg_noleap_5x5.nc'
p2 = './data/'
f = Dataset(p2+n,'r')
f.set_auto_mask(False)

nino34_v = f.variables['v'][:]
print(nino34_v.shape) #(1,1,36,72) -> 이전에 저장할 때 표준에 따라 (time, level, lat, lon)으로 저장
nino34_v = np.squeeze(nino34_v)
print("cor shape: ", nino34_v.shape) #(36,72)
nino34_v = np.where(nino34_v==-9.99e+08,np.nan,nino34_v)       

# Plot
cmap = plt.colormaps['bwr']
clevs = np.arange(-1.0, 1.01, 0.25) 
x, y  = np.meshgrid(
    np.arange(0, 360.0, 5.0),     
    np.arange(-87.5, 90.0, 5.0)   
)

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
cf = ax.contourf(x, y, nino34_v, levels=clevs, cmap=cmap, transform=ccrs.PlateCarree(),extend='both', zorder=4) 


plt.title('(a) COR (v & nino3.4)',loc='left',fontsize=10)


cax_cbar = fig.add_axes([0.1, 0.2, 0.8, 0.013]) 
cbar = plt.colorbar(cf, cax=cax_cbar, orientation='horizontal') 
cbar.ax.tick_params(labelsize=5,direction='in',length=6,width=0.4,color='black',zorder=6)
cbar.set_label('Correlation') 

plt.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,hspace=0.3,wspace=0.3)
plt.savefig(p+n+'.png', dpi=300)
plt.close()










