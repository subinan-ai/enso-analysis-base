from netCDF4 import Dataset
import numpy as np
import os, sys
import pathlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

i_path="./data/"
i_name="ERA.sst.195001_201912_1x1.nc"

i_path2="./data/"
i_name2="ERA.mtpr.195001_201912_1x1.nc"

#계절별 sst climatology를 한 장에 그리기
f=Dataset(i_path+i_name, 'r')
sst=f.variables['sst'][:]
f.close()
print("sst shape: ", sst.shape)
sst=np.array(sst)
sst=np.where(sst==-9.99e+08, np.nan, sst)
sst=sst.reshape(-1, 12, 180, 360)   #(70, 12, 180, 360)
print(sst.shape)
#sst=(sst-np.nanmean(sst, axis=0))/np.nanstd(sst, axis=0)    #sst normalization
sst_2yr=np.append(sst[:-1,:,:,:], sst[1:,:,:,:], axis=1)   #(69, 24, 180, 360)
# --- 연도별 계절 평균 (sst) ---
sst_MAM_yr = np.mean(sst[:,2:5,:,:], axis=1)       # (70,180,360)
sst_JJA_yr = np.mean(sst[:,5:8,:,:], axis=1)
sst_SON_yr = np.mean(sst[:,8:11,:,:], axis=1)
sst_DJF_yr = np.mean(sst_2yr[:,11:14,:,:], axis=1) # (69,180,360)

# --- climatology (°C) ---
clim_MAM = np.nanmean(sst_MAM_yr, axis=0) - 273
clim_JJA = np.nanmean(sst_JJA_yr, axis=0) - 273
clim_SON = np.nanmean(sst_SON_yr, axis=0) - 273
clim_DJF = np.nanmean(sst_DJF_yr, axis=0) - 273

'''
print(clim_MAM.shape)
print(clim_JJA.shape)
print(clim_SON.shape)
print(clim_DJF.shape)
'''
#Plot
x, y=np.meshgrid(
    np.arange(0, 360.0, 1.0),
    np.arange(-89.5, 90.0, 1.0)
)
cmap=plt.colormaps['RdBu_r']
all=np.stack([clim_MAM, clim_JJA, clim_SON, clim_DJF], axis=0)
vmin=np.nanmin(all)
vmax=np.nanmax(all)
clevs=np.linspace(vmin, vmax, 20)

fig=plt.figure(figsize=(11, 8.5))
fields=[clim_MAM, clim_JJA, clim_SON, clim_DJF]
titles=['(a) MAM', '(b) JJA', '(c) SON', '(d) DJF']

for i, (dat, ttl) in enumerate(zip(fields, titles), start=1):
    ax=fig.add_subplot(2, 2, i, projection=ccrs.PlateCarree(central_longitude=180.0))
    grd = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                       linewidth=0.5, linestyle="dashed",
                       color="gray", alpha=0.7, zorder=7)
    lon_ticks = np.arange(-180, 181, 60)
    lat_ticks = np.arange(-90, 91, 30)
    grd.ylocator = mticker.FixedLocator(lat_ticks)
    ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())
    lon_fmt = LongitudeFormatter(number_format=".0f", degree_symbol=r"$^{\circ}$")
    lat_fmt = LatitudeFormatter(number_format=".0f", degree_symbol=r"$^{\circ}$")
    ax.xaxis.set_major_formatter(lon_fmt)
    ax.yaxis.set_major_formatter(lat_fmt)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)

    # 해안선
    ax.coastlines(linewidth=0.7, color='dimgray', zorder=7)

    # 필드 플로팅
    cf = ax.contourf(
        x, y, dat,
        levels=clevs,
        cmap=cmap,
        transform=ccrs.PlateCarree(),
        extend='both',
        zorder=4
    )

    # 제목
    ax.set_title(ttl, loc='left', fontsize=10)
#공통 컬러바
cax_cbar=fig.add_axes([0.1, 0.08, 0.8, 0.013])
cbar=plt.colorbar(cf, cax=cax_cbar, orientation='horizontal')
cbar.ax.tick_params(labelsize=7, direction='in', length=6, width=0.4, color='black', zorder=6)
cbar.set_label('SST Climatology (°C)', fontsize=8)

plt.subplots_adjust(
    bottom=0.15, top=0.93, left=0.06, right=0.97, hspace=0.25, wspace=0.1
)
out_name='ERA_sst_seasonal_clim.png'
plt.savefig(i_path+out_name, dpi=300)
plt.close()
 
#계절별 mtpr climatology를 한 장에 그리기
f2=Dataset(i_path2+i_name2, 'r')
mtpr=f2.variables['mtpr'][:]
f2.close()
print("mtpr shape: ", mtpr.shape)
mtpr=np.array(mtpr)
mtpr=np.where(mtpr==-9.99e+08, np.nan, mtpr)
mtpr=mtpr.reshape(-1, 12, 180, 360)

mtpr_2yr=np.append(mtpr[:-1,:,:,:], mtpr[1:,:,:,:], axis=1)
# --- 연도별 계절 평균 (mtpr) ---
mtpr_MAM_yr = np.mean(mtpr[:,2:5,:,:], axis=1) 
mtpr_JJA_yr = np.mean(mtpr[:,5:8,:,:], axis=1)
mtpr_SON_yr = np.mean(mtpr[:,8:11,:,:], axis=1)
mtpr_DJF_yr = np.mean(mtpr_2yr[:,11:14,:,:], axis=1)   # (69,180,360)

# --- climatology (mm/day) ---
clim_mtpr_MAM = np.nanmean(mtpr_MAM_yr, axis=0) * 86400 
clim_mtpr_JJA = np.nanmean(mtpr_JJA_yr, axis=0) * 86400
clim_mtpr_SON = np.nanmean(mtpr_SON_yr, axis=0) * 86400
clim_mtpr_DJF = np.nanmean(mtpr_DJF_yr, axis=0) * 86400

'''
print(clim_mtpr_MAM.shape)
print(clim_mtpr_JJA.shape)
print(clim_mtpr_SON.shape)
print(clim_mtpr_DJF.shape)
'''
#Plot
x, y=np.meshgrid(
    np.arange(0, 360.0, 1.0),
    np.arange(-89.5, 90.0, 1.0)
)
cmap=plt.colormaps['RdBu_r']
all=np.stack([clim_mtpr_MAM, clim_mtpr_JJA, clim_mtpr_SON, clim_mtpr_DJF], axis=0)
vmin=np.nanmin(all)
vmax=np.nanmax(all)
clevs=np.linspace(vmin, vmax, 20)

fig=plt.figure(figsize=(11, 8.5))
fields=[clim_mtpr_MAM, clim_mtpr_JJA, clim_mtpr_SON, clim_mtpr_DJF]
titles=['(a) MAM', '(b) JJA', '(c) SON', '(d) DJF']

for i, (dat, ttl) in enumerate(zip(fields, titles), start=1):
    ax=fig.add_subplot(2, 2, i, projection=ccrs.PlateCarree(central_longitude=180.0))
    grd = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                       linewidth=0.5, linestyle="dashed",
                       color="gray", alpha=0.7, zorder=7)
    lon_ticks = np.arange(-180, 181, 60)
    lat_ticks = np.arange(-90, 91, 30)
    grd.ylocator = mticker.FixedLocator(lat_ticks)
    ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())
    lon_fmt = LongitudeFormatter(number_format=".0f", degree_symbol=r"$^{\circ}$")
    lat_fmt = LatitudeFormatter(number_format=".0f", degree_symbol=r"$^{\circ}$")
    ax.xaxis.set_major_formatter(lon_fmt)
    ax.yaxis.set_major_formatter(lat_fmt)
    plt.setp(ax.get_xticklabels(), fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)

    # 해안선
    ax.coastlines(linewidth=0.7, color='dimgray', zorder=7)

    # 필드 플로팅
    cf = ax.contourf(
        x, y, dat,
        levels=clevs,
        cmap=cmap,
        transform=ccrs.PlateCarree(),
        extend='both',
        zorder=4
    )

    # 제목
    ax.set_title(ttl, loc='left', fontsize=10)
#공통 컬러바
cax_cbar=fig.add_axes([0.1, 0.08, 0.8, 0.013])
cbar=plt.colorbar(cf, cax=cax_cbar, orientation='horizontal')
cbar.ax.tick_params(labelsize=7, direction='in', length=6, width=0.4, color='black', zorder=6)
cbar.set_label('mtpr Climatology (mm/day)', fontsize=8)

plt.subplots_adjust(
    bottom=0.15, top=0.93, left=0.06, right=0.97, hspace=0.25, wspace=0.1
)
out_name='ERA_mtpr_seasonal_clim.png'
plt.savefig(i_path+out_name, dpi=300)
plt.close()


#sst와 mtpr climatology의 계절별 pattern correlation
#MAM
mask_MAM = ~np.isnan(clim_MAM) & ~np.isnan(clim_mtpr_MAM)

A = clim_MAM[mask_MAM].flatten()          # SST
B = clim_mtpr_MAM[mask_MAM].flatten()     # mtpr

A_mean = np.mean(A)
B_mean = np.mean(B)

A_anom = A - A_mean
B_anom = B - B_mean

num = np.sum(A_anom * B_anom)
den = np.sqrt(np.sum(A_anom**2) * np.sum(B_anom**2))

pc_MAM = num / den
print("MAM pattern corr =", pc_MAM)
#JJA
mask_JJA= ~np.isnan(clim_JJA) & ~np.isnan(clim_mtpr_JJA)

A = clim_JJA[mask_JJA].flatten()
B = clim_mtpr_JJA[mask_JJA].flatten()

A_mean = np.mean(A)
B_mean = np.mean(B)

A_anom = A - A_mean
B_anom = B - B_mean

num = np.sum(A_anom * B_anom)
den = np.sqrt(np.sum(A_anom**2) * np.sum(B_anom**2))

pc_JJA = num / den
#SON
mask_SON= ~np.isnan(clim_SON) & ~np.isnan(clim_mtpr_SON)

A = clim_SON[mask_SON].flatten()
B = clim_mtpr_SON[mask_SON].flatten()

A_mean = np.mean(A)
B_mean = np.mean(B)

A_anom = A - A_mean
B_anom = B - B_mean

num = np.sum(A_anom * B_anom)
den = np.sqrt(np.sum(A_anom**2) * np.sum(B_anom**2))

pc_SON = num / den
#DJF
mask_DJF= ~np.isnan(clim_DJF) & ~np.isnan(clim_mtpr_DJF)

A = clim_DJF[mask_DJF].flatten()
B = clim_mtpr_DJF[mask_DJF].flatten()

A_mean = np.mean(A)
B_mean = np.mean(B)

A_anom = A - A_mean
B_anom = B - B_mean

num = np.sum(A_anom * B_anom)
den = np.sqrt(np.sum(A_anom**2) * np.sum(B_anom**2))

pc_DJF = num / den

print("MAM (sst & mtpr): ", pc_MAM)
print("JJA (sst & mtpr): ", pc_JJA)
print("SON (sst & mtpr): ", pc_SON)
print("DJF (sst & mtpr): ", pc_DJF)
