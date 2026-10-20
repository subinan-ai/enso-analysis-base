from netCDF4 import Dataset
import numpy as np
import os, sys
import pathlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#빗금친 곳이 유의미한 곳 / 빗금이 쳐지지 않은 곳은 상관성이 높게 나와도 유의미하지 않은 지역이라고 해석
i_path1 = './data/ERA.sst_anom.195001_201912.nc'    # Nino3.4  
i_path2 = './data/u10m_mon_anom.199301_201812.nc'    # u anom 

f = Dataset(i_path1,'r')
nino34 = f['sst'][:]
f.close()

f = Dataset(i_path2, 'r')
u = f['u'][:]
f.close()

nino34 = np.ma.masked_equal(nino34, -9.99e+08)

nino34 /= np.std(nino34, axis=0)
u /= np.std(u, axis=0)

## area average nino3.4 ##
dat2=np.array(nino34)
dat2=dat2.reshape(-1, 36, 72)
dat2=np.where(dat2==-9.99e+08, np.nan, dat2)
lat=np.arange(-87.5, 90, 5)
lat=np.cos(np.radians(lat))
lat=np.full_like(dat2, lat[np.newaxis,...,np.newaxis])
lat=np.where(np.isnan(dat2)==True, np.nan, lat)
pick_dat=dat2[:, 18:20, 38:49]
pick_lat=lat[:, 18:20, 38:49]
nino34 = np.nansum(pick_dat*pick_lat, axis=(1,2)) / np.nansum(pick_lat, axis=(1,2))

nino34 = nino34.reshape(70,12)
nino34 = np.append(nino34[:-1,:],nino34[1:,:], axis=1).reshape(69,24)

nino34 = nino34[:,11:14]
nino34 = np.nanmean(nino34, axis=1) # 69
djf_nino34 = nino34[43:43+25]       # 1993-2018 year crop
print(djf_nino34)

## u10 anom ## 
u = u.reshape(26,12,36,72)   

u_2yr = np.append(u[:-1, :, :, :],u[1:,  :, :, :],  axis=1)            # (25, 24, 36, 72)

# DJF 
u_djf = np.mean(u_2yr[:, 11:14, :, :], axis=1)    # (25, 36, 72)
# JJA
u_jja = np.mean(u_2yr[:, 5:8, :, :], axis=1)   
# SON
u_son = np.mean(u_2yr[:, 8:11, :, :], axis=1)    
# MAM(+1)
u_mam = np.mean(u_2yr[:, 14:17, :, :], axis=1)   

djf_nino34 = djf_nino34.reshape(25,1,1)
djf_nino34 = np.repeat(djf_nino34, 36,axis=1)
djf_nino34 = np.repeat(djf_nino34, 72, axis=2)

## corr ## 
xy = np.mean(djf_nino34 * u_djf, axis=0)     # (36,72) 
xx = np.mean(djf_nino34**2, axis=0)
yy = np.mean(u_djf**2, axis=0)
cor_u = xy / np.sqrt(xx * yy)
print(np.min(cor_u))
print(np.max(cor_u))

xy_jja = np.mean(djf_nino34 * u_jja, axis=0) 
xx_jja = np.mean(djf_nino34**2, axis=0)
yy_jja = np.mean(u_jja**2, axis=0)
cor_u_jja = xy_jja / np.sqrt(xx_jja * yy_jja)
print(np.min(cor_u_jja))
print(np.max(cor_u_jja))

xy_son = np.mean(djf_nino34 * u_son, axis=0)  
xx_son = np.mean(djf_nino34**2, axis=0)
yy_son = np.mean(u_son**2, axis=0)
cor_u_son = xy_son / np.sqrt(xx_son * yy_son)
print(np.min(cor_u_son))
print(np.max(cor_u_son))

xy_mam = np.mean(djf_nino34 * u_mam, axis=0)    
xx_mam = np.mean(djf_nino34**2, axis=0)
yy_mam = np.mean(u_mam**2, axis=0)
cor_u_mam = xy_mam / np.sqrt(xx_mam * yy_mam)
print(np.min(cor_u_mam))
print(np.max(cor_u_mam))

##regression##
xy = np.mean(djf_nino34*u_djf, axis=0)   # Covariance term E[XY]
xx = np.mean(djf_nino34**2, axis=0)  # Variance term E[X²]
reg_djf = xy / xx

xy_jja = np.mean(djf_nino34*u_jja, axis=0)   
xx_jja = np.mean(djf_nino34**2, axis=0)  
reg_jja = xy_jja / xx_jja

xy_son = np.mean(djf_nino34*u_son, axis=0)  
xx_son = np.mean(djf_nino34**2, axis=0) 
reg_son = xy_son / xx_son

xy_mam = np.mean(djf_nino34*u_mam, axis=0)  
xx_mam = np.mean(djf_nino34**2, axis=0)  
reg_mam = xy_mam / xx_mam
#########correlation map

p = './data/'
n = 'cor_DJF_nino_u10m_1993-2018_2.5deg_noleap_5x5.nc'
'''
name = 'DJF_NINO34_u10m_idx_1993_2018.nc'
p2 = './data/'
f = Dataset(p2+n,'r')
f.set_auto_mask(False)
N=25
t_crit=2.06

djf_u = f.variables['u'][:]
print(djf_u.shape) #(1,1,36,72) 
djf_u = np.squeeze(djf_u)
print("cor shape: ", djf_u.shape) #(36,72)
djf_u = np.where(djf_u==-9.99e+08,np.nan,djf_u)
'''
'''
f_u=Dataset(p2+name, 'r')
f_u.set_auto_mask(False)
f_u=f_u.variables['u_nino34'][:]
f_u = np.where(f_u==-9.99e+08,np.nan,f_u)   
print("f_u shape: ", f_u.shape)

djf_u_std=np.std(f_u, axis=0)

sig_pos=np.zeros(djf_u.shape)
sig_pos=(djf_u * np.sqrt(25))/djf_u_std
print("sig_pos shape: ", sig_pos.shape) #(36, 72)
'''
# === cor mask 계산 ===
N = 25
t_crit = 2.06
tval_djf = cor_u * np.sqrt(N-2) / np.sqrt(1 - cor_u**2)
djf_mask = np.abs(tval_djf) >= t_crit
mask1 = np.squeeze(djf_mask)
#jja
tval_jja = cor_u_jja * np.sqrt(N-2) / np.sqrt(1 - cor_u_jja**2)
jja_mask = np.abs(tval_jja) >= t_crit
mask2 = np.squeeze(jja_mask)
#son
tval_son = cor_u_son * np.sqrt(N-2) / np.sqrt(1 - cor_u_son**2)
son_mask = np.abs(tval_son) >= t_crit
mask3 = np.squeeze(son_mask)
#mam
tval_mam = cor_u_mam * np.sqrt(N-2) / np.sqrt(1 - cor_u_mam**2)
mam_mask = np.abs(tval_mam) >= t_crit
mask4 = np.squeeze(mam_mask)
# === reg mask 계산 ===
N = 25
t_crit = 2.06
tval_djf_reg = reg_djf * np.sqrt(N-2) / np.sqrt(1 - reg_djf**2)
djf_mask_reg = np.abs(tval_djf_reg) >= t_crit
mask1_reg = np.squeeze(djf_mask)
#jja
tval_jja_reg = reg_jja * np.sqrt(N-2) / np.sqrt(1 - reg_jja**2)
jja_mask_reg = np.abs(tval_jja_reg) >= t_crit
mask2_reg = np.squeeze(jja_mask)
#son
tval_son_reg = reg_son * np.sqrt(N-2) / np.sqrt(1 - reg_son**2)
son_mask_reg = np.abs(tval_son_reg) >= t_crit
mask3_reg = np.squeeze(son_mask)
#mam
tval_mam_reg = reg_mam * np.sqrt(N-2) / np.sqrt(1 - reg_mam**2)
mam_mask_reg = np.abs(tval_mam_reg) >= t_crit
mask4_reg = np.squeeze(mam_mask)
# Plot
# --- (1) 유의성 마스크가 이미 있다고 가정 ---
# N = 25
# t_crit = 2.06
# tval_djf = cor_u * np.sqrt(N-2) / np.sqrt(1 - cor_u**2)
# djf_mask = np.abs(tval_djf) >= t_crit
# mask1 = np.squeeze(djf_mask)   # (36,72)

seasons = {
    'DJF' : (cor_u,   mask1),
    'JJA' : (cor_u_jja, mask2),
    'SON' : (cor_u_son, mask3),
    'MAM' : (cor_u_mam, mask4)
}

titles = ['DJF ',
          'JJA ',
          'SON ',
          'MAM ']

cmap = plt.colormaps['bwr']
clevs = np.arange(-1.0, 1.01, 0.25)

x, y = np.meshgrid(
    np.arange(0, 360.0, 5.0),
    np.arange(-87.5, 90.0, 5.0)
)

fig = plt.figure(figsize=(14, 10))

# ----------------------
#   2 × 2 subplot loop
# ----------------------
for i, (season, (cor_field, mask_field)) in enumerate(seasons.items()):

    ax = fig.add_subplot(2, 2, i+1,
                         projection=ccrs.PlateCarree(central_longitude=180.0))

    # Gridlines
    grd = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                       linewidth=0.5, linestyle="dashed",
                       color="gray", alpha=0.7, zorder=7)
    lon_ticks = np.arange(-180, 181, 60)
    lat_ticks = np.arange(-90, 91, 30)
    grd.ylocator = mticker.FixedLocator(lat_ticks)
    ax.set_xticks(lon_ticks,crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks,crs=ccrs.PlateCarree())

    lon_fmt = LongitudeFormatter(number_format=".0f", degree_symbol=r"$^{\circ}$")
    lat_fmt = LatitudeFormatter(number_format=".0f", degree_symbol=r"$^{\circ}$")
    ax.xaxis.set_major_formatter(lon_fmt)
    ax.yaxis.set_major_formatter(lat_fmt)
    plt.setp(ax.get_xticklabels(), fontsize=6)
    plt.setp(ax.get_yticklabels(), fontsize=6)

    # Coastline
    ax.coastlines(linewidth=0.7, color='dimgray', zorder=7)

    # ---- Main correlation field ----
    cf = ax.contourf(
        x, y, cor_field,
        levels=clevs, cmap=cmap,
        transform=ccrs.PlateCarree(),
        extend='both', zorder=4
    )

    # ---- Significance hatch ----
    sig = np.where(mask_field, 1.0, 0.0)
    ax.contourf(
        x, y, sig,
        levels=[0.5, 1.5],
        colors='none',
        hatches=['///'],
        transform=ccrs.PlateCarree(),
        zorder=8
    )

    ax.set_title(titles[i], fontsize=9, loc='left')

# -------------------------------
#     공통 Colorbar (아래쪽)
# -------------------------------
cax = fig.add_axes([0.1, 0.05, 0.8, 0.015])   # left, bottom, width, height
cbar = plt.colorbar(cf, cax=cax, orientation='horizontal')
cbar.ax.tick_params(labelsize=7, direction='in', length=4, width=0.4)
cbar.set_label('Correlation(u10m & DJF NINO3.4)')

plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.12, 
                    wspace=0.15, hspace=0.20)

plt.savefig(p + 'correlation_4panel.png', dpi=300)
plt.close()

seasons = {
    'DJF' : (reg_djf,   mask1_reg),
    'JJA' : (reg_jja, mask2_reg),
    'SON' : (reg_son, mask3_reg),
    'MAM' : (reg_mam, mask4_reg)
}

titles = ['DJF ',
          'JJA ',
          'SON ',
          'MAM ']

cmap = plt.colormaps['bwr']
clevs = np.arange(-1.0, 1.01, 0.25)

x, y = np.meshgrid(
    np.arange(0, 360.0, 5.0),
    np.arange(-87.5, 90.0, 5.0)
)

fig = plt.figure(figsize=(14, 10))

# ----------------------
#   2 × 2 subplot loop
# ----------------------
for i, (season, (cor_field, mask_field)) in enumerate(seasons.items()):

    ax = fig.add_subplot(2, 2, i+1,
                         projection=ccrs.PlateCarree(central_longitude=180.0))

    # Gridlines
    grd = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                       linewidth=0.5, linestyle="dashed",
                       color="gray", alpha=0.7, zorder=7)
    lon_ticks = np.arange(-180, 181, 60)
    lat_ticks = np.arange(-90, 91, 30)
    grd.ylocator = mticker.FixedLocator(lat_ticks)
    ax.set_xticks(lon_ticks,crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks,crs=ccrs.PlateCarree())

    lon_fmt = LongitudeFormatter(number_format=".0f", degree_symbol=r"$^{\circ}$")
    lat_fmt = LatitudeFormatter(number_format=".0f", degree_symbol=r"$^{\circ}$")
    ax.xaxis.set_major_formatter(lon_fmt)
    ax.yaxis.set_major_formatter(lat_fmt)
    plt.setp(ax.get_xticklabels(), fontsize=6)
    plt.setp(ax.get_yticklabels(), fontsize=6)

    # Coastline
    ax.coastlines(linewidth=0.7, color='dimgray', zorder=7)

    # ---- Main correlation field ----
    cf = ax.contourf(
        x, y, cor_field,
        levels=clevs, cmap=cmap,
        transform=ccrs.PlateCarree(),
        extend='both', zorder=4
    )

    # ---- Significance hatch ----
    sig = np.where(mask_field, 1.0, 0.0)
    ax.contourf(
        x, y, sig,
        levels=[0.5, 1.5],
        colors='none',
        hatches=['///'],
        transform=ccrs.PlateCarree(),
        zorder=8
    )

    ax.set_title(titles[i], fontsize=9, loc='left')

# -------------------------------
#     공통 Colorbar (아래쪽)
# -------------------------------
cax = fig.add_axes([0.1, 0.05, 0.8, 0.015])   # left, bottom, width, height
cbar = plt.colorbar(cf, cax=cax, orientation='horizontal')
cbar.ax.tick_params(labelsize=7, direction='in', length=4, width=0.4)
cbar.set_label('Regression(u10m & DJF NINO3.4)')

plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.12, 
                    wspace=0.15, hspace=0.20)

plt.savefig(p + 'regression_4panel.png', dpi=300)
plt.close()
