from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.signal import detrend

i_path = './data/'
i_name = 'ERA.sst_anom.195001_201912_1x1.nc'

f = Dataset(i_path + i_name, 'r')
sst = f['sst'][:]
f.close()

sst = np.where(sst == -9.99e+08, np.nan, sst)
sst = sst.reshape(-1, 12, 180, 360)

sst_76_97 = sst[48:]          # 1976–1997 (sst[26:48]) -> 22 years / 1998-2019 (sst[49:]) -> 21 years
ond = np.nanmean(sst_76_97[:, 9:12], axis=1)   # OND 평균 (22, lat, lon)
lat = np.arange(-89.5, 90, 1)
lon = np.arange(0.5, 360, 1)

lat_idx = np.where((lat >= -35) & (lat <=20))[0]
lon_idx = np.where((lon >= 30) & (lon <= 120))[0]

lat_box = lat[lat_idx]
lon_box = lon[lon_idx]

dat = ond[:, lat_idx][:, :, lon_idx]    # (22, ny, nx)
nt, ny, nx = dat.shape
w = np.sqrt(np.cos(np.radians(lat_box)))
dat_w = dat * w[None, :, None]

X = dat_w.reshape(nt, ny*nx)
valid_mask = ~np.any(np.isnan(X), axis=0)
X = X[:, valid_mask]                 

X = detrend(X, axis=0, type='linear')
X = X - X.mean(axis=0)

U, s, Vt = np.linalg.svd(X, full_matrices=False)
pc = U * s

eofs = np.full((Vt.shape[0], ny*nx), np.nan)
eofs[:, valid_mask] = Vt
eofs = eofs.reshape(-1, ny, nx)
print("eofs shape: ", eofs.shape) #(22, 55, 90)
var = (s**2) / np.sum(s**2)

eofs_scaled = np.empty_like(eofs)
for k in range(3):
    eofs_scaled[k] = eofs[k] * s[k] / np.sqrt(nt - 1)


#Plot======================================
p = './data/'
n = 'EOF_PC_Var_1998_2019_FINAL'
fig = plt.figure(figsize=(14, 11))
years = np.arange(1998, 2020)
plot_levels = np.linspace(-0.9, 0.9, 13)

lon_min, lon_max = lon_box.min(), lon_box.max()
lat_min, lat_max = lat_box.min(), lat_box.max()
for k in range(3):
    
    # EOF
    ax = fig.add_subplot(3, 2, 2*k+1, projection=ccrs.PlateCarree())
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m')
    
    cf = ax.contourf(lon_box, lat_box, eofs_scaled[k],
                     levels=plot_levels, cmap='RdBu_r',
                     extend='both', transform=ccrs.PlateCarree())
    
    gl = ax.gridlines(draw_labels=True, linestyle=':', alpha=0.4)
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    
    ax.set_title(f"EOF{k+1} ({var[k]*100:.2f}%)", fontweight='bold')

    # PC
    ax2 = fig.add_subplot(3, 2, 2*k+2)
    ax2.plot(years, pc[:, k], 'r.-')
    ax2.axhline(0, color='black', linewidth=0.8)
    ax2.grid(True, linestyle=':', alpha=0.5)
    ax2.set_title(f"PC{k+1}", fontweight='bold')
    if k != 2:
        ax2.set_xticklabels([])

# 컬러바
cbar_ax = fig.add_axes([0.31, 0.05, 0.33, 0.02])
plt.colorbar(cf, cax=cbar_ax, orientation='horizontal',
             label='SST Anomaly (K)')

plt.tight_layout(rect=[0, 0.09, 1, 1])
plt.savefig(p + n + '.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()






