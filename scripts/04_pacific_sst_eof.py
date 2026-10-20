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
sst = sst.reshape(-1,12,180,360) #(70, 12, 180, 360)
sst = np.where(sst==-9.99e+08,np.nan,sst)

dat2=np.array(sst)
dat2=dat2.reshape(-1, 180, 360)
dat2=np.where(dat2==-9.99e+08, np.nan, dat2)
lat = np.arange(-89.5, 90, 1.0)   
lon = np.arange(0.5, 360.0, 1.0)    
coslat=np.cos(np.radians(lat))
lat2=np.full_like(dat2, coslat[np.newaxis,...,np.newaxis])
lat2=np.where(np.isnan(dat2)==True, np.nan, lat2)

lat_idx = np.where((lat >= -20) & (lat <= 20))[0]         #태평양 지역
lon_idx = np.where((lon >= 120) & (lon <= 240))[0]

print(lat_idx) 
print(lon_idx)  
#선택한 위도/경도
lat_box = lat[lat_idx.min():lat_idx.max()+1]
lon_box= lon[lon_idx.min():lon_idx.max()+1]
# 태평양 박스 SST anom
pick_dat = dat2[:, lat_idx.min():lat_idx.max()+1,
                   lon_idx.min():lon_idx.max()+1]
nt, ny, nx = pick_dat.shape
print("pick_dat shape:", pick_dat.shape)  # (time, ny, nx) -> (840, 40, 120)
#위도 가중치 적용
coslat=np.cos(np.radians(lat_box))
w=np.sqrt(coslat)
w2=w[:,np.newaxis] #(ny,1)
print("w2 shape: ", w2.shape) #(40, 1)
X=pick_dat*w2[np.newaxis, :, :] # 태평양 지역에 위도 가중치 적용

X=X.reshape(nt, ny*nx)
valid_mask=~np.all(np.isnan(X), axis=0)
X=X[:,valid_mask]
#시간 평균 제거
X_mean=np.nanmean(X, axis=0)
X=X-X_mean
X=np.nan_to_num(X) #NaN을 0으로 채워서 SVD 가능하도록
#SVD로 EOF/PC/분산 기여율 계산=====================
# X = U S Vt
U, s, Vt=np.linalg.svd(X, full_matrices=False)
#EOF 패턴 복원(전체 격자로 채우기)
eofs=np.full((Vt.shape[0], ny*nx), np.nan)
eofs[:, valid_mask] =Vt
eofs=eofs.reshape(-1, ny, nx)
#PC = U*s
pc=U*s
#분산 기여율
var=(s**2)/np.sum(s**2)
print("variance fraction for 1-3: ", var[:3]*100)
#===================================
#EV1-3 & PC1-3 plot
p = './data/'
n = 'EOF_PC_Var'

fig = plt.figure(figsize=(14, 10))

eofs_scaled = np.empty_like(eofs) 
for k in range(eofs.shape[0]):
    eofs_scaled[k] = eofs[k] * (s[k] / np.sqrt(nt-1))

for k in range(3): 

    ax = fig.add_subplot(3, 2, 2*k+1)
    vmin = np.nanmin(eofs_scaled[k])
    vmax = np.nanmax(eofs_scaled[k])
    vabs = max(abs(vmin), abs(vmax))
    levels = np.linspace(-vabs, vabs, 21)
    cf = ax.contourf(
        lon_box,
        lat_box,
        eofs_scaled[k, :, :],  
        levels=levels,
        cmap='RdBu_r',
        extend='both'
    )

    ax.set_xlim(120, 240)
    ax.set_ylim(-20, 20)

    ax.set_title(f"EV{k+1} (EOF{k+1}) — {var[k]*100:.2f}% Var", fontsize=13)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")

    ax.set_xticks([120, 150, 180, 210, 240])
    ax.set_yticks([-20, -10, 0, 10, 20])

    cbar = plt.colorbar(cf, ax=ax, orientation='vertical', pad=0.02)
    cbar.set_label('SST anomaly (°C)')

    ax2 = fig.add_subplot(3, 2, 2*k+2)
    ax2.plot(pc[:, k], 'k', linewidth=0.8)
    ax2.set_title(f"PC{k+1}", fontsize=13)
    ax2.set_xlabel("Time Index")
    ax2.set_ylabel("Amplitude")
    ax2.grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig(p + n + '.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()



