from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import detrend
#저위도/중위도/고위도별 강수 아노말리의 트렌드 제거 전과 후의 시계열
#강수는 트렌드를 가지는 변수가 아니라서 트렌드 제거 전과 후가 비슷하게 나옴
# =========== 데이터 읽기 ===========
i_path = './data/'
i_name = 'ERA.mtpr.195001_201912_1x1.nc'
f = Dataset(i_path + i_name, 'r')
f.set_auto_mask(False)

mtpr = f.variables['mtpr'][:]       
lat  = f.variables['lat'][:]      
lon  = f.variables['lon'][:]     
f.close()
#원본 mtpr 데이터 -> 단위 변환 -> anomaly -> 정규화
print("mtpr raw shape:", mtpr.shape)  # (840, 1, 180, 360)
mtpr = np.squeeze(mtpr)             # (840, 180, 360)
mtpr = np.array(mtpr)
mtpr = np.where(mtpr == -9.99e+08, np.nan, mtpr)
mtpr = mtpr * 86400     #mm/day 
mtpr=mtpr.reshape(-1, 12, 180, 360)
clim=np.nanmean(mtpr, axis=0)
anom=mtpr-clim
std=np.nanstd(anom, axis=0)
anom_std= anom/std 
mtpr=anom_std.reshape(840, 180, 360)

nt, ny, nx = mtpr.shape
lat_weight=np.cos(np.radians(lat))  #위도 가중치 (cos(lat))
#경도 (80-140E)
lon_idx=np.where((lon >= 80) & (lon <= 140))[0]

#(a) 저위도 (80–140E, 0–20N)
lat_idx_a=np.where((lat>=0) & (lat<=20))[0]
#(b) 중위도 (80–140E, 20–40N)
lat_idx_b=np.where((lat>=20) & (lat<=40))[0]
#(c) 고위도 (80-140E, 40-60N)
lat_idx_c=np.where((lat>=40) & (lat<=60))[0]

#(a) 저위도 시계열
a=mtpr[:,lat_idx_a,:][:,:,lon_idx]  #데이터 슬라이싱
w_a=lat_weight[lat_idx_a] #해당 영역의 위도에 해당하는 가중치 -> (nlat,)
w_a=np.repeat(w_a[:,np.newaxis], len(lon_idx), axis=1)  #lon 방향 복제
w_a=w_a[np.newaxis,:,:] #시간축에 대해 가중치는 동일하므로 한 번만 넣고 broadcasting으로 적용 -> 자동으로 (840,nlat,nlon)처럼 적용됨

num_a=np.nansum(a*w_a, axis=(1,2))    #각각의 그리드에 대해 dat(t,lat,lon)*cos(lat) -> 월별 면적가중 총량
num_b=np.nansum(np.where(np.isnan(a), np.nan, w_a ), axis=(1,2))  #결측값 있는 그리드 가중치 제외

mtpr_a=num_a/num_b
mtpr_a_year=mtpr_a.reshape(-1,12).mean(axis=1)
mtpr_a_detrend=detrend(mtpr_a, type='linear')
mtpr_a_year_detrend=detrend(mtpr_a_year, type='linear')
#(b) 중위도 시계열
b=mtpr[:,lat_idx_b,:][:,:,lon_idx]
w_b=lat_weight[lat_idx_b]
w_b=np.repeat(w_b[:,np.newaxis], len(lon_idx), axis=1)
w_b=w_b[np.newaxis,:,:]

num_c=np.nansum(b*w_b, axis=(1,2))
num_d=np.nansum(np.where(np.isnan(b), np.nan, w_b), axis=(1,2))

mtpr_b=num_c/num_d
mtpr_b_year=mtpr_b.reshape(-1, 12).mean(axis=1)
mtpr_b_detrend=detrend(mtpr_b, type='linear')
mtpr_b_year_detrend=detrend(mtpr_b_year, type='linear')
#(c) 고위도 시계열
c=mtpr[:,lat_idx_c,:][:,:,lon_idx]
w_c=lat_weight[lat_idx_c]
w_c=np.repeat(w_c[:,np.newaxis], len(lon_idx), axis=1)
w_c=w_c[np.newaxis,:,:]

num_e=np.nansum(c*w_c, axis=(1,2))
num_f=np.nansum(np.where(np.isnan(c), np.nan, w_c), axis=(1,2))

mtpr_c=num_e/num_f
mtpr_c_year=mtpr_c.reshape(-1,12).mean(axis=1)
mtpr_c_detrend=detrend(mtpr_c, type='linear')
mtpr_c_year_detrend=detrend(mtpr_c_year, type='linear')

#==================
#월별
t=np.arange(nt)
years=nt//12
start_year=1950
years=np.arange(start_year, start_year+years)
tick_years=years[::5]
tick_idx=(tick_years-start_year)*12
#연도별
years_2=np.arange(1950, 1950+len(mtpr_a_year))
t2=years_2

#======Plot==============
fig=plt.figure(figsize=(11,9))
#(a) 저위도
ax1=fig.add_subplot(3, 1, 1)
ax1.plot(t2, mtpr_a_year, 'k', label='Low lat', zorder=5)
ax1.plot(t2, mtpr_a_year_detrend, 'r', label="Low lat Detrended", zorder=6)
ax1.set_ylim(-1,1)
ax1.set_xlim(1950, 2019)    #연도별에서 추가
ax1.set_xticks(np.arange(1950,2020,5))    #ax1.set_xticks(tick_idx) ->월별 / ax1.set_xticks(np.arange(1950,2020,5))
#ax1.set_xticklabels(tick_years)
ax1.set_ylabel('mtpr anomaly', fontsize=11)
ax1.set_title('(a) Low lat (80-140E, 0-20N)', loc='left')
ax1.grid(alpha=0.3)
ax1.legend(loc='upper left', fontsize=9)
#(b) 중위도
ax1=fig.add_subplot(3, 1, 2)
ax1.plot(t2, mtpr_b_year, 'k', label='Mid lat', zorder=5)
ax1.plot(t2, mtpr_b_year_detrend, 'r', label="Mid lat Detrended", zorder=6)
ax1.set_ylim(-1,1)
ax1.set_xlim(1950, 2019)    #연도별에서 추가
ax1.set_xticks(np.arange(1950,2020,5)) 
#ax1.set_xticklabels(tick_years)
ax1.set_ylabel('mtpr anomaly', fontsize=11)
ax1.set_title('(b) Mid lat (80-140E, 20-40N)', loc='left')
ax1.grid(alpha=0.3)
ax1.legend(loc='upper left', fontsize=9)
#(c) 고위도
ax1=fig.add_subplot(3, 1, 3)
ax1.plot(t2, mtpr_c_year, 'k', label='High lat', zorder=5)
ax1.plot(t2, mtpr_c_year_detrend, 'r', label="High lat Detrended", zorder=6)
ax1.set_ylim(-1,1)
ax1.set_xlim(1950, 2019)    #연도별에서 추가
ax1.set_xticks(np.arange(1950,2020,5)) 
#ax1.set_xticklabels(tick_years)
ax1.set_ylabel('mtpr anomaly', fontsize=11)
ax1.set_title('(c) High lat (80-140E, 40-60N)', loc='left')
ax1.grid(alpha=0.3)
ax1.legend(loc='upper left', fontsize=9)

plt.suptitle('mtpr anomaly (1950–2019) [mm/day]', fontsize=14, x=0.1, ha='left')
plt.subplots_adjust(bottom=0.15, top=0.9, left=0.08, right=0.97, hspace=0.25)

plt.savefig(i_path + 'ERA_mtpr_anom_yearly_ts_1950_2019.png', dpi=300)
plt.close()

