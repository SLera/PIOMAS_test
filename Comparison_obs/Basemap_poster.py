import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import sys
from matplotlib import font_manager
font_manager._rebuild()
import matplotlib
matplotlib.rcParams['font.family'] = 'Montserrat'
sys.path.append('/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test/functions/')
from matplotlib.patches import Polygon


def draw_screen_poly( line, fill, lats, lons, m, hatch):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, edgecolor= line, facecolor=fill, hatch= hatch, alpha=0.4)
    plt.gca().add_patch(poly)

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

INDIR = '/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/'
OUTDIR = '/home/valeria/NIERSC/Scripts/IceVolume/PIOMAS_test_results/Comparison_obs/figs/Greenland/Hi_mean_monthly/'

chi_months = np.load(INDIR+'CryosatHi_full_sorted_monthly')
chi_u_months = np.load(INDIR+'CryosatUnc_full_sorted_monthly')
phi_months = np.load(INDIR+'PIOMAS_full_sorted_monthly')
clat = np.load(INDIR+'clat_full')
clon = np.load(INDIR+'clon_full')

cor_array = np.load(INDIR+'PIOMAS_Cryosat_correlation_array_full')
months = np.array([10,11,12,1,2,3,4])

cor = cor_array[220:370,80:255]
clat = clat[220:370,80:255]
clon = clon[220:370,80:255]

#PLOR REGION OF INTEREST

fsize = cm2inch((20,20))
plt.figure(figsize=fsize)

title = 'test_map'
outfname = OUTDIR + 'Regions_Greenland_Irminger.png'

m = Basemap(resolution="i",
            projection='laea', lat_ts=90, lat_0=90., lon_0=0.,
            llcrnrlon=clon[-1, 0], llcrnrlat=clat[-1, 0],
            urcrnrlon=clon[0, -1], urcrnrlat=clat[0, -1])
m.drawmeridians(np.arange(-180, 180, 10), labels=[0,0,0,1], fontsize=12)
m.drawparallels(np.arange(30, 80, 10), labels=[1,0,1,0], fontsize=12)


#Greenland Sea polygon
#PIOMAS GRID
# [66,83,320,25]
lat0 =66
lat1 = 83
lon0 = -40
lon1 = 25
resolution = 100

lats1 = np.linspace( lat0, lat0, resolution )
lons1 = np.linspace( lon0, lon1, resolution )

lats2 = np.linspace( lat0, lat1, resolution )
lons2 = np.linspace( lon1, lon1, resolution )

lats3 = np.linspace( lat1, lat1, resolution )
lons3 = np.linspace( lon1, lon0, resolution )

lats4 = np.linspace( lat1, lat0, resolution )
lons4 = np.linspace( lon0, lon0, resolution )

lats = np.concatenate((lats1,lats2,lats3,lats4))
lons = np.concatenate((lons1,lons2,lons3,lons4))
line = 'red'
fill = 'red'
hatch = ''
draw_screen_poly(line, fill, lats, lons, m, hatch )

#Irminger_Labrador Sea polygon
#PIOMAS GRID
# [49,66,300,340]

lat0 = 55
lat1 = 66
lon0 = -60
lon1 = -20
resolution = 100

lats1 = np.linspace( lat0, lat0, resolution )
lons1 = np.linspace( lon0, lon1, resolution )

lats2 = np.linspace( lat0, lat1, resolution )
lons2 = np.linspace( lon1, lon1, resolution )

lats3 = np.linspace( lat1, lat1, resolution )
lons3 = np.linspace( lon1, lon0, resolution )

lats4 = np.linspace( lat1, lat0, resolution )
lons4 = np.linspace( lon0, lon0, resolution )

lats = np.concatenate((lats1,lats2,lats3,lats4))
lons = np.concatenate((lons1,lons2,lons3,lons4))
line = 'blue'
fill = 'blue'
hatch = ''
draw_screen_poly(line, fill, lats, lons, m, hatch )

# FRAM
lat0 =82
lat1 = 80.6
lon0 = -12
lon1 = 20
resolution = 100

lats1 = np.linspace( lat0, lat0, resolution )
lons1 = np.linspace( lon0, lon1, resolution )

lats2 = np.linspace( lat0, lat1, resolution )
lons2 = np.linspace( lon1, lon1, resolution )

lats3 = np.linspace( lat1, lat1, resolution )
lons3 = np.linspace( lon1, lon0, resolution )

lats4 = np.linspace( lat1, lat0, resolution )
lons4 = np.linspace( lon0, lon0, resolution )

lats = np.concatenate((lats1,lats2,lats3,lats4))
lons = np.concatenate((lons1,lons2,lons3,lons4))
line = '#3F5D7D'
fill = '#3F5D7D'
hatch = ''
draw_screen_poly(line, fill, lats, lons, m, hatch )

m.fillcontinents(color='#3F5D7D', lake_color='#3F5D7D')
plt.title(title)
plt.savefig(outfname)
plt.close()