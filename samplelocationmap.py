import pandas as pd
import geopandas
import matplotlib.pyplot as plt
import contextily as ctx


df = pd.DataFrame(
    {'Sample': ['P', 'C', 'E', 'J', 'V', 'U', '1', '2', '6'],
     'Latitude': [-2.358740, -16.462776,  -4.704314,  -8.921115,   -7.674723,  -6.720369,  -4.583105,   -7.820663, -16.853987],
     'Longitude': [-74.122498, -65.349432, -70.227259, -73.734416, -70.557612, -69.781242, -71.430290, -70.048701, -64.773992]})

gdf = geopandas.GeoDataFrame(
    df, geometry=geopandas.points_from_xy(df.Longitude, df.Latitude))
gdf = gdf.set_crs(epsg =4326)
#gdf = gdf.to_crs(epsg=3857)
world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
#world = world.to_crs(epsg=3857)
fig, ax = plt.subplots(1,1)
world[world.continent == 'South America'].plot(ax=ax,color='white', alpha = 0)
ctx.add_basemap(ax, source=ctx.providers.Stamen.Terrain, crs = 4326, attribution = False)

gdf.plot(ax=ax, color='black', markersize = 3, label = "_nolegend_")
gdf.plot(ax=ax, color='red', markersize = 1.5, label = "Sample Reaches")
ax.legend(loc = "upper right", fancybox = False, fontsize = "medium")
#ax.set_ylim(bottom = -20)
plt.savefig("samplelocation.png", dpi = 1500)
plt.close()

fig, (ax1, ax2) = plt.subplots(1,2)
ctx.add_basemap(ax, source=ctx.providers.Esri.OceanBasemap, crs = 4326, attribution = False)

