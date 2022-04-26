# Visualization Examples

Geodata also provides the users with different methods to visualize outputs. 

To start, import the geodata package with a logger for detailed debugging.


```python
import geodata
import logging
logging.basicConfig(level=logging.INFO)
```

We also import the `geopandas` and `cartopy` libraries to retrieve and show geospatial [shapefiles](https://en.wikipedia.org/wiki/Shapefile) on the plot, and the `IPython` library to download generated animation as HTML file. These libaries are helpful, but not required to use geodata for visualization.


```python
import geopandas as gpd
import cartopy.io.shapereader as shpreader
from IPython.display import HTML
```

Download example datasets and create cutouts. We will get the hourly aerosol data and the hourly radiation data.


```python
#Download aerosol hourly data
aerosol_hourly_data = geodata.Dataset(module="merra2",
                                      years=slice(2020, 2020),
                                      months=slice(1,12),
                                      weather_data_config = "surface_aerosol_hourly")

#Download radiation hourly data
slv_hourly_data = geodata.Dataset(module = "merra2", 
                                  years = slice(2011, 2011),
                                  months = slice(1,1),
                                  weather_data_config = "slv_radiation_hourly")

if aerosol_hourly_data.prepared == False:
    aerosol_hourly_data.get_data()
    
#Download radiation hourly data only on 2011/01/01
if slv_hourly_data.prepared == False:
    slv_hourly_data.get_data(testing = True) 

#Create northern china aerosol Cutout
cutout_pm25 = geodata.Cutout(name="beijing19",
                             module="merra2",
                             weather_data_config="surface_aerosol_hourly",
                             xs=slice(105, 123),
                             ys=slice(27, 43),
                             years=slice(2019, 2019),
                             months=slice(1,12))

#Create china solar Cutout
cutout_solar = geodata.Cutout(name = "china-2011-slv-hourly-test",
                        module = "merra2",
                        weather_data_config = "slv_radiation_hourly",
                        xs = slice(73, 136), 
                        ys = slice(18, 54), 
                        years = slice(2011, 2011), 
                        months = slice(1,1))

cutout_solar.prepare()
cutout_pm25.prepare()
```

Generate pm25 and solar PV Outputs.


```python
ds_pm25 = geodata.convert.pm25(cutout_pm25)
ds_solar = geodata.convert.pv(cutout_solar, panel = "KANEKA", orientation = "latitude_optimal")
```

## Time Series Visualization

### Default time series method call

We can use `geodata.plot.time_series` to visualize time series data from the output xarray DataArray, such as `ds_pm25` or `ds_solar`. Its minimal method call find the mean value of all grid cell for every time point in the dataset. For example, with `ds_solar`, we can visualize the spatially aggregated averages AC power over time.


```python
geodata.plot.time_series(ds_solar)
```


![png](https://github.com/east-winds/geodata/blob/master/images/visualization/output_12_0.png)


### Spatial and temporal aggregation

The `time_series` method can take in tuple parameters `lat_slice` and `lon_slice` to select grid cells within that range (inclusive). For example, if we want to find the aggregated value for all grid cells between latitude 35 degree and 36 degree, we set `lat_slice` to be (35, 36). The `agg_slice_method` parameter will specify the aggregation method for aggregating grid cells sliced by `lat_slice` or `lon_slice`. By default, `agg_slice_method` is set to mean aggregation. 

We use the latitude-sliced time-series visualization on the PM2.5 output below. Note that since we have hourly data for the year 2019, we will have 24 * 365 = 8760 timepoints for each hour. However, we can reduce the number of timepoints by taking in a `time_factor` parameter that tells the method how many timepoints to aggregate on. Here, we take 24 * 7 as the `time_factor` so that we will aggregate the data by week, as there are 24 * 7 hours in a week. The `agg_time_method` parameter will specify the aggregation method for time aggregation. By default, `agg_time_method` is set to mean aggregation. 

For example, below we visualize the weekly averages of sum of PM2.5 for region within latitude slice (35, 36).


```python
geodata.plot.time_series(ds_pm25, 
                         lat_slice = (35, 36), agg_slice_method = 'sum', 
                         time_factor = 24 * 7)
```


![png](https://github.com/east-winds/geodata/blob/master/images/visualization/output_16_0.png)


If we have `lat_slice` or `lon_slice` inputs, and want to plot the time series for every single grid cell without aggregating them, they can specify `agg_slice = False`. This will generate one line for each grid cell.

The method also takes in user-defined title with the `title` parameter.


```python
geodata.plot.time_series(ds_pm25, lat_slice = (35, 36), lon_slice = (110, 111), 
                         agg_slice = False, 
                         time_factor = 24 * 7,
                         title = "PM2.5 Time Series - lat(35-36) lon(110-111) weekly average")
```


![png](https://github.com/east-winds/geodata/blob/master/images/visualization/output_18_0.png)


### Multiple coordinate points

You can also use a dictionary of name-coordinate pairs to plot different grid cells. The coordinates value of this `coord_dict` does not have to be exact, as the method can automatically find the grid cell containing the coordinate input.


```python
coord_d = {'Beijing': (30.9, 116.4),
           'Shanghai': (31.2, 121.47),
           'Xi\'an': (34.2, 108.9)}

geodata.plot.time_series(ds_pm25, coord_dict = coord_d, time_factor = 24 * 7)
```

    INFO:geodata.plot:Find grid cell containing coordinate for Beijing at lat = 30.5, lon = 116.25.
    INFO:geodata.plot:Find grid cell containing coordinate for Shanghai at lat = 31.0, lon = 121.25.
    INFO:geodata.plot:Find grid cell containing coordinate for Xi'an at lat = 34.0, lon = 108.75.
    


![png](https://github.com/east-winds/geodata/blob/master/images/visualization/output_21_1.png)


## Heatmap Visualization

### Default method call

Geodata can plot a spatial heatmap of output values. Since the output is a time-series containing more than 2 dimensions, this method will aggregate the values by mean at different timepoints for each grid cells by default. For example, to see the annual mean PM2.5 in our Cutout region, we use the following method call:


```python
geodata.plot.heatmap(ds_pm25)
```


![png](https://github.com/east-winds/geodata/blob/master/images/visualization/output_25_0.png)


### Add shapefiles to the plot

The `heatmap` method can also take in a `shape` parameter, which takes in a `geopandas` dataframe or series of shape objects. Let us use the province shapes from `cartopy` shape-reader and save the path as `prov_path`. This can also be the path to user-supplied shape files. 


```python
prov_path = shpreader.natural_earth(resolution='10m', category='cultural', 
                                    name = 'admin_1_states_provinces')
shapes = gpd.read_file(prov_path, encoding = 'utf-8')
```


```python
geodata.plot.heatmap(ds_pm25, shape = shapes)
```


![png](https://github.com/east-winds/geodata/blob/master/images/visualization/output_29_0.png)


### Selecting timepoint

If we do not want the temporally aggregated plot, we can specify the exact time point or its index in the dataArray. In the following method call, `t = 0` uses index to select the first time point in `ds_pm25`.


```python
geodata.plot.heatmap(ds_pm25, t = 0, shape = shapes)
```


![png](https://github.com/east-winds/geodata/blob/master/images/visualization/output_32_0.png)


We can also take in the exact time point from `ds_pm25` as a string. We can also change the map type from the default `colormesh` to `contour`, and customize the title text like the following:


```python
geodata.plot.heatmap(ds_pm25, t = '2019-01-01T00:30:00', 
                     map_type = 'contour', 
                     shape = shapes, 
                     title = 'Contour plot', title_size = 20)
```


![png](https://github.com/east-winds/geodata/blob/master/images/visualization/output_34_0.png)


Let's use the `heatmap` method on the solar PV output xarray `ds_solar`. Below we select the 7th time point for the `ds_solar` dataArray with the provincial shapes on the same plot.

Note that the default map color of the method is `bone_r`, which is not ideal for visualizing solar PV. Therefore, we switch the `cmap` parameter to `Wistia`. You can view a complete list of matplotlib map color [here](https://matplotlib.org/stable/gallery/color/colormap_reference.html).



```python
geodata.plot.heatmap(ds_solar, t = 6, shape = shapes, shape_width = 0.25, shape_color = 'navy',
                     map_type = 'contour', cmap = 'Wistia')
```


![png](https://github.com/east-winds/geodata/blob/master/images/visualization/output_36_0.png)


### Animation

The drawback of plotting a static heatmap with `heatmap` is that we cannot see the changes over time like the `time_series` plots. However, the `heatmap_animation` method can create an animation of heatmap with time as another dimension in the plot.

The parameters of the heatmap_animation is very similar to the ones for `heatmap`. You can use `time_factor` to find aggregated mean or sum. Here, we create the animation with averages for every two hours in the day.


```python
geodata.plot.heatmap_animation(ds_solar, cmap = 'Wistia', 
                               time_factor = 2, 
                               shape = shapes, shape_width = 0.25, shape_color = 'navy')
```
![png](https://github.com/east-winds/geodata/blob/master/images/visualization/pv_animation.gif)


The users can save the animation to a file, which requires the `HTML` method from the `IPython` package we imported earlier. It also requires the users to use the Jupyter Notebook in a browser, and have already generated the heatmap animation in the notebook, because `geodata.plot.save_animation` will extract the javascript content string from the animation in the Jupyter Notebook, and use HTML() method to enable the browser to download the file.

Save the animation above as a file named `solar_pv_2011_01_01_animation.html`.


```python
HTML(geodata.plot.save_animation("solar_pv_2011_01_01_animation.html"))
```




