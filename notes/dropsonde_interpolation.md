# How to create a vertical transect from dropsonde data
While the code that performs dropsonde gridding is too cumbersome and should be significantly refactored and documented, in the meantime this guide can be useful to someone who works with dropsonde or radiosonde observations.
 
1) As shown in cell 29 of [Fig5_Fig6_ds_vs_um.ipynb](Fig5_Fig6_ds_vs_um.ipynb) notebook, the 2D cross-section is produced by `ds_cross_sect()` function from `misc_utils` module.
 
2) For that I need `leg_grid` and `leg_dist` arrays, which are produced on the previous line by `get_leggrid()` function.
That [function](../misc_utils.py#L186) takes a list of "dropsonde instances". In that old code I used custom classes from `faamtools.avaps` [package](https://github.com/dennissergeev/faamtools/blob/master/faamtools/avaps.py#L56), but if I were to write this code from scratch today, I would definitely use pandas DataFrames. You can notice that in `get_leggrid()` function (and others) I use **not** the "raw" data, but the "fil" data (short for "filtered"), which are columns of dropsonde data with rows that have at least one missing value, removed. Again, I kind of reinvented the wheel here, because nowadays pandas package allows you to do it much easier.
Anyway, back to the `get_leggrid()` function. Essentially, it calculates a 1D array of distances (in km) along a line of dropsonde locations. The locations are taken as mean longitude and latitude, and then a simple great circle formula is used to calculate distance on a sphere.
 
So if I had 6 dropsondes, the outputs would look like e.g.
```python
leg_dist = [0., 12.,  34., 56., 77., 89.]  # numpy array of shape (6, ). So the distance between first and second dropsonde is 12 km, etc.
leg_grid = [0., 1., 2., 3., ... 88., 89]  # numpy array from min to max, evenly spaced by `xstep`, by default 1 km.
```
The `leg_grid` array is the target x-grid that you want to interpolate to.
 
3) These two arrays are passed to `ds_cross_sect()` function. This function's signature looks like this:
   ```python
   ds_cross_sect(xgrid, zz, leg_dist, ds_id, data, coordname, varname, filt_kw):
       ...
   ```
  When it's called in my notebook,
   - `leg_grid` is passed as `xgrid`,
   - `zz` is a 1D array of target altitudes that you want to interpolate to (z-grid),
   - `leg_dist`,
   - `ds_id` is index of dropsonde,
   - `data` is list of dropsonde instances,
   - `coordname` is the name of vertical coordinate for interpolation (I'm using "hgt", height),
   - `varname` is the variable to interpolate, e.g. temperature
   - `filt_kw` is a dictionary of filter keywords that are passed to `misc_utils.blf()` function, which is built on top of `scipy.signal.butter()` function (for smoothing).
 
The first step is to call `ds_regz()` [function](../misc_utils.py#L168), which creates a 2D array of dropsonde variable. This array `regz` is regular in height (defined by `zz`) and has shape of (*number of height levels X number of dropsondes*). Also, the Butterworth lowpass filter is applied on the line 157, if filtering dictionary is given.
The interpolation itself is done using scipy.interpolate.interp1d() function.
 
OK, back to ds_cross_sect(), line 170. Now a 2D array of zeros is prepared for the final cross-section. The remaining step is to interpolate each slice of the second axis of the regz array to regular x-grid, i.e. go from array with (*number of height levels X number of dropsondes*) shape to the array (*number of height levels X number of x-axis points*). This is done by iterating number-of-height-levels times on lines [172-174](../misc_utils.py#L172).
 
 
### TL;DR
Your steps would probably be:
- create regular x-grid (1D)
- calculate distances of dropsondes along the cross-section
- create regular z-grid (1D)
- interpolate each of the dropsondes to z-grid (and optionally smooth out the data by butterworth-like filter)
- Interpolate each regular-height slice to x-grid points
All you need is `scipy.interpolate.interp1d()` function. With pandas DataFrames all this could be done quite quickly I think.
