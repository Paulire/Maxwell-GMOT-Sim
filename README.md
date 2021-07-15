# Maxwell Simulation using Meep for modelling GMOTs

# Class Refrance

```class linear_gmot( num_period, period, grating_width, grating_height, **kwarg )```

This simulation is initialsed by the usere here

 * ```num_period``` - number of periods to be modelled
 * ```period``` - peak-to-peak length of the grating
 * ```grating_width``` - length of grating trench
 * ```grating_hieght``` - height of the grating

There are also the following optional arguments

 * ```grating_material``` - Sets the grating matiral (default ```meep.Al```). Must be a ```meep.Medium``` to be accepted
 * ```res``` - Set the simulation's resolution (default 50 pixels/micron)
 * ```greating_func``` - the gratings are build by a function, the user may define their own (see billow at some point)
 * ```run_2D``` - Set True to run the simulation in 2D rather than the default 3D
 * ```polarization``` - Set the polerization (not added)
 * ```fname``` - Name of a json file (without the .json suffix) which contains existing simulation data

```def run( **kwarg )```

This is function builds and runs the simulation enviroment - work in process

 * ```symmetries``` - The symmetries by default are true, but it can be turn off (false) which will tripple the resorce use (double in the case of 2D simulations)
 * ```animate``` - Makes a 2D crossection animation of a single pulse 
 * ```n2f_file``` - Load existing n2f monitor object file from a HDF5 file (without the .h5 suffix) and does not run the simulation

```def get_far_field( ff_dist, ff_pnt, theta)```

Must be called after the run function, (unless an existing near2far monitor has been loaded with ```run( n2f_file="filename" )```) this funtion builds the far field data.

 * ```ff_dist``` - Distance to the far field
 * ```ff_pnt``` - The number of points which are sampled at the far field
 * ```theta``` - The maximum angle from normal incidence with the grating which is sampled a the far field

```def save_data( fname )```

Save the all the relevent simulation data for analysis and replicating the simulation settp.

 * ```fname``` - Name of a json file (without the .json suffix) to save the data to

```def save_n2f_obj( fname )```
 * ```fname``` - Name of HDF5 file (without the .h5 suffix) to save the near2far monitor to

```def gmot_efficacy( **kwarg )```

Returns the efficancy of the GMOT

```def load_data( fname, **kwarg )```
This is not appart of ```class linear_gmot```. This loads the json file containing the simulation data (including outpouts) are returns it as a ```dict```

* ```fname``` - Name of a json file (without the .json suffix) to load the data from
