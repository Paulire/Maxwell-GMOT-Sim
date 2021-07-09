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

```def run( **kwarg )```

This is function builds and runs the simulation enviroment - work in process

 * ```symmetries``` - The symmetries by default are true, but it can be turn off (false) which will tripple the resorce use (double in the case of 2D simulations)
 * ```animate``` - Makes a 2D crossection animation of a single pulse 
