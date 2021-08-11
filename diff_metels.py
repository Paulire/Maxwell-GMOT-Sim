#!/usr/bin/env ipython3

import lib.gmot as gmot
import meep.materials as mat
import meep as mp
import numpy as np
from matplotlib import pyplot as plt

# Define matirials which will be used
meterials = [ mp.Medium( index=1e20 ),          # High refractive index (reflectivity≈1)
              mat.Al,                           # Aluminium
              mat.Au,                           # Gold
              mat.Pd ]                          # Paladium
names = [ "high_n",
          "Al",
          "Au",
          "Pd" ]

# Loop for each Metal
for i in range( len( meterials ) ):
    # Set up the simulation enviroment
    A = gmot.linear_gmot( num_period=20,
                          period=1.080, #
                          grating_width=1.080*0.6,
                          grating_height=0.78/4,
                          wvl=0.675, #
                          dwvl=0.225, # 
                          nwvl=101, #
                          grating_material = meterials[i],
                          res=50, run_2D=True )
    
    # Run the simulation proper
    A.run()
    A.save_n2f_obj( names[i]  )

    # Generate the far field
    A.get_far_field( ff_dist=5e3, ff_pnt=1000, theta=4*np.pi/9 )

    # Get far field efficancy
    A.get_diffraction_efficacy( )
    A.plot_diffraction_efficacy( plot_total=True, fname="diff_metals/" + "/eff_" + str( names[i]  ) + ".pdf" )
    A.save_data( fname=str( names[i] )  )


