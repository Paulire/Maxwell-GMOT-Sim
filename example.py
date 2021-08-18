#!/usr/bin/env ipython3

import lib.gmot as gmot
import meep.materials as mat
import meep as mp
from sys import exit
import numpy as np
from matplotlib import pyplot as plt

# Solid metal
A = gmot.linear_gmot( num_period=4,             # Number of periods to be sampled
                      period=1.080,             # The period
                      grating_width=1.080*0.6,  # The width of the etching
                      grating_height=0.78/4,    # Etch depth
                      coating_height=0.05,      # Thickness of metal coating
                      wvl=0.700,                # Mean wavlength
                      dwvl=0.200,               # Wavelength range ( so 600-800nm )
                      nwvl=101,                 # Number of wavlengths sampled
                      grating_material=mat.Al,  # The matirial the grating is made of (without this line it is Al by default)
                      res=50, run_2D=True )

A.run( )#plot_settup=True )
A.get_far_field( ff_dist=5e3, ff_pnt=825, theta=4*np.pi/9)
A.get_diffraction_efficacy(  )
A.save_data( 'out' )
A.default_output( "temp/data", include_animation=True, include_settup=True ) 

####################################################################################################
# Example of sim with silicon and metal layer gmot
A = gmot.linear_gmot( num_period=4,             
                      period=1.080,             
                      grating_width=1.080*0.6,  
                      grating_height=0.78/4,    
                      coating_height=0.05,      
                      wvl=0.700,                
                      dwvl=0.200,
                      nwvl=101,                 
                      greating_func=gmot.silicon_metal_coated, # Swtiches to the built in Si Metal mot
                      grating_material=mat.cSi,                # This ensures the chip is silicon
                      res=50, run_2D=True )
