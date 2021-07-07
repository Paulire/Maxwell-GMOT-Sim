import meep as mp
from meep import materials as mat
import numpy as np
from matplotlib import pyplot as plt

def __default_greating__( self, cen, size_x, size_y ):
    geo = [ mp.Block( size=mp.Vector3( self.period - self.grating_width, size_y, self.grating_height ),
                      center=mp.Vector3( cen[0] + 0.5*( self.period - self.grating_width ), cen[1], 0.5*( 2*cen[2]+self.grating_height ) ),
                      material=self.grating_material ) ]

    return geo

class linear_gmot:
    def __init__( self,
                  num_period,
                  period,
                  grating_width,
                  grating_height,
                  **kwarg ):

        # Set the verlibles which must be defined
        self.num_period = int( num_period )              # Number of gratings
        self.period = float( period )                      # Size of grating
        self.grating_width = float( grating_width )        # Width of the deep part of the grating
        self.grating_height = float( grating_height )      # The hight of the grating

        # Set the optional verlibles
        # The grating matirial
        try:
            self.grating_material = kwarg["grating_material"]
        except:
            self.grating_material = mat.Al

        # Simulation resolution
        try:
            self.res = int( kwarg["res"] )
        except:
            self.res = 50

        # Simulation greating builder
        try:
            self.greating_func = kwarg["greating_func"]
        except:
            self.greating_func = __default_greating__

        # Simulation Polarization
        try:
            if kwarg["polarization"] == 'X':
                self.polarization = mp.Ex
            elif kwarg["polarization"] == 'Y':
                self.polarization = mp.Ey
            elif kwarg["polarization"] == 'LEFT':
                self.polarization = mp.Ey
            elif kwarg["polarization"] == 'RIGHT':
                self.polarization = mp.Ey
            else:
                self.polarization = None
        except:
            self.polarization = mp.Ey

        if self.polarization == None:
            raise ValueError("'polarization' must be either 'X', 'Y', 'LEFT' or 'RIGHT'")


        # Check if input data is correct
        if grating_width >= period:
            raise ValueError("'grating_width' cannot be equal to or grater than 'period'")
        elif grating_width <= 0:
            raise ValueError("'grating_width' must be a positive nonzero value")
        elif period <= 0: 
            raise ValueError("'period' must be a positive nonzero value")
        elif num_period <= 0:
            raise ValueError("'num_period' must be a positive nonzero value")
        elif grating_height <= 0:
            raise ValueError("'grating_height' must be a positive nonzero value")
        elif type( self.grating_material ) != mp.geom.Medium:
            raise TypeError("'grating_material' must be of type 'meep.geom.Medium'")
        elif type( self.greating_func ) != type( __default_greating__ ):
            raise TypeError("'greating_func' must be a function not " + str( type( __default_greating__ ) ))

    def run( self ):
        frq = 1
        dfrq = 1

        ### Build cell ###
        # chip size - the size of the chip in the x/y direction
        # padding - the distance between the sides of the chip and the PML layer
        #################
        dpml = 1
        chip_size = self.num_period*self.period + ( self.period - self.grating_width )
        padding = 2*self.period
        plate_thickness = 1

        sx = dpml + padding + chip_size + padding + dpml
        sy = sx
        sz = dpml + plate_thickness + self.grating_height + 3 + dpml
        cell = mp.Vector3( sx, sy, sz )
        pml_layer = [ mp.PML( dpml ) ]

        # Create source
        source = [ mp.Source( mp.GaussianSource( frq, dfrq ),
                              component=self.polarization, 
                              center=mp.Vector3( z=0.5*sz-dpml ),
                              size=mp.Vector3( chip_size, chip_size ) ) ]

        # Build greating
        geometry = [ mp.Block( size=mp.Vector3( chip_size, chip_size, plate_thickness ),
                               center=mp.Vector3( z=-0.5*( sz + plate_thickness ) + dpml ),
                               material=self.grating_material ) ]

        x_cen = np.linspace( -0.5*chip_size, 0.5*chip_size, self.num_period, endpoint=False) + 0.5*self.period
        z_cen = -0.5*sz + dpml + plate_thickness

        for i in range( self.num_period ):
            geometry.extend( self.greating_func( self, 
                                                 [ x_cen[i], 0, z_cen ],
                                                 self.period,
                                                 chip_size ) )


        symmetries = [ mp.Mirror( mp.Y ) ]
        sim = mp.Simulation( cell, self.res, geometry, source, boundary_layers=pml_layer, symmetries=symmetries )
        sim.run(until=10)


