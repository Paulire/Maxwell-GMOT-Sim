import meep as mp
from meep import materials as mat
import numpy as np
from matplotlib import pyplot as plt

def __default_greating__( self, cen, size_x, size_z ):
    slit_len = self.period - self.grating_width

    geo = [ mp.Block( size=mp.Vector3( slit_len, self.grating_height, size_z ),
                      center=mp.Vector3( cen[0],  0.5*( 2*cen[1]+self.grating_height), cen[2] ),
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

        # The simulation object will be here and can be accessed anywhere
        sim = []

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

        # Should this be a 2D simulation
        try:
            self.run_2D = kwarg["run_2D"]
        except:
            self.run_2D = False

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
            self.polarization = mp.Ez

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
        elif type( self.run_2D ) != bool:
            raise TypeError("'run_2D' must a 'bool'")

    def run( self, **kwarg ):
        wvl = 1.0
        frq = 1/wvl
        dfrq = 0.5

        ### Build cell ###
        # chip size - the size of the chip in the x/y direction
        # padding - the distance between the sides of the chip and the PML layer
        # The x-z plane where the chip sits, the x direction give a cross section of the gratings
        # The y direction is the field propergation direction
        #################
        dpml = 0.5*wvl*3
        chip_size_x = self.num_period*self.period
        chip_size_z = 0 if self.run_2D == True else 1
        padding = 2*self.period*0
        plate_thickness = 1

        sx = dpml + padding + chip_size_x + padding + dpml
        sy = dpml + plate_thickness + self.grating_height + frq*5.0 + dpml
        sz = 0 if self.run_2D == True else chip_size_z + 2*padding + 2*dpml
        cell = mp.Vector3( sx, sy, sz )
        pml_layer = [ mp.PML( dpml ) ]

        # Create source
        source = [ mp.Source( mp.GaussianSource( frq, dfrq ),
                              component=self.polarization, 
                              center=mp.Vector3( y=0.5*sy-dpml ),
                              size=mp.Vector3( chip_size_x, z=chip_size_z ) ) ]

        # The greating base chip
        geometry = [ mp.Block( size=mp.Vector3( chip_size_x, plate_thickness, chip_size_z),
                               center=mp.Vector3( y=-0.5*( sy - plate_thickness ) + dpml ),
                               material=self.grating_material ) ]

        x_cen = np.linspace( -0.5*chip_size_x, 0.5*chip_size_x, self.num_period, endpoint=False)
        x_cen += 0.5*(x_cen[1] - x_cen[0])
        y_cen = -0.5*sy + dpml + plate_thickness

        # The gratings are built here
        for i in range( self.num_period ):
            geometry.extend( self.greating_func( self, 
                                                 [ x_cen[i], y_cen, 0 ],
                                                 self.period,
                                                 chip_size_z ) )


        # Add symitries if not stated or is stated
        try:
            if kwarg["symmetries"] == True:
                symmetries = [ mp.Mirror( mp.X ) ]
                symmetries.append( mp.Mirror( mp.Z ) ) if self.run_2D == False else None
            else:
                symmetries = []
        except:
            symmetries = [ mp.Mirror( mp.X ) ]
            symmetries.append( mp.Mirror( mp.Z ) ) if self.run_2D == False else None

        self.sim = mp.Simulation( cell, self.res, geometry, source, boundary_layers=pml_layer, symmetries=symmetries )

        # Check if the user has specified an animation
        try:
            if kwarg["animate"] == True:
                animate = True
            else:
                animate == False
        except:
            animate = False
        # If so, then run the animation
        if animate == True:
            self.__animate_func__( sx, sy )
            return 0

        # If no animation is requested, the a normal run will comence
        # First the n2f monitor is added
        n2f_point = mp.Vector3( y=-0.5*sy + dpml + plate_thickness + 1.05*self.grating_height )
        n2f_region = mp.Near2FarRegion( center=n2f_point, size=mp.Vector3( chip_size_x ) )
        n2f_obj = self.sim.add_near2far( frq, dfrq, 100, n2f_region )
        

        self.sim.run(until=7)
        self.sim.plot2D(fields=mp.Ez,
                   field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'none' },
                   boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3},
                   output_plane=mp.Volume( size=mp.Vector3( sx, sy ) ))
        plt.show()

    # The function builds the animation
    def __animate_func__( self, sx, sy, **kwarg ):
        animate = mp.Animate2D(self.sim,
                fields=mp.Ez,
                normilize=True,
                field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'none'},
                boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3},
                eps_parameters={'cmap':'binary'},
                output_plane=mp.Volume( size=mp.Vector3( sx, sy ) ))

        self.sim.run(mp.at_every(0.5,animate), until_after_sources=mp.stop_when_fields_decayed( 5,mp.Ez, mp.Vector3(), 1e-6 ))
        
        animate.to_mp4( 6, 'anm.mp4' )

