import meep as mp
from meep import materials as mat
import numpy as np
from matplotlib import pyplot as plt
from sys import exit

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
        self.sim = []
        self.n2f_obj = None
        self.ff_data = None
        self.incidence_flux_obj = None
        self.flux_frq = None
        self.nfrq = 11
        self.net_loss_flux_obj = None
        self.net_loss_flux_data = None
        self.ff_points = None 
        self.ff_angles = None

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

        # Wavelength selection (default 780nm)
        try:
            self.wvl = kwarg["wvl"] 
            if type( self.wvl ) == int:
                self.wvl = float( self.wvl )
        except:
            self.wvl = 0.780

        # Wavelength range selection (default 10nm)
        try:
            self.dwvl = kwarg["dwvl"] 
            if type( self.dwvl ) == int:
                self.dwvl = float( self.dwvl )
        except:
            self.dwvl = 0.01

        # The output file can have an identifiy ID atached to it
        try: 
            self.file_ID = kwarg["file_ID"]
        except:
            self.file_ID = ""

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
        elif type( self.wvl ) != float:
            raise TypeError("'wvl' must be a float")
        elif type( self.dwvl ) != float:
            raise TypeError("'dwvl' must be a float")
        elif type( self.file_ID ) != str:
            raise TypeError("'file' must be a string")

    def run( self, **kwarg ):
        frq = 1/self.wvl
        dfrq = 1/( self.wvl - 0.5*self.dwvl ) - 1/( self.wvl + 0.5*self.dwvl )
        dfrq = 0.5

        ## Build cell ###
        # chip size - the size of the chip in the x/y direction
        # padding - the distance between the sides of the chip and the PML layer
        # The x-z plane where the chip sits, the x direction give a cross section of the gratings
        # The y direction is the field propergation direction
        #################
        # Define the basic units for the simulation
        dpml = 0.5*self.wvl*3                                    # PML thickness
        chip_size_x = self.num_period*self.period           # Chip length in x
        chip_size_z = 0 if self.run_2D == True else 1       
        padding = 0                                         # Padding between the wall and the side of the chip
        plate_thickness = 1                                 # Thickness of the chip

        sx = dpml + padding + chip_size_x + padding + dpml
        sy = dpml + plate_thickness + self.grating_height + frq*5.0 + dpml
        sz = 0 if self.run_2D == True else chip_size_z + 2*padding + 2*dpml
        cell = mp.Vector3( sx, sy, sz )
        pml_layer = [ mp.PML( dpml ) ]

        # Create source
        source = [ mp.Source( mp.GaussianSource( frq, dfrq, is_integrated=True ),
                              component=self.polarization, 
                              center=mp.Vector3( y=0.5*sy-dpml ),
                              size=mp.Vector3( chip_size_x*0.98, z=chip_size_z ) ) ]

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
                raise
            else:
                symmetries = []
        except:
            symmetries = [ mp.Mirror( mp.X ) ]
            symmetries.append( mp.Mirror( mp.Z ) ) if self.run_2D == False else None


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
            self.sim = mp.Simulation( cell, self.res, geometry, source, boundary_layers=pml_layer, symmetries=symmetries )
            self.__animate_func__( sx, sy )
            return 0

        # If no animation is requested, the a normal run will comence
        # The simulation must be run twice, once with no geomitry - in order to remove the incoming 
        # field data from the n2f monitor
        self.sim = mp.Simulation( cell_size=cell,
                                  resolution=self.res,
                                  geometry=[],
                                  sources=source,
                                  boundary_layers=pml_layer,
                                  symmetries=symmetries )

        # This is the n2f for no geomitry
        n2f_point = mp.Vector3( y=-0.5*sy + dpml + plate_thickness + 1.05*self.grating_height )
        n2f_region = mp.Near2FarRegion( center=n2f_point, size=mp.Vector3( chip_size_x ), direction=mp.Y )
        n2f_obj_no_chip = self.sim.add_near2far( frq, dfrq, self.nfrq, n2f_region )

        # This is the incoming flux, it is idenical in shape to the n2f 
        incidence_flux_region = mp.Near2FarRegion( center=n2f_point, size=mp.Vector3( chip_size_x ), direction=mp.Y )
        self.incidence_flux_obj = self.sim.add_flux( frq, dfrq, self.nfrq, incidence_flux_region )

        # First run
        self.sim.run( until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,n2f_point,1e-12 ) )

        # Get the incoming n2f data
        n2f_data_no_chip = self.sim.get_near2far_data( n2f_obj_no_chip )

        # Get the near field flux and save if specified
        self.incidence_flux_data = mp.get_fluxes( self.incidence_flux_obj )
        
        self.sim.reset_meep()

        # Second run with the chip
        self.sim = mp.Simulation( cell_size=cell,
                                  resolution=self.res,
                                  geometry=geometry,
                                  sources=source,
                                  boundary_layers=pml_layer,
                                  symmetries=symmetries )
        

        # Add the near2far monitor then set to remove the incoming data
        self.n2f_obj = self.sim.add_near2far( frq, dfrq, self.nfrq, n2f_region )
        self.sim.load_minus_near2far_data( self.n2f_obj, n2f_data_no_chip )

        # Add flux monitor which is the net loss
        self.net_loss_flux_obj = self.sim.add_flux( frq, dfrq, self.nfrq, incidence_flux_region )

        # Run for a second time
        self.sim.run( until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,n2f_point,1e-12 ) )

        # Get the net flux data
        self.net_loss_flux_data = mp.get_fluxes( self.net_loss_flux_obj )

        return 0

    # The function builds the animation
    def __animate_func__( self, sx, sy, **kwarg ):
        animate = mp.Animate2D(self.sim,
                fields=mp.Ez,
                normilize=True,
                field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'none'},
                boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3},
                eps_parameters={'cmap':'binary'},
                #realtime=True,
                output_plane=mp.Volume( size=mp.Vector3( sx, sy ) ))

        self.sim.run(mp.at_every(0.5,animate), until_after_sources=mp.stop_when_fields_decayed( 5,mp.Ez, mp.Vector3(), 1e-6 ))
        
        animate.to_mp4( 6, 'anm.mp4' )
        self.sim.plot2D(fields=mp.Ez,
                   field_parameters={'alpha':0.8, 'cmap':'RdBu', 'interpolation':'none' },
                   boundary_parameters={'hatch':'o', 'linewidth':1.5, 'facecolor':'y', 'edgecolor':'b', 'alpha':0.3},
                   output_plane=mp.Volume( size=mp.Vector3( sx, sy ) ))
        plt.show()

    # Users invoke this request computaion of the far fields
    def get_far_field( self, ff_dist=5e3, ff_pnt=500, theta=np.pi/4, **kwarg ):
        if self.n2f_obj == None:
            raise RuntimeError( "Can't generate far fields without near field data first. Use 'run' first." )
        
        # Get the size of the far field and calculate is's resolution in pixles per micron
        ff_size = 2*abs(ff_dist)*np.tan( theta )
        ff_res = ff_pnt/ff_size

        # Generate the far field data
        self.ff_data = self.sim.get_farfields( self.n2f_obj,
                                               resolution=ff_res,
                                               center=mp.Vector3( y=ff_dist ),
                                               size=mp.Vector3( x=ff_size ) )

        # Save the position and angle data
        self.ff_points = np.linspace( -0.5*ff_size, 0.5*ff_size, ff_pnt )
        # fix line bellow
        self.ff_angles = np.arctan( self.ff_points/ff_dist )
        print(ff_size)

        return 0
    
    def save_ff_data( self, fname=None, **kwarg ):
        if fname == None or type( fname ) != str:
            raise TypeError( "'fname' must be a string" )

        import json

        try:
            data_file = open( fname, "w" )
        except:
            print("Could not save file")
            return 1

        output_data = self.ff_data
        output_data.update( { "angles": self.ff_angles,"points": self.ff_points } )

        output_data = { k:str(n.tolist()) for k,n in output_data.items() }

        output_data.update( {
            "num_period": self.num_period,
            "period": self.period, 
            "grating_width": self.grating_width,
            "grating_height": self.grating_height,
            "wvl": self.wvl } )

        json_data = json.dumps( output_data )
        del( output_data )

        data_file.write( json_data )
        del( json_data )

        data_file.close()
        del( data_file )

        return 0 

    # Plots the far field
    def plot_far_field( self, x_axis="angle", fname=None, dpi=300, **kwarg ):
        ff_p_vector = []
        ff_frq = mp.get_flux_freqs( self.n2f_obj )
        index = np.where( np.array( ff_frq ) == 1/self.wvl )[0][0]

        ff_p_vector = [ 
                        np.cross( np.array( [ self.ff_data['Ex'][i,index],
                                              self.ff_data['Ey'][i,index],
                                              self.ff_data['Ez'][i,index] ] ),
                                  np.array( [ self.ff_data['Hx'][i,index],
                                              self.ff_data['Hy'][i,index],
                                              self.ff_data['Hz'][i,index] ] ),
                                ) for i in range( len( self.ff_data['Ex'][:,0] ) ) ] 

        fig, axs = plt.subplots()
        axs.plot( self.far_field_angles, np.abs( self.ff_data )**2, '-k' )
        axs.tick_params( direction="in" )
        axs.grid( which="both" )
        axs.set_ylabel( "Poynting vector", size="x-large")
        axs.set_xlabel( "Angle (rad)", size="x-large") if x_axis == "angle" else axs.set_xlabel( "Far field position (μm)", size="x-large")

        # Save the file unless show plot if it is named, else just show
        plt.tight_layout()
        plt.savefig( fname, dpi=dpi ) if fname == None else plt.show()

    # Determines the efficny of the GMOT in regards to it's flux (power)
    def gmot_efficacy( self, **kwarg ):
        # Find the frequncy location for the wavelength
        flux_frq = mp.get_flux_freqs( self.incidence_flux_obj )
        index = np.where( np.array( flux_frq ) == 1/self.wvl )[0][0]

        # fix, (a-b)/a not b/a
        return abs( self.net_loss_flux_data[ index ]/self.incidence_flux_data[ index ] )

    """def n2f_efficacy( self, **kwarg ):
        # Compute the far field Poynting flux
        ff_flux = []

        for i in range( self.nfrq ):
            ff_flux.append( np.sum( [ 
                            np.real( np.cross( np.conj( np.array( [ self.ff_data['Ex'][j,i],
                                                  self.ff_data['Ey'][j,i],
                                                  self.ff_data['Ez'][j,i] ] )),
                                      np.array( [ self.ff_data['Hx'][j,i],
                                                  self.ff_data['Hy'][j,i],
                                                  self.ff_data['Hz'][j,i] ] ),
                                    ))[1] for j in range( len( self.ff_data['Ex'][:,0] ) ) ] ) )

        ff_frq = mp.get_near2far_freqs( self.n2f_obj )
        self.flux_frq = mp.get_flux_freqs( self.incidence_flux_obj )

        ff_index = np.where( np.array( ff_frq ) == 1/self.wvl )[0][0]
        flux_index = np.where( np.array( self.flux_frq ) == 1/self.wvl )[0][0]

        effic = np.abs( ff_flux[ ff_index ]/(self.incidence_flux_data[ flux_index ] ) )
        coff = self.ff_points[1] - self.ff_points[0]
        print( self.ff_points[1] - self.ff_points[0] )
        print( effic )

        return effic*coff"""

