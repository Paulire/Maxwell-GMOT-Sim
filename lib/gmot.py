import meep as mp
from meep import materials as mat
import numpy as np
from matplotlib import pyplot as plt
from sys import exit
import json

def __default_greating__( self, cen, size_x, size_z ):
    slit_len = self.period - self.grating_width

    geo = [ mp.Block( size=mp.Vector3( slit_len, self.grating_height, size_z ),
                      center=mp.Vector3( cen[0],  0.5*( 2*cen[1]+self.grating_height), cen[2] ),
                      material=self.grating_material ) ]

    return geo


# The loading system is defined here in order to allow the user to use it indipendently 
# of the simulation enviroment
def load_data( fname=None, **kwarg ):
    # Check if the file name is a string
    if fname == None or type( fname ) != str:
        raise TypeError( "'fname' must be a string" )

    # Load data from file
    try:
        input_file = open( fname )
    except:
        raise IOError( "Could not find the file " + fname )
    input_data = json.load( input_file )
    input_file.close()

    # Convert arrays to the correct format (str->list)
    input_data['Ez'] = eval( input_data['Ez'] )
    input_data['Ey'] = eval( input_data['Ey'] )
    input_data['Ex'] = eval( input_data['Ex'] )
    input_data['Hz'] = eval( input_data['Hz'] )
    input_data['Hy'] = eval( input_data['Hy'] )
    input_data['Hx'] = eval( input_data['Hx'] )
    input_data['angles'] = eval( input_data['angles'] )
    input_data['points'] = eval( input_data['points'] )

    return input_data

class linear_gmot:
    def __init__( self,
                  num_period=None,
                  period=None,
                  grating_width=None,
                  grating_height=None,
                  fname=None,
                  **kwarg ):

        # If a file is specified, the simulation will attempt to read it
        if fname != None:
            # Load the file and get the dict containing the sim data
            sim_data = load_data( fname=fname )

            # Set the basic data
            self.num_period = int( sim_data["num_period"] )
            self.period = float( sim_data["period"] )
            self.grating_width = float( sim_data["grating_width"] )
            self.grating_height = float( sim_data["grating_height"] )
            self.wvl = float( sim_data["wvl"] )
            self.nwvl = int( sim_data["nwvl"] )
            self.dwvl = float( sim_data["dwvl"] )
            self.res = int( sim_data["res"] )
            self.run_2D = bool( sim_data["run_2D"] )

            # Get far field data from loaded file
            self.ff_points = np.array( sim_data['points'] )
            self.ff_angles = np.array( sim_data['angles'] )
            self.ff_data = { key:np.array(sim_data[key]) for key in ['Ex','Ey','Ez','Hx','Hy','Hz'] }

            self.diff_efficacy = np.array( sim_data["diff_efficacy"] )

            self.incidence_flux_data = np.array( sim_data["incidence_flux_data"] )

            self.ff_dist = sim_data["ff_dist"]

            self.frq_values = np.array( sim_data["frq_values"] )

            self.incidence_flux_data = np.array( sim_data[ "incidence_flux_data" ] )
            self.net_loss_flux_data = np.array( sim_data[ "net_loss_flux_data" ] )

        # If no file is specifed, then the code shall use the input arguments
        else:
            # Set the verlibles which must be defined
            self.num_period = int( num_period )              # Number of gratings
            self.period = float( period )                      # Size of grating
            self.grating_width = float( grating_width )        # Width of the deep part of the grating
            self.grating_height = float( grating_height )      # The hight of the grating

            # The simulation object will be here and can be accessed anywhere
            self.ff_data = None
            self.ff_points = None 
            self.ff_angles = None
            self.diff_efficacy = np.array( [] )
            self.incidence_flux_data = None
            self.ff_dist = None
            self.frq_values = None
            self.net_loss_flux_data = None

            # Simulation resolution
            try:
                self.res = int( kwarg["res"] )
            except:
                self.res = 50

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

            # The number of wavelenths to be sampled
            try:
                self.nwvl = kwarg["nwvl"]
                if type( self.nwvl ) == float:
                    self.nwvl = int( self.nwvl )
            except:
                self.nwvl = 11

        # Simulation greating builder
        try:
            self.greating_func = kwarg["greating_func"]
        except:
            self.greating_func = __default_greating__

        # Set the optional verlibles
        # The grating matirial
        try:
            self.grating_material = kwarg["grating_material"]
        except:
            self.grating_material = mat.Al

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


        self.sim = []
        self.n2f_obj = None
        self.incidence_flux_obj = None
        self.flux_frq = None
        self.net_loss_flux_obj = None

        # Check if input data is correct
        if self.grating_width >= self.period:
            raise ValueError("'grating_width' cannot be equal to or grater than 'period'")
        elif self.grating_width <= 0:
            raise ValueError("'grating_width' must be a positive nonzero value")
        elif self.period <= 0: 
            raise ValueError("'period' must be a positive nonzero value")
        elif self.num_period <= 0:
            raise ValueError("'num_period' must be a positive nonzero value")
        elif self.grating_height <= 0:
            raise ValueError("'grating_height' must be a positive nonzero value")
        elif type( self.grating_material ) != mp.geom.Medium:
            raise TypeError("'grating_material' must be of type 'meep.geom.Medium'")
        elif type( self.greating_func ) != type( __default_greating__ ):
            raise TypeError("'greating_func' must be a function not " + str( type( __default_greating__ ) ))
        elif type( self.run_2D ) != bool:
            raise TypeError("'run_2D' must a 'bool'")
        elif type( self.wvl ) != float:
            try:
                self.wvl = float( self.wvl )
            except:
                raise TypeError("'wvl' must be a float")
        elif type( self.dwvl ) != float:
            raise TypeError("'dwvl' must be a float")
        elif type( self.nwvl  ) != int:
            raise TypeError("'nwvl' must be an int")
        elif self.polarization == None:
            raise ValueError("'polarization' must be either 'X', 'Y', 'LEFT' or 'RIGHT'")

    def run( self, **kwarg ):
        frq = 1/self.wvl
        wvl_max = self.wvl + 0.5*self.dwvl
        wvl_min = self.wvl - 0.5*self.dwvl
        frq_max = 1/wvl_min
        frq_min = 1/wvl_max
        dfrq = frq_max - frq_min
        frq_cen = 0.5*( frq_max + frq_min )

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
        #pml_layer = [ mp.Absorber( dpml ) ]

        # Create source
        try:
            # If the user want't a monochromatic source, it can be requested here
            # The monochromatic source is a continus pulse which emmits for 5 periods
            if kwarg['mono_chrome'] == True:
                source = [ mp.Source( mp.ContinuousSource( frq, end_time=5*self.wvl),
                                      component=self.polarization, 
                                      #center=mp.Vector3( y=0.5*sy-dpml ),
                                      center=mp.Vector3( y=0 ),
                                      size=mp.Vector3( chip_size_x*0.98, z=chip_size_z ) ) ]
            else:
                raise TypeError("")
        except:
            # If monochrome is not specified then a guassian source is used
            source = [ mp.Source( mp.GaussianSource( wavelength=self.wvl, fwidth=self.dwvl, is_integrated=True ),
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

        # Check if the user wants to build the enviroment for n2f analysis form an od life
        try:
            if kwarg["n2f_file"]:
                n2f_only = True
        except:
            n2f_only = False

        if n2f_only == True:
            if type( kwarg["n2f_file"] ) != str:
                raise TypeError( "'n2f_file must be a sting" )
            
            self.sim = mp.Simulation( cell_size=cell,
                                      resolution=self.res,
                                      geometry=geometry,
                                      sources=source,
                                      boundary_layers=pml_layer,
                                      symmetries=symmetries )

            n2f_point = mp.Vector3( y=-0.5*sy + dpml + plate_thickness + 1.10*self.grating_height )
            n2f_region = mp.Near2FarRegion( center=n2f_point, size=mp.Vector3( chip_size_x ), direction=mp.Y )
            self.n2f_obj = self.sim.add_near2far( frq, dfrq, self.nwvl, n2f_region )
            self.sim.load_near2far( kwarg["n2f_file"], self.n2f_obj )
            
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
        n2f_point = mp.Vector3( y=-0.5*sy + dpml + 1.05*( plate_thickness + self.grating_height ) )
        n2f_region = mp.Near2FarRegion( center=n2f_point, size=mp.Vector3( chip_size_x ), direction=mp.Y )
        n2f_obj_no_chip = self.sim.add_near2far( frq, dfrq, self.nwvl, n2f_region )

        # This is the incoming flux, it is idenical in shape to the n2f 
        incidence_flux_region = mp.Near2FarRegion( center=n2f_point, size=mp.Vector3( chip_size_x ), direction=mp.Y )
        self.incidence_flux_obj = self.sim.add_flux( frq, dfrq, self.nwvl, incidence_flux_region )

        # First run
        self.sim.run( until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,n2f_point,1e-12 ) )

        # Get the incoming n2f data
        n2f_data_no_chip = self.sim.get_near2far_data( n2f_obj_no_chip )

        # Get the near field flux and save if specified
        self.incidence_flux_data = np.array( mp.get_fluxes( self.incidence_flux_obj ) )
        
        self.sim.reset_meep()

        # Second run with the chip
        self.sim = mp.Simulation( cell_size=cell,
                                  resolution=self.res,
                                  geometry=geometry,
                                  sources=source,
                                  boundary_layers=pml_layer,
                                  symmetries=symmetries )

        # Add the near2far monitor then set to remove the incoming data
        self.n2f_obj = self.sim.add_near2far( frq, dfrq, self.nwvl, n2f_region )
        self.sim.load_minus_near2far_data( self.n2f_obj, n2f_data_no_chip )

        # Add flux monitor which is the net loss
        self.net_loss_flux_obj = self.sim.add_flux( frq, dfrq, self.nwvl, incidence_flux_region )

        # Run for a second time
        self.sim.run( until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,n2f_point,1e-12 ) )

        # Get the net flux data
        self.net_loss_flux_data = np.array( mp.get_fluxes( self.net_loss_flux_obj ) )

        # Store the flux frequncies
        self.frq_values =  np.array( mp.get_near2far_freqs( self.n2f_obj ) ) 

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

        self.ff_dist = ff_dist
        
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
        self.ff_angles = np.arctan( self.ff_points/ff_dist )

        return 0

    def get_far_field_efficacy( self, order=1, **kwarg ):
        # Validate input values
        if type( order ) != int:
            return TypeError( "'order' must be an int" )
        elif order < 0:
            return ValueError("'order must be greater than or equle to zero'")

        flux = []

        # Find the efficancy for each wavelength
        for i in range( len( np.array( self.frq_values  ) ) ):
            # Manually find's the pynting vector and uses it's magnitude as a diffraction finder
            Px = np.real( np.conj( self.ff_data['Ey'][:,i] )*self.ff_data['Hz'][:,i] - np.conj( self.ff_data['Ez'][:,i] )*self.ff_data['Hy'][:,i] ) 
            Py = np.real( np.conj( self.ff_data['Ez'][:,i] )*self.ff_data['Hx'][:,i] - np.conj( self.ff_data['Ex'][:,i] )*self.ff_data['Hz'][:,i] ) 
            Pz = np.real( np.conj( self.ff_data['Ex'][:,i] )*self.ff_data['Hy'][:,i] - np.conj( self.ff_data['Ey'][:,i] )*self.ff_data['Hx'][:,i] ) 

            Pv = np.sqrt( Px**2 + Py**2 + Pz**2 )

            flux.append( np.sum( Py )*(self.ff_points[1]-self.ff_points[0] ) )
            
        flux = np.array(flux)
        wvl = 1000*np.divide(1, self.frq_values )

        plt.plot( wvl, np.abs( flux)/( np.abs( self.incidence_flux_data ) - np.abs( self.net_loss_flux_data ) ), '-r' )


        """meep_flux = self.n2f_obj.flux( mp.Y,
                           where=mp.Volume( center = mp.Vector3( y=self.ff_dist ),
                                            size = mp.Vector3( self.ff_points[-1] - self.ff_points[1] ) ),
                           resolution=1
                           )

        meep_flux = np.array( meep_flux )
        plt.plot( wvl, np.abs( meep_flux)/( np.abs( self.incidence_flux_data ) - np.abs( self.net_loss_flux_data ) ), '-.g' )"""


        plt.xlabel("Wavelength (nm)", size="x-large")
        plt.xticks(fontsize='large')
        plt.ylabel("Far and Near flux ratio", size="x-large")
        plt.yticks(fontsize='large')
        plt.tick_params( direction='in', length=4 )
        plt.minorticks_on()
        plt.tick_params( which='minor', length=2, direction='in'  )
        plt.savefig( "near_far_flux_v2.pdf", dpi=300 )
        plt.show()

            



        return 0

    def get_diffraction_efficacy( self, order=1, **kwarg ):
        # Validate input values
        if type( order ) != int:
            return TypeError( "'order' must be an int" )
        elif order < 0:
            return ValueError("'order must be greater than or equle to zero'")


        # Retreaves the frequncy values for the far fields
        index = np.abs( self.frq_values - 1/self.wvl ).argmin()

        # Diffraction efficacies are held hear for each frequncy
        # 0 = order zero, 1 = order one, 2 = order minus one
        self.diff_efficacy = np.zeros( ( 3, len( np.array( self.frq_values  ) ) ) )

        # Set up loop counters
        count_frq = 1

        # Find the efficancy for each wavelength
        for i in range( len( np.array( self.frq_values  ) ) ):
            print( "get_diffraction_efficacy working on frequency " + str( count_frq ) + " of " + str( len( np.array( self.frq_values  )  ) + 1 ) + " (" + str( int( 100*count_frq/( len( np.array( self.frq_values  ) ) )  )) + "% done)" ) 
            count_frq += 1


            # Calulates each maximas angl
            wvl = 1/self.frq_values[i]
            diffraction_angles_mean = [ np.arcsin( j*wvl/self.period ) for j in range( -order, order+1 ) ]
            position_mean = self.ff_dist*np.tan( diffraction_angles_mean )

            # Find the efficancy for each maxima at this wavelength
            for j in range( len( diffraction_angles_mean ) ):
                # Manually find's the pynting vector and uses it's magnitude as a diffraction finder
                Px = np.real( np.conj( self.ff_data['Ey'][:,i] )*self.ff_data['Hz'][:,i] - np.conj( self.ff_data['Ez'][:,i] )*self.ff_data['Hy'][:,i] ) 
                Py = np.real( np.conj( self.ff_data['Ez'][:,i] )*self.ff_data['Hx'][:,i] - np.conj( self.ff_data['Ex'][:,i] )*self.ff_data['Hz'][:,i] ) 
                Pz = np.real( np.conj( self.ff_data['Ex'][:,i] )*self.ff_data['Hy'][:,i] - np.conj( self.ff_data['Ey'][:,i] )*self.ff_data['Hx'][:,i] ) 

                Pv = np.sqrt( Px**2 + Py**2 + Pz**2 )
                field_magnitude = np.abs( Py )**2

                # Find the position of this cenre
                mean_index = int( ( position_mean[j] - self.ff_points[0] )/( self.ff_points[1] - self.ff_points[0] ) )

                # This loop will find the minimum field to the right of the first order diffraction
                # It records the index of the field at that point
                current_index_right = mean_index
                while field_magnitude[ current_index_right ]/field_magnitude[ mean_index ] > 0.01:
                    current_index_right +=1
                current_index_left = mean_index
                while field_magnitude[ current_index_left ]/field_magnitude[ mean_index ] > 0.01:
                    current_index_left -= 1
                
                # Determine the width of the maxima and it's mean position
                maxima_width = self.ff_points[ current_index_right ] - self.ff_points[ current_index_left ]
                maxima_mean = (self.ff_points[ current_index_right ] + self.ff_points[ current_index_left ])/2

                # Calculate flux
                self.diff_efficacy[j,i] =  np.sum( Py[ current_index_left:current_index_right ] )*( self.ff_points[1] - self.ff_points[0] )


        # Convert to effiancy
        self.diff_efficacy = self.diff_efficacy/self.incidence_flux_data

        return 0

    def plot_diffraction_efficacy( self, plot_total=False, fname=None, dpi=300, **kwarg ):
        # Get the frequncies
        # Convert to wavelengths
        ff_wvl = 1e3*np.divide( 1, self.frq_values )

        # Pre set style for plot
        line_style = [ '-r','-b','-.g','-k' ]
        leg = [ '$\eta_{-1}$', '$\eta_0$', '$\eta_1$', '$3\eta_1/(1-\eta_0)$' ]

        # Need to be renamebed, this is the analytical relation for the efficancies
        combine = 3*self.diff_efficacy[2,:]/( 1 - self.diff_efficacy[1,:] )

        # Plot each of the efficancies
        [ plt.plot( ff_wvl, np.abs( self.diff_efficacy[i,:]), line_style[i], label=leg[i] ) for i in range( len( self.diff_efficacy[:,0] ) ) ]
        plt.plot( ff_wvl, np.abs( combine ), line_style[-1], label=leg[-1] )
        if plot_total == True:
            total = np.array( [ np.abs( self.diff_efficacy[i,:]) for i in range( len( self.diff_efficacy[:,0] ) ) ] )
            total = np.sum( total, axis=0)
            plt.plot( ff_wvl, total, '-g', label='$\Sigma_i \eta_i$' )
        plt.legend(frameon=False)#, loc='left center' )
        plt.xlabel("Wavelength (nm)", size="x-large")
        plt.xticks(fontsize='large')
        plt.ylabel("Efficiencies", size="x-large")
        plt.yticks(fontsize='large')
        plt.tick_params( direction='in', axis='both', length=4, right=True, top=True, which="both" )
        plt.minorticks_on()
        plt.tick_params( which='minor', length=2, direction='in'  )
        plt.xlim( ff_wvl[-1], ff_wvl[0] )
        plt.ylim( 0, 1.50 )

        if fname == None:
            plt.show()
        else:
            plt.savefig( fname, dpi=dpi )

        return 0

    # Plots a colour map of |E|^2 with wavelength on the x-axis and angle on the y-axis
    def plot_wavelength_pattern( self,fname=None, trace_braggs=False, dpi=150, **kwarg ):
        # Get the frequncies
        ff_frq = mp.get_flux_freqs( self.n2f_obj )

        # Convert to wavelengths
        ff_wvl = np.divide( 1, ff_frq )

        # Get the field values and normilse
        field = np.abs( self.ff_data['Ez'] )**2
        for i in range( len( ff_wvl ) ):
            field[:,i] = field[:,i]/np.max(field[:,i])

        # Plot the field relitve to angles and wavelngths
        fig, axs = plt.subplots(  )
        axs.pcolormesh( ff_wvl, self.ff_angles, field, cmap='Blues',shading='flat' )

        if trace_braggs == True:
            x =np.linspace( ff_wvl[0], ff_wvl[-1], 1000 )
            [ axs.plot( x, np.arcsin( i*x/self.period ), '-r' ) for i in [-1,0,1] ]

        axs.tick_params( direction="in" )
        axs.set_xlabel( "Wavelength ($\mu m$)", size="x-large")
        axs.set_ylabel( "Angle (rad)", size="x-large") 
        axs.set_ylim( self.ff_angles[0], self.ff_angles[-1] )

        if fname == None:
            plt.show()
        else:
            plt.savefig( fname, dpi=dpi )

    
    # Allows simulation data to be saved to a JSON file
    def save_data( self, fname=None, **kwarg ):
        # Check if the file name is a string
        if fname == None or type( fname ) != str:
            raise TypeError( "'fname' must be a string" )
        # Attempt to open, else error
        try:
            data_file = open( fname + str('.json'), "w" )
        except:
            print("Could not save file")
            return 1

        # Load the far fied data
        if self.ff_data == None:
            output_data = { 'Ex': '[]',
                            'Ey': '[]',
                            'Ez': '[]',
                            'Hx': '[]',
                            'Hy': '[]',
                            'Hz': '[]',
                            'angles': '[]',
                            'points': '[]'
                          }
        else:
            output_data = self.ff_data
            output_data.update( { "angles": self.ff_angles,"points": self.ff_points } )

            # Convert the far field arrays to strings to preserve complex numbers
            output_data = { k:str(n.tolist()) for k,n in output_data.items() }

        # Add all the other data
        output_data.update( {
            "num_period": self.num_period,
            "period": self.period, 
            "grating_width": self.grating_width,
            "grating_height": self.grating_height,
            "wvl": self.wvl,
            "nwvl": self.nwvl,
            "dwvl": self.dwvl,
            "res": self.res,
            "run_2D": int( self.run_2D ),
            "diff_efficacy": self.diff_efficacy.tolist(),
            "incidence_flux_data": self.incidence_flux_data.tolist(),
            "ff_dist": self.ff_dist,
            "frq_values": self.frq_values.tolist(),
            "incidence_flux_data": self.incidence_flux_data.tolist(),
            "net_loss_flux_data": self.net_loss_flux_data.tolist() } )

        # Dump the data to the json file
        json_data = json.dumps( output_data )
        data_file.write( json_data )
        data_file.close()

        del( output_data )
        del( json_data )
        del( data_file )

        return 0 

    def save_n2f_obj( self, fname, **kwarg ):
        self.sim.save_near2far( fname, self.n2f_obj )
        return 0

    # Plots the far field
    def plot_far_field( self, x_axis="angle", wvl=None, fname=None, dpi=300, **kwarg ):
        ff_p_vector = []

        if wvl == None:
            wvl = self.wvl
        index = np.abs( self.frq_values - 1/wvl ).argmin()

        Px = np.real( np.conj( self.ff_data['Ey'][:,index] )*self.ff_data['Hz'][:,index] - np.conj( self.ff_data['Ez'][:,index] )*self.ff_data['Hy'][:,index] ) 
        Py = np.real( np.conj( self.ff_data['Ez'][:,index] )*self.ff_data['Hx'][:,index] - np.conj( self.ff_data['Ex'][:,index] )*self.ff_data['Hz'][:,index] ) 
        Pz = np.real( np.conj( self.ff_data['Ex'][:,index] )*self.ff_data['Hy'][:,index] - np.conj( self.ff_data['Ey'][:,index] )*self.ff_data['Hx'][:,index] ) 

        Pv = np.sqrt( Px**2 + Py**2 + Pz**2 )
        Pv_norm_abs = np.abs( Pv )**2
        Pv_norm_abs /= np.max( Pv_norm_abs )

        fig, axs = plt.subplots()
        if x_axis == "angle":
            axs.plot( self.ff_angles,  Pv_norm_abs , '-k' )
        else:
            axs.plot( self.ff_points,  Pv_norm_abs , '-k' )
        axs.tick_params( direction="in" )
        axs.grid( which="both" )
        axs.set_ylabel( "Poynting vector", size="x-large")
        axs.set_xlabel( "Angle (rad)", size="x-large") if x_axis == "angle" else axs.set_xlabel( "Far field position (μm)", size="x-large")

        # Save the file unless show plot if it is named, else just show
        plt.tight_layout()
        if fname == None:
            plt.show()
        else:
            plt.savefig( fname, dpi=dpi )

        return 0

    # Determines the efficny of the GMOT in regards to it's flux (power)
    def gmot_efficacy( self, **kwarg ):
        # Find the frequncy location for the wavelength

        # fix, (a-b)/a not b/a
        return 1-abs( self.net_loss_flux_data/self.incidence_flux_data )

