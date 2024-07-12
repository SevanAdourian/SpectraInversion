# Imports
import numpy as np
from scipy.interpolate import interp1d

#
# ---
#

_header_fields = ['name',
                  'if_anis', 't_ref', 'if_deck',
                  'num_layers', 'index_icb', 'index_cmb', 'num_ocean_layers']

parameter_column_map = {'r': 0, 'rho': 1, 'vpv': 2, 'vsv': 3,
                        'qk': 4, 'qm': 5, 'vph': 6, 'vsh': 7,
                        'eta': 8}

class Model1D:
    """
    A python class for manipulating MINOS / yannos style tabular models
    """

    def __init__(self, model_file_name = []):
        # set the model file name
        self.file_name = model_file_name
        
        # initialize the list of interpolators 
        self.interpolate = []

        # set the populated flag to False
        self.populated = False

        # export the header field names
        self.header_fields = _header_fields

    def get_file_name(self):
        # return the model file name
        return self.file_name

    def set_file_name(self, model_file_name):
        # set the model file name
        self.file_name = model_file_name

    def get_header(self):
        # check for model population
        if not self.populated:
            raise Exception('Call to get_header() when the model has not been populated')

        # return the header fields as a tuple
        return (self.name,
                self.if_anis, self.t_ref, self.if_deck,
                self.num_layers, self.index_icb, self.index_cmb, self.num_ocean_layers)

    def set_header(self, header_tuple):
        # make sure this is a tabular model
        if header_tuple[_header_fields.index('if_deck')] == 0:
            raise Exception('Non-tabular model detected in set_model()')
        
        # populate the header
        (self.name,
         self.if_anis, self.t_ref, self.if_deck,
         self.num_layers, self.index_icb, self.index_cmb, self.num_ocean_layers) = header_tuple

    def get_table(self):
        # check for model population
        if not self.populated:
            raise Exception('Call to get_table() when model has not been populated')

        # return the model table
        return self.table

    def set_table(self, model_table):
        # populate the model table
        self.table = model_table

        # find discontinuity indices
        self.index_discon = np.array([ l for l in range(self.num_layers - 1)
                                       if self.table[l,0] == self.table[l+1,0] ])
        self.num_discon = self.index_discon.size

        # save discontinuity radii as well
        self.radius_discon = self.table[self.index_discon,0]

        # set the populated flag to True
        self.populated = True

    def get_parameter(self, parameter_name):
        # check for model population
        if not self.populated:
            raise Exception('Call to get_parameter() when model has not been populated')

        # return the column of the model table
        return self.table[:, parameter_column_map[parameter_name]]

    def set_parameter(self, parameter_name, values):
        # check for model population
        if not self.populated:
            raise Exception('Call to get_parameter() when model has not been populated')

        # set column of the table
        self.table[:, parameter_column_map[parameter_name]] = values

    def get_discons(self):
        # check for model population
        if not self.populated:
            raise Exception('Call to get_discons() when model has not been populated')

        # return the model table
        return self.radius_discon, self.index_discon

    def set_model(self, header_tuple, model_table):
        # populate the header
        (self.name,
         self.if_anis, self.t_ref, self.if_deck,
         self.num_layers, self.index_icb, self.index_cmb, self.num_ocean_layers) = header_tuple

        # make sure this is a tabular model
        if self.if_deck == 0:
            raise Exception('Non-tabular model detected in set_model()')

        # populate the model table
        self.table = model_table

        # find discontinuity indices
        self.index_discon = np.array([ l for l in range(self.num_layers - 1)
                                       if self.table[l,0] == self.table[l+1,0] ])
        self.num_discon = self.index_discon.size

        # save discontinuity radii as well
        self.radius_discon = self.table[self.index_discon,0]

        # set the populated flag to True
        self.populated = True

    def load_from_file(self):
        # check for initialization of file name
        if not self.file_name:
            raise Exception('Call to load_from_file() when file_name not initialized')

        # open the previously specified model file
        f_in = open(self.file_name, "r");

        # begin with the header ...
        
        #  read in the model name
        self.name = f_in.readline().strip()

        #  read in if_anis, t_ref, and if_deck
        line = f_in.readline().strip().split()
        self.if_anis = int(line[0])
        self.t_ref = float(line[1])
        self.if_deck = int(line[2])

        #  read in number of layers in the model and the ocean only, as
        #  well as the indices (from below) of the ICB and CMB discons
        line = f_in.readline().strip().split()
        self.num_layers = int(line[0])
        self.index_icb = int(line[1]) - 1
        self.index_cmb = int(line[2]) - 1
        self.num_ocean_layers = int(line[3])

        # read in the remainder of the file (the table) as numpy ndarray
        self.table = np.loadtxt(f_in)

        # close the file
        f_in.close()

        # find discontinuity indices (again from below) 
        self.index_discon = np.array([ l for l in range(self.num_layers - 1)
                                       if self.table[l,0] == self.table[l+1,0] ])
        self.num_discon = self.index_discon.size

        # save discontinuity radii as well
        self.radius_discon = self.table[self.index_discon,0]

        # set the populated flag to True
        self.populated = True

    def write_to_file(self):
        # check for initialization of file name
        if not self.file_name:
            raise Exception('Call to write_to_file() when file_name not initialized')
        
        # check for model population
        if not self.populated:
            raise Exception('Call to write_to_file() when model has not been populated')

        # open the output file
        f_out = open(self.file_name, "w")

        # write the header
        f_out.write('%s\n' % (self.name))
        f_out.write('%i %f %i\n' % (self.if_anis, self.t_ref, self.if_deck))
        f_out.write('%i %i %i %i\n' % (self.num_layers, self.index_icb + 1, self.index_cmb + 1, self.num_ocean_layers))

        # write the model table
        for l in range(self.num_layers):
            f_out.write('%8.0f%9.2f%9.2f%9.2f%9.1f%9.1f%9.2f%9.2f%9.5f\n' % tuple(self.table[l,:]))

        # close the file
        f_out.close()

    def __init_interpolation(self):
        # initialize interpolators for each continuous model segment ...
        #  store lower bounds on each segment
        self.radius_segment = np.zeros((self.num_discon+1,))
        self.radius_segment[0] = np.array([ self.table[0,0] ])
        self.radius_segment[1:] = self.radius_discon

        #  segment 0: center to first discon
        lower_bound = 0
        upper_bound = self.index_discon[0]
        self.interpolate = [ interp1d(self.table[lower_bound:upper_bound+1,0].T,
                                      self.table[lower_bound:upper_bound+1,1:].T) ]

        #  segment i: discon i to i+1
        for d in range(self.index_discon.size - 1):
            lower_bound = self.index_discon[d]+1
            upper_bound = self.index_discon[d+1]
            self.interpolate.append(interp1d(self.table[lower_bound:upper_bound+1,0].T,
                                             self.table[lower_bound:upper_bound+1,1:].T))

        #  segment num_discon: final discon to the surface (this will likely be the ocean)
        lower_bound = self.index_discon[-1]
        upper_bound = self.num_layers - 1
        self.interpolate.append(interp1d(self.table[lower_bound:upper_bound+1,0].T,
                                         self.table[lower_bound:upper_bound+1,1:].T))


    def get_values(self, radii, eval_discon_above = [], parameter = None):
        # if interpolation has not already been initialized, do so
        if not self.interpolate:
            self.__init_interpolation()

        # evaluate from below by default
        if len(eval_discon_above) == 0:
            eval_discon_above = np.array([ False for ir in range(radii.size) ])
            
        # initialize the result ndarray to zero
        result = np.zeros((radii.size,9))
        
        # loop over radii ...
        for ir in range(radii.size):
            # test whether the radius is a discontinuity
            if radii[ir] in self.radius_discon:
                # if so, evaluate from below unless eval_discon_above is True
                if eval_discon_above[ir]:
                    result[ir,:] = self.table[self.index_discon[radii[ir] == self.radius_discon] + 1,:]
                else:
                    result[ir,:] = self.table[self.index_discon[radii[ir] == self.radius_discon],:]
            else:
                # interpolate the model table fields to the desired radius
                if radii[ir] == 0.0:
                    below = [0]
                else:
                    below, = np.nonzero(self.radius_segment < radii[ir])
                result[ir,0] = radii[ir]
                result[ir,1:] = self.interpolate[below[-1]](radii[ir])

        # possibly strip the result down to a single parameter
        if parameter:
            result = result[:, parameter_column_map[parameter]]

        # return the assembled result
        return result

    
    def perturb_layer(self, parameter_name, start_layer, end_layer, pert_value):
        param = self.get_parameter(parameter_name)
        param[start_layer:end_layer] *= (1+pert_value)
        self.set_parameter(parameter_name, param)

        return
