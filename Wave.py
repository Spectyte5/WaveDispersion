from array import array
import numpy as np
from dataclasses import dataclass, field
from scipy.optimize import bisect
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from Material import Material

@dataclass
class Wave: 
    material : Material
    modes_nums = {
        'symmetric': None,
        'antisymmetric': None
    }
    freq_thickness_max : int # Maximum value of frequency x thickness 
    freq_thickness_points : int # Number of frequency x thickness points 
    cp_step : int 
    cp_max : int 
    structure_mode : str 
    structure_freq : array
    rows : int 
    columns : int
    modes_functions : dict = field(init=False)
    velocities_dict : dict = field(init=False)
    increasing_mode : str = field(init=False)
    tolerance = 1e-3
    structure_cp = {}
    velocites_symmetric = None 
    velocites_antisymmetric = None
    structure_result = None 
    get_converted_mode = lambda _, key: int(key[2:]) * 2 if key.startswith('S') else int(key[2:]) * 2 + 1

    def __post_init__(self):
        self.modes_functions = {
        'symmetric': self.calculate_symmetric,
        'antisymmetric': self.calculate_antisymmetric }
        self.velocities_dict = { 'C_R' : self.material.rayleigh_wave_velocity,
                                 'C_S' : self.material.shear_wave_velocity, 
                                 'C_L' : self.material.longitudinal_wave_velocity } 
        
        # Solve for symmetric and antisymmetric modes
        self.velocites_symmetric = self.solve_freq_equations('symmetric')
        self.velocites_antisymmetric = self.solve_freq_equations('antisymmetric')

        # Calculate wave structure
        self.structure_result = self.calculate_wave_structure()

    def calculate_dispersion_components(self, phase_velocity, freq_thickness):
        angular_freq = 2 * np.pi * (freq_thickness / self.material.thickness) 
        k = angular_freq / phase_velocity
        p = np.sqrt((angular_freq / self.material.longitudinal_wave_velocity)**2 - k**2, dtype=np.complex128)
        q = np.sqrt((angular_freq / self.material.shear_wave_velocity)**2 - k**2, dtype=np.complex128)
        return k, p, q

    def interpolate_result(self, result, kind='cubic'):
        # Function to interpolate result and update the dictionary
        interp_result = {}
        for key, values in result.items():

            if len(values) > 3:
                fd, cp = np.array(values).T
                # Perform interpolation using cubic spline
                interp_func = interp1d(fd, cp, kind=kind)

                # Generating finer fd values for smoother plot
                fd_finer = np.linspace(min(fd), max(fd), 2000)
                cp_interp = interp_func(fd_finer)

                # Store interpolated data in the new dictionary
                interp_result[key] = list(zip(fd_finer, cp_interp))

                # Store cp as function of fd for wave structure
                self.structure_cp[key] = interp_func

        return interp_result

    def calculate_group_wavenumber(self, result):
        # Function to calculate group velocity and wave number and update the dictionary
        for key, values in result.items():
            fd_val = [point[0] for point in values]
            cp_val = [point[1] for point in values]
            k_val, cg_val = [], []
            # Estimating the derivative of cp_val with respect to fd_val using interpolation
            univ_s = InterpolatedUnivariateSpline(fd_val, cp_val)
            cp_prime, n = univ_s.derivative(), int(key[-1])
            # Calculate cg using the estimated derivative
            for i in range(len(fd_val)):
                cg = self.get_group_velocity_eq(cp_val[i], fd_val[i], cp_prime, key)
                k = (fd_val[i]*2*np.pi/self.material.thickness)/cp_val[i]
                cg_val.append(cg)
                k_val.append(k)
            # Add cg_val as the third value in each tuple
            result[key] = [(fd_val[i], cp_val[i], cg_val[i], k_val[i]) for i in range(len(fd_val))]
        # return updated dict
        return result

    def calculate_wave_structure(self, samples_x=100):

        fd_values, cp_values,  u_w_array = [], [], {}

        result = self.velocites_symmetric if self.structure_mode.startswith('S') else self.velocites_antisymmetric

        for _, values in result.items():
            for tuple_value in values:
                fd_values.append(tuple_value[0])
                cp_values.append(tuple_value[1])

        # Create array between -d/2 and d/2
        x = np.linspace(-self.material.half_thickness, self.material.half_thickness, samples_x) 

        for fd in self.structure_freq:
            cp = self.structure_cp[self.structure_mode](fd)          
            u, w = self.calculate_wavestructure_components(x, cp, fd)
            if fd not in u_w_array:
                u_w_array[fd] = []
            # Append result to array
            u_w_array[fd].append([u, w, x])
            
        return u_w_array

@dataclass
class Shearwave(Wave):
    def __init__(self, material, modes_nums, freq_thickness_max, freq_thickness_points, cp_step, cp_max, structure_mode: str, structure_freq: array, rows: int, columns: int):
        self.increasing_mode = 'S_0'
        self.modes_nums['symmetric'], self.modes_nums['antisymmetric'] = modes_nums
        super().__init__(material, freq_thickness_max, freq_thickness_points, cp_step, cp_max, structure_mode, structure_freq, rows, columns)

    def calculate_symmetric(self, phase_velocity, freq_thickness):
        _,_,q = self.calculate_dispersion_components(phase_velocity, freq_thickness)
        return np.real(np.sin(q*self.material.half_thickness))

    def calculate_antisymmetric(self, phase_velocity, freq_thickness):
        _,_,q = self.calculate_dispersion_components(phase_velocity, freq_thickness)
        return  np.real(np.cos(q*self.material.half_thickness))

    def get_group_velocity_eq(self, cp, fd, cp_prime, key):
        n = self.get_converted_mode(key)
        cg = self.material.shear_wave_velocity * np.sqrt(1 - ((n/2)**2) / ((fd/self.material.shear_wave_velocity)**2))
        return cg

    def calculate_wavestructure_components(self, x, cp, fd):
        k,_,_ = self.calculate_dispersion_components(cp, fd)

        n = self.get_converted_mode(self.structure_mode)

        if self.structure_mode.startswith('S'):
            B = 1  
            u = B * np.cos(n * np.pi * x / self.material.half_thickness) * np.exp(-1j * k * x)
            w = None
        else:
            A = 1
            u = A * np.cos(n * np.pi * x / self.material.half_thickness) * np.exp(-1j * k * x)
            w = None

        return u, w

    def solve_freq_equations(self, mode_type):
        fd_range = np.linspace(0, self.freq_thickness_max, self.freq_thickness_points)
        result = {}
        prefix = mode_type[0].upper() + '_'
        mode = 0

        for mode in range(0, self.modes_nums[mode_type]):
            for fd in fd_range:                       
                key = prefix + str(mode)
                n = self.get_converted_mode(key)
                # calculate phase velocity [m/s]
                cp = 2 * self.material.shear_wave_velocity * fd / np.sqrt(4 * fd**2 - (n * self.material.shear_wave_velocity)**2)
                if not np.isnan(cp) and cp < self.cp_max:
                    if key not in result:
                        result[key] = []
                    result[key].append([fd, cp])
        
        # Interpolate for smoother plot
        result = self.interpolate_result(result)

        # Calculate group velocity and wave number
        result = self.calculate_group_wavenumber(result)
        return result

@dataclass
class Lambwave(Wave):
    def __init__(self, material, modes_nums, freq_thickness_max, freq_thickness_points, cp_step, cp_max, structure_mode: str, structure_freq: array, rows: int, columns: int):
        self.increasing_mode = 'A_0'
        self.modes_nums['symmetric'], self.modes_nums['antisymmetric'] = modes_nums
        super().__init__(material, freq_thickness_max, freq_thickness_points, cp_step, cp_max, structure_mode, structure_freq, rows, columns)
    
    def calculate_symmetric(self, phase_velocity, freq_thickness):
        k,p,q = self.calculate_dispersion_components(phase_velocity, freq_thickness)
        return np.real(np.tan(q*self.material.half_thickness)/q + (4*(k**2)*p*np.tan(p*self.material.half_thickness))/(q**2 - k**2)**2)

    def calculate_antisymmetric(self, phase_velocity, freq_thickness):
        k,p,q = self.calculate_dispersion_components(phase_velocity, freq_thickness)
        return np.real(q * np.tan(q*self.material.half_thickness) + (((q**2 - k**2)**2)*np.tan(p*self.material.half_thickness))/(4*(k**2)*p))

    def get_group_velocity_eq(self, cp, fd, cp_prime, key):
        cg = cp**2 * (cp - fd * cp_prime(fd))**-1
        return cg

    def calculate_wavestructure_components(self, x, cp, fd):
        k,p,q = self.calculate_dispersion_components(cp, fd)

        if self.structure_mode.startswith('S'):
            A, B = -2*k*q*np.cos(q*self.material.half_thickness) / ((k**2 - q**2) * np.cos(p*self.material.half_thickness)), 1   
            u = 1j*(k*A*np.cos(p*x) + q*B*np.cos(q*x))
            w = -p*A*np.sin(p*x) + k*B*np.sin(q*x)
        else:
            A, B = 2*k*q*np.sin(q*self.material.half_thickness) / ((k**2 - q**2) * np.sin(p*self.material.half_thickness)), 1
            u = 1j*(k*A*np.sin(p*x) - q*B*np.sin(q*x))
            w = p*A*np.cos(p*x) + k*B*np.cos(q*x)

        return u, w

    def solve_freq_equations(self, mode_type):
        fd_range = np.linspace(0, self.freq_thickness_max, self.freq_thickness_points)
        result = {}
        modes_func = self.modes_functions.get(mode_type)
        prefix = mode_type[0].upper() + '_'
        initial_fd = []

        for fd in fd_range:
            # Initial phase velocity estimation
            mode = 0
            init_cp = 0
            init_cp2 = init_cp + self.cp_step

            while init_cp2 < self.cp_max:
                # Evaluate the function values
                cp_0 = modes_func(init_cp, fd)
                cp_1 = modes_func(init_cp2, fd)

                # Check if in the num_of modes range:
                if mode < self.modes_nums[mode_type]:
                    # Check if sign changes
                    if not (np.isnan(cp_0) or np.isnan(cp_1)) and np.sign(cp_0) != np.sign(cp_1):
                        # Using bisection method to find the root
                        cp_root = bisect(f=modes_func, a=init_cp, b=init_cp2, args=(fd,))
                        if np.abs(modes_func(cp_root, fd)) < self.tolerance:
                            key = prefix + str(mode)

                            # if key is not in result, create new list
                            if key not in result:
                                result[key] = []

                            # check if value increased for modes different then increasing mode
                            while key != self.increasing_mode and result[key] and cp_root > result[key][-1][1]:
                                mode += 1
                                key = prefix + str(mode)
                                if key not in result:
                                    result[key] = []

                            if not result[key]:
                                if fd in initial_fd:
                                    continue  # Skip if fd exists in initial_fd
                                initial_fd.append(fd)

                            result[key].append([fd, cp_root])
                            mode += 1

                # Update init_cps for the next iteration
                init_cp += self.cp_step
                init_cp2 += self.cp_step
        
        # Interpolate for smoother plot
        result = self.interpolate_result(result)

        # Calculate group velocity and wave number
        result = self.calculate_group_wavenumber(result)
        return result