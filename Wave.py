import numpy as np
from dataclasses import dataclass, field
from scipy.optimize import bisect
import matplotlib.pyplot as plt
from Plot import Plot
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline

# refine solver for shearwaves
# wave structure
@dataclass
class Wave:  
    density = 2700            
    youngs_modulus = 68.9e9  # Young's modulus in Pascal
    poissons_ratio = 0.33  # Poisson's ratio (unitless)
    longitudinal_wave_velocity = np.sqrt(youngs_modulus*(1-poissons_ratio) / (density*(1+poissons_ratio)*(1-2*poissons_ratio))) # Longitudinal wave velocity of the material, in m/s r
    shear_wave_velocity = np.sqrt(youngs_modulus / (2*density*(1+poissons_ratio))) # Shear wave velocity of the material, in m/s r
    rayleigh_wave_velocity = shear_wave_velocity * ((0.862+1.14*poissons_ratio) / (1+poissons_ratio)) # Rayleigh wave velocity of the material, in m/s r  
    thickness = 1e-2 # plate thickness, in m
    half_thickness = thickness / 2 # Half thickness, in m
    freq_thickness_max = 10000 # Maximum value of frequency x thickness 
    freq_thickness_points = 100 # Number of frequency x thickness points 
    modes_functions : dict = field(init=False)
    velocities_dict : dict = field(init=False)
    increasing_mode : str = field(init=False)
    cp_step = 100
    cp_max = 15000

    def __post_init__(self):
        self.modes_functions = {
        'symmetric': self.calculate_symmetric,
        'antisymmetric': self.calculate_antisymmetric }
        self.velocities_dict = { 'C_R' : self.rayleigh_wave_velocity,
                                 'C_S' : self.shear_wave_velocity, 
                                 'C_L' : self.longitudinal_wave_velocity } 

    def calculate_dispersion_components(self, phase_velocity, freq_thickness):
        angular_freq = 2 * np.pi * (freq_thickness / self.thickness) 
        k = angular_freq / phase_velocity
        p = np.sqrt((angular_freq / self.longitudinal_wave_velocity)**2 - k**2, dtype=np.complex128)
        q = np.sqrt((angular_freq / self.shear_wave_velocity)**2 - k**2, dtype=np.complex128)
        return k, p, q

    def calculate_symmetric(self):
        pass

    def calculate_antisymmetric(self):
        pass

    def interpolate_result(self, result, kind='cubic'):
        # Function to interpolate result and update the dictionary
        interp_result = {}
        for key, values in result.items():

            if len(values) > 3:
                fd, cp = np.array(values).T
                # Perform interpolation using cubic spline
                interp_func = interp1d(fd, cp, kind=kind)

                # Generating finer fd values for smoother plot
                fd_finer = np.linspace(min(fd), max(fd), 1000)
                cp_interp = interp_func(fd_finer)

                # Store interpolated data in the new dictionary
                interp_result[key] = list(zip(fd_finer, cp_interp))

        return interp_result

    def calculate_group_wavenumber(self, result):
        # Function to calculate group velocity and wave number and update the dictionary
        for key, values in result.items():
            fd_val = [point[0] for point in values]
            cp_val = [point[1] for point in values]
            k_val, cg_val = [], []
            # Estimating the derivative of cp_val with respect to fd_val using interpolation
            univ_s = InterpolatedUnivariateSpline(fd_val, cp_val)
            cp_prime = univ_s.derivative()
            # Calculate cg using the estimated derivative
            for i in range(len(fd_val)):
                cg = cp_val[i]**2 * (cp_val[i] - (fd_val[i]) * cp_prime(fd_val[i]))**-1
                k = (fd_val[i]*2*np.pi/self.thickness)/cp_val[i]
                cg_val.append(cg)
                k_val.append(k)
            # Add cg_val as the third value in each tuple
            result[key] = [(fd_val[i], cp_val[i], cg_val[i], k_val[i]) for i in range(len(fd_val))]
        # return update dict
        return result

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

                # Check if sign changes
                if not (np.isnan(cp_0) or np.isnan(cp_1)) and np.sign(cp_0) != np.sign(cp_1):
                    # Using bisection method to find the root
                    cp_root = bisect(f=modes_func, a=init_cp, b=init_cp2, args=(fd,))
                    if np.abs(modes_func(cp_root, fd)) < 1e-3:
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

@dataclass
class Shearwave(Wave):
    def __init__(self):
        super().__init__()
        self.increasing_mode = "S_0"
    def calculate_symmetric(self, phase_velocity, freq_thickness):
        _,_,q = self.calculate_dispersion_components(phase_velocity, freq_thickness)
        return np.real(np.sin(q*self.half_thickness))

    def calculate_antisymmetric(self, phase_velocity, freq_thickness):
        _,_,q = self.calculate_dispersion_components(phase_velocity, freq_thickness)
        return  np.real(np.cos(q*self.half_thickness))

@dataclass
class Lambwave(Wave):
    def __init__(self):
        super().__init__()
        self.increasing_mode = "A_0"
    
    def calculate_symmetric(self, phase_velocity, freq_thickness):
        k,p,q = self.calculate_dispersion_components(phase_velocity, freq_thickness)
        return np.real(np.tan(q*self.half_thickness)/q + (4*(k**2)*p*np.tan(p*self.half_thickness))/(q**2 - k**2)**2)

    def calculate_antisymmetric(self, phase_velocity, freq_thickness):
        k,p,q = self.calculate_dispersion_components(phase_velocity, freq_thickness)
        return np.real(q * np.tan(q*self.half_thickness) + (((q**2 - k**2)**2)*np.tan(p*self.half_thickness))/(4*(k**2)*p))


if __name__ == "__main__":
    lamb = Lambwave()

    # Solve for symmetric and antisymmetric modes
    resultsym = lamb.solve_freq_equations('symmetric')
    resultanti = lamb.solve_freq_equations('antisymmetric')
    plot = Plot(lamb.freq_thickness_max, resultsym, resultanti, lamb.velocities_dict)
    plot.show_plots()
    