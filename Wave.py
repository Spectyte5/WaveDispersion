from array import array
import numpy as np
from dataclasses import dataclass, field
from scipy.optimize import bisect
from scipy.interpolate import UnivariateSpline, InterpolatedUnivariateSpline
from Material import Material
from scipy.special import jv, iv, kv, yv

@dataclass
class Wave: 
    """
    A class representing a base wave type.

    Dataclass storing parameters connected to general dispersion of the waves,
    like max value of phase velocity and frequency times thickness,step between the values and
    mode amount for symmetric and antisymmetric modes. It also has all the dictionaries for
    storing the results, for velocity and wave structure.


    Attributes:
    material (Material) : Material class object providing the material information.
    modes_nums (dict) : Number of symmetric and antisymmetric modes to find.
    freq_thickness_max (int) : Max value of Frequency x Thickness [kHz x mm].
    cp_max (int) : Max value of Frequency x Thickness [m/s].
    structure_mode (str) : Which mode should be used for wavestructure plot.
    structure_freq (array) : Frequencies at which to check Wavestructure.
    rows (int) : Number of rows for Wavestructure plot.
    colsumns (int) : Number of columns for Wavestructure plot.
    freq_thickness_points (int) : Number of frequency x thickness points to find.
    cp_step (int) : Step between phase velocity points checked.
    kind (int) : Kind of interpolation, 3 is used for cubic.
    smoothen (int) : Order of smoothening used for interpolation.
    tolerance (float) : Tolerance for the root-finding.
    modes_functions (dict) : Maps mode_type with functions calculating symmetric and antisymmetric modes.
    velocities_dict (dict) : Maps material bulk velocities.
    increasing_mode (str) : Unique mode that increases instead of decreasing with frequency x thickness.
    structure_cp (dict) : Phase velocity used for wavestructure calculation.
    structure_result (dict) : Result obtained for wavestructure.

    Methods:
        __post_init__():
            Post-initialization method that sets default step values and plate velocites values.
        calculate_dispersion_components(phase_velocity, freq_thickness):
            Calculates dispersion componetns k, p and q.
        interpolate_result(result):
            Interpolate result and update the dictionary
        calculate_group_wavenumber(result):
            Calculates group velocity and wave number and updates the dictionary
        calculate_wave_structure(samples_x):
            Calculates wavestructure components u and w.
    """
    material : Material
    modes_nums = {
        'symmetric': None,
        'antisymmetric': None
    }
    freq_thickness_max : int 
    cp_max : int 
    structure_mode : str 
    structure_freq : array
    rows : int 
    columns : int
    freq_thickness_points : int
    cp_step : int
    kind = 3
    smoothen = 0
    tolerance = 1e-3
    modes_functions : dict = field(init=False)
    velocities_dict : dict = field(init=False)
    increasing_mode : str = field(init=False)
    structure_cp = {}
    structure_result = None 
    set_converted_mode = lambda x: f"{'S_' if int(x.split('_')[1]) % 2 == 0 else 'A_'}{int(x.split('_')[1]) // 2}"
    get_converted_mode = lambda _, key: int(key[2:]) * 2 if key.startswith('S') else int(key[2:]) * 2 + 1

    def __post_init__(self):
        """
        Post-initialization method that sets default step values and plate velocites values.

        This method is automatically called after the dataclass instance is initialized. 
        It sets default step values and plate velocites values.
        """
        self.velocities_dict = { 'C_R' : self.material.rayleigh_wave_velocity,
                                 'C_S' : self.material.shear_wave_velocity, 
                                 'C_L' : self.material.longitudinal_wave_velocity } 
        
        if not self.cp_step:
            self.cp_step = self.cp_max // 100
        if not self.freq_thickness_points:
            self.freq_thickness_points = self.freq_thickness_max // 100

    def calculate_dispersion_components(self, phase_velocity: float, freq_thickness: float) -> float | float | float:
        """
        Calculates dispersion componetns k, p and q.

        Parameters:
            phase_velocity (float): Current phase velocity value.
            freq_thickness (float): Current frequency and thickness product value.

        Returns:
            float: wavenumber
            float: p-component
            float: q-component
        """
        angular_freq = 2 * np.pi * (freq_thickness / self.material.thickness) 
        k = angular_freq / phase_velocity
        p = np.sqrt((angular_freq / self.material.longitudinal_wave_velocity)**2 - k**2, dtype=np.complex128)
        q = np.sqrt((angular_freq / self.material.shear_wave_velocity)**2 - k**2, dtype=np.complex128)
        return k, p, q

    def interpolate_result(self, result: dict) -> dict :
        """
        Interpolate result and update the dictionary

        Performs interpolation using cubic spline, generates finer fd values for smoother plot 
        and stores interpolated data in the new dictionary. Additionaly cp is also stored
        as function of frequency and thickness product for wave structure.

        Parameters:
            result (dict) : Dictionary with results to interpolate.

        Returns:
            dict: interpolated result
        """
        interp_result = {}
        for key, values in result.items():

            if len(values) > 3:
                fd, cp = np.array(values).T
                interp_func = UnivariateSpline(fd, cp, k=self.kind, s=self.smoothen)
                fd_finer = np.linspace(min(fd), max(fd), 2000)
                cp_interp = interp_func(fd_finer)
                interp_result[key] = list(zip(fd_finer, cp_interp))
                self.structure_cp[key] = interp_func

        return interp_result

    def calculate_group_wavenumber(self, result: dict) -> dict:
        """
        Calculates group velocity and wave number and updates the dictionary

        Estimates the derivative of phase velocity with respect to fd using interpolation,
        then calculates group_velocity using the estimated derivative.
        Adds cg_val as the third value in each tuple and returns updated dict.

        Parameters:
            result (dict) : Dictionary with dispersion results.

        Returns:
            dict: Dictionary with extended results
        """
        for key, values in result.items():
            fd_val = [point[0] for point in values]
            cp_val = [point[1] for point in values]
            k_val, cg_val = [], []
            univ_s = InterpolatedUnivariateSpline(fd_val, cp_val)
            cp_prime = univ_s.derivative()
            for i in range(len(fd_val)):
                cg = self.get_group_velocity_eq(fd_val[i], key) if isinstance(self, Shearwave) \
                   else self.get_group_velocity_eq(cp_val[i], fd_val[i], cp_prime)
                k = (fd_val[i]*2*np.pi/self.material.thickness)/cp_val[i]
                cg_val.append(cg)
                k_val.append(k)
            result[key] = [(fd_val[i], cp_val[i], cg_val[i], k_val[i]) for i in range(len(fd_val))]

        return result

    def calculate_wave_structure(self, samples_x: int=100) -> dict:
        """
        Calculates wave_structure in and out of plane components.

        Creates array between -thickness/2 and thickness/2. 
        Uses the formula for wave structure to calculate u and w compontents.
        Depending on the cases appends results to array.

        Parameters:
            samples_x (int, optional) : Number of samples for the thickness array.

        Returns:
            list: List of in and out of plane components
        """
        fd_values, cp_values,  u_w_array = [], [], {}

        result = self.velocites_symmetric if self.structure_mode.startswith('S') else self.velocites_antisymmetric

        for _, values in result.items():
            for tuple_value in values:
                fd_values.append(tuple_value[0])
                cp_values.append(tuple_value[1])

        x = np.linspace(-self.material.half_thickness, self.material.half_thickness, samples_x) 

        for fd in self.structure_freq:
            cp = self.structure_cp[self.structure_mode](fd)          
            u, w = self.calculate_wavestructure_components(x) if isinstance(self, Shearwave) else \
                self.calculate_wavestructure_components(x, cp, fd)
            if fd not in u_w_array:
                u_w_array[fd] = []
            u_w_array[fd].append([u, w, x])
            
        return u_w_array

@dataclass
class Shearwave(Wave):
    """
    A class representing a Shearwave type.

    Dataclass inheriting Wave class, providing results and parameters specific to Shearwaves.
    Class has different methods for calculation of dispersion and wavestructure for Shearwaves.

    Attributes:
    set_converted_key (lambda) : sets key from Shearwave to Lambwave style key.
    get_converted_mode (lambda) : gets mode number for Shearwave from Lambwave style key.

    Methods:
    __init__(material, modes_nums, freq_thickness_max, cp_max, structure_mode, structure_freq, rows, columns, freq_thickness_points, cp_step):
        Initialization method that calls function calculating results.
    __post_init__():
        Post-initialization method that calls function calculating results.
    calculate_symmetric(phase_velocity, freq_thickness):
        Calcultes symmetric modes dispersion equation.
    calculate_antisymmetric(phase_velocity, freq_thickness):
        Calcultes antisymmetric modes dispersion equation.
    get_group_velocity_eq(cp, fd, cp_prime, key):
        Calculates group velocity from current phase velocity and fd.
    calculate_wavestructure_components(x, cp, fd):
        Calculates wavestructure components u and w.
    solve_freq_equations(mode_type):
        Solves frequency equations and returns function of cp with respect to fd.
    """
    set_converted_key = lambda _, x: f"{'S_' if int(x.split('_')[1]) % 2 == 0 else 'A_'}{int(x.split('_')[1]) // 2}"
    get_converted_mode = lambda _, key: int(key[2:]) * 2 if key.startswith('S') else int(key[2:]) * 2 + 1
    
    def __init__(self, material, modes_nums, freq_thickness_max, cp_max, structure_mode=None, structure_freq=None, rows=None, columns=None, freq_thickness_points=None, cp_step=None):
        """
        Initialization method that calls function calculating results.

        This method is automatically called when the dataclass instance is initialized. 
        It calls init method of the base class Wave and sets additional parameters.
        """
        self.increasing_mode = 'S_0'
        self.modes_nums['symmetric'], self.modes_nums['antisymmetric'] = modes_nums
        super().__init__(material, freq_thickness_max, cp_max, structure_mode, structure_freq, rows, columns, freq_thickness_points, cp_step)

    def __post_init__(self):
        """
        Post-initialization method that calls function calculating results.

        This method is automatically called after the dataclass instance is initialized. 
        It solves dispersion equation for symmetric and antisymmetric modes and calculates wavestructure.
        """
        self.modes_functions = {
            'symmetric': self.calculate_symmetric,
            'antisymmetric': self.calculate_antisymmetric }

        super().__post_init__()
        self.velocites_symmetric = self.solve_freq_equations('symmetric')
        self.velocites_antisymmetric = self.solve_freq_equations('antisymmetric')
        if self.structure_freq and self.structure_mode:
            self.structure_mode=self.set_converted_key(self.structure_mode)
            self.structure_result = self.calculate_wave_structure()

    def calculate_symmetric(self, phase_velocity, freq_thickness) -> float:
        """
        Calcultes symmetric modes dispersion equation.

        Uses formula for symmetric modes dispersion equation and returns its real component.

        Parameters:
            phase_velocity (float): Current phase velocity value.
            freq_thickness (float): Current frequency and thickness product value.

        Returns:
            float: Dispersion result.
        """
        _,_,q = self.calculate_dispersion_components(phase_velocity, freq_thickness)
        return np.real(np.sin(q*self.material.half_thickness))

    def calculate_antisymmetric(self, phase_velocity, freq_thickness):
        """
        Calcultes antisymmetric modes dispersion equation.

        Uses formula for antisymmetric modes dispersion equation and returns its real component.

        Parameters:
            phase_velocity (float): Current phase velocity value.
            freq_thickness (float): Current frequency and thickness product value.

        Returns:
            float: Dispersion result.
        """
        _,_,q = self.calculate_dispersion_components(phase_velocity, freq_thickness)
        return  np.real(np.cos(q*self.material.half_thickness))

    def get_group_velocity_eq(self, fd: float, key: str) -> float:
        """
        Calculates group velocity from current phase velocity and fd.

        Parameters:
            freq_thickness (float): Current frequency and thickness product value.
            key (str, not used): Full mode name

        Returns:
            float: group velocity
        """
        n = self.get_converted_mode(key)
        cg = self.material.shear_wave_velocity * np.sqrt(1 - ((n/2)**2) / ((fd/self.material.shear_wave_velocity)**2))
        return cg

    def calculate_wavestructure_components(self, x: np.ndarray) -> float | float:
        """
        Calculates wavestructure components u and w.

        Uses formula for wavestructure to obtain in-plane and out-of-plane components.

        Parameters:
            x (np.ndarray) : Array of points between -thickness/2 and thickness/2.

        Returns:
            float: u-component
            float: w-component
        """
        n = self.get_converted_mode(self.structure_mode)

        if self.structure_mode.startswith('S'):
            B = 1  
            u = B * np.cos(n * np.pi * x / self.material.thickness) 
            w = None
        else:
            A = 1
            u = A * np.sin(n * np.pi * x / self.material.thickness) 
            w = None

        return u, w

    def solve_freq_equations(self, mode_type: str) -> dict: 
        """
        Solves frequency equations and returns function of cp with respect to fd.
        The algorithm is based on the one presented in J.L.Roses Ultrasonic Guided Waves in Solid Media - chapter 12. 
        Use relation between material's shear_wave_velocity and frequency x thickness to obtain the phase velocity 

        Parameters:
            mode_type (str) : Types of the modes shown on the plot: sym, anti or both.

        Returns:
            dict: Dictionary with dispersion results.
        """
        fd_range = np.linspace(0, self.freq_thickness_max, self.freq_thickness_points)
        result = {}
        prefix = mode_type[0].upper() + '_'
        mode = 0

        for mode in range(0, self.modes_nums[mode_type]):
            for fd in fd_range:                       
                key = prefix + str(mode)
                n = self.get_converted_mode(key)
                cp = 2 * self.material.shear_wave_velocity * fd / np.sqrt(4 * fd**2 - (n * self.material.shear_wave_velocity)**2)
                if not np.isnan(cp) and cp < self.cp_max:
                    if key not in result:
                        result[key] = []
                    result[key].append([fd, cp])
        
        result = self.interpolate_result(result)

        result = self.calculate_group_wavenumber(result)
        return result

@dataclass
class Lambwave(Wave):
    """
    A class representing a Lambwave type.

    Dataclass inheriting Wave class, providing results and parameters specific to Lambwaves.
    Class has different methods for calculation of dispersion and wavestructure for Lambwaves.

    Methods:
    __init__(material, modes_nums, freq_thickness_max, cp_max, structure_mode, structure_freq, rows, columns, freq_thickness_points, cp_step):
        Initialization method that calls function calculating results.
    __post_init__():
        Post-initialization method that calls function calculating results.
    calculate_symmetric(phase_velocity, freq_thickness):
        Calcultes symmetric modes dispersion equation.
    calculate_antisymmetric(phase_velocity, freq_thickness):
        Calcultes antisymmetric modes dispersion equation.
    get_group_velocity_eq(cp, fd, cp_prime, key):
        Calculates group velocity from current phase velocity and fd.
    calculate_wavestructure_components(x, cp, fd):
        Calculates wavestructure components u and w.
    solve_freq_equations(mode_type):
        Solves frequency equations and returns function of cp with respect to fd.
    """
    def __init__(self, material, modes_nums, freq_thickness_max, cp_max, structure_mode=None, structure_freq=None, rows=None, columns=None, freq_thickness_points=None, cp_step=None):
        """
        Initialization method that calls function calculating results.

        This method is automatically called when the dataclass instance is initialized. 
        It calls init method of the base class Wave and sets additional parameters.
        """
        self.increasing_mode = 'A_0'
        self.modes_nums['symmetric'], self.modes_nums['antisymmetric'] = modes_nums
        super().__init__(material, freq_thickness_max, cp_max, structure_mode, structure_freq, rows, columns, freq_thickness_points, cp_step)
    
    def __post_init__(self):
        """
        Post-initialization method that calls function calculating results.

        This method is automatically called after the dataclass instance is initialized. 
        It solves dispersion equation for symmetric and antisymmetric modes and calculates wavestructure.
        """
        super().__post_init__()
        self.modes_functions = {
            'symmetric': self.calculate_symmetric,
            'antisymmetric': self.calculate_antisymmetric }
        self.velocites_symmetric = self.solve_freq_equations('symmetric')
        self.velocites_antisymmetric = self.solve_freq_equations('antisymmetric')
        if self.structure_freq and self.structure_mode:
            self.structure_result = self.calculate_wave_structure()

    def calculate_symmetric(self, phase_velocity: float, freq_thickness: float) -> float:
        """
        Calcultes symmetric modes dispersion equation.

        Uses formula for symmetric modes dispersion equation and returns its real component.

        Parameters:
            phase_velocity (float): Current phase velocity value.
            freq_thickness (float): Current frequency and thickness product value.

        Returns:
            float: Dispersion result.
        """
        k,p,q = self.calculate_dispersion_components(phase_velocity, freq_thickness)
        return np.real(np.tan(q*self.material.half_thickness)/q + (4*(k**2)*p*np.tan(p*self.material.half_thickness))/(q**2 - k**2)**2)

    def calculate_antisymmetric(self, phase_velocity: float, freq_thickness: float) -> float:
        """
        Calcultes antisymmetric modes dispersion equation.

        Uses formula for antisymmetric modes dispersion equation and returns its real component.

        Parameters:
            phase_velocity (float): Current phase velocity value.
            freq_thickness (float): Current frequency and thickness product value.

        Returns:
            float: Dispersion result.
        """
        k,p,q = self.calculate_dispersion_components(phase_velocity, freq_thickness)
        return np.real(q * np.tan(q*self.material.half_thickness) + (((q**2 - k**2)**2)*np.tan(p*self.material.half_thickness))/(4*(k**2)*p))

    def get_group_velocity_eq(self, cp: float, fd: float, cp_prime: float) -> float:
        """
        Calculates group velocity from current phase velocity and fd.

        Parameters:
            phase_velocity (float): Current phase velocity value.
            freq_thickness (float): Current frequency and thickness product value.
            cp_prime (float): Current phase velocity derivative.

        Returns:
            float: Group velocity
        """
        cg = cp**2 * (cp - fd * cp_prime(fd))**-1
        return cg

    def calculate_wavestructure_components(self, x: np.ndarray, cp: float, fd: float) -> float | float:
        """
        Calculates wavestructure components u and w.

        Uses formula for wavestructure to obtain in-plane and out-of-plane components.

        Parameters:
            x (np.ndarray) : Array of points between -thickness/2 and thickness/2.
            phase_velocity (float): Current phase velocity value.
            freq_thickness (float): Current frequency and thickness product value.

        Returns:
            float: u-component
            float: w-component
        """
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

    def solve_freq_equations(self, mode_type: str) -> dict:
        """
        Solves frequency equations and returns a function of cp with respect to fd.

        The algorithm is based on the one presented in J.L. Rose's "Ultrasonic Guided Waves in Solid Media," 
        Chapter 6. The algorithm proceeds as follows:
        
        1. Choose a frequency-thickness product fd_0.
        2. Make an initial estimate of the phase velocity cp_0.
        3. Evaluate the signs of the left-hand sides of the frequency equations.
        4. Choose another phase velocity cp_1 > cp_0 and re-evaluate the signs of the frequency equations.
        5. Repeat steps (3) and (4) until a sign change occurs, assuming this happens between cp_n and cp_n+1.
        6. Use bisection to precisely locate the phase velocity in the interval cp_n < cp < cp_n+1 
           where the left-hand side of the equation is close enough to zero.
        7. After finding the root, continue searching at this fd for other roots according to steps (2) through (6).
        8. Choose another fd product and repeat steps (2) through (7).

        Parameters:
            mode_type (str): The type of modes shown on the plot: 'sym', 'anti', or 'both'.

        Returns:
            dict: The results of the frequency equation solutions, indexed by mode type.
        """
        fd_range = np.linspace(0, self.freq_thickness_max, self.freq_thickness_points)
        result = {}
        modes_func = self.modes_functions.get(mode_type)
        prefix = mode_type[0].upper() + '_'
        initial_fd = []

        for fd in fd_range:
            mode = 0
            init_cp = 0
            init_cp2 = init_cp + self.cp_step

            while init_cp2 < self.cp_max:
                cp_0 = modes_func(init_cp, fd)
                cp_1 = modes_func(init_cp2, fd)

                if mode < self.modes_nums[mode_type]:
                    if not (np.isnan(cp_0) or np.isnan(cp_1)) and np.sign(cp_0) != np.sign(cp_1):
                        cp_root = bisect(f=modes_func, a=init_cp, b=init_cp2, args=(fd,))
                        if np.abs(modes_func(cp_root, fd)) < self.tolerance:
                            key = prefix + str(mode)
                            if key not in result:
                                result[key] = []

                            while key != self.increasing_mode and result[key] and cp_root > result[key][-1][1]:
                                mode += 1
                                key = prefix + str(mode)
                                if key not in result and mode:
                                    result[key] = []
                            
                            if not result[key]:
                                if fd in initial_fd:
                                    break
                                initial_fd.append(fd)

                            if mode < self.modes_nums[mode_type]:
                                result[key].append([fd, cp_root])
                                mode += 1
                            else:
                                result.pop(key, None)

                init_cp += self.cp_step
                init_cp2 += self.cp_step

        result = self.interpolate_result(result)

        result = self.calculate_group_wavenumber(result)
        return result

@dataclass
class Axialwave(Wave):
    """
    A class representing a Axialwave type.

    Dataclass inheriting Wave class, providing results and parameters specific to Axialwaves.
    Class has different methods for calculation of dispersion and wavestructure for Axialwaves.

    Methods:
    __init__(material, modes_nums, freq_thickness_max, cp_max, structure_mode, structure_freq, rows, columns, freq_thickness_points, cp_step):
        Initialization method that calls function calculating results.
    __post_init__():
        Post-initialization method that calls function calculating results.
    calculate_dispersion_components(phase_velocity, freq_thickness):
        Calculates dispersion components k, alpha_sq and beta_sq.
    calculate_bessel_functions(n, r, alpha_sq, beta_sq):
        Calculates bessel or modified bessel functions depending on the input.
    calculate_characteristic_equation(phase_velocity, freq_thickness, n):
        Calculates bessel or modified bessel functions depending on the input.
    get_group_velocity_eq(cp, fd, cp_prime, key):
        Calculates group velocity from current phase velocity and fd.
    solve_freq_equations(mode_type):
        Solves frequency equations and returns function of cp with respect to fd.
    """
    modes_nums = {
        'wavenumber': None,
        'circumfential_order': None,
    }

    def __init__(self, material, modes_nums, freq_thickness_max, cp_max, structure_mode=None, structure_freq=None, rows=None, columns=None, freq_thickness_points=None, cp_step=None):
        """
        Initialization method that calls function calculating results.

        This method is automatically called when the dataclass instance is initialized. 
        It calls init method of the base class Wave and sets additional parameters.
        """
        self.increasing_mode = 'T_(0,1)'
        self.modes_nums['wavenumber'], self.modes_nums['circumfential_order'] = modes_nums
        super().__init__(material, freq_thickness_max, cp_max, structure_mode, structure_freq, rows, columns, freq_thickness_points, cp_step)

    def __post_init__(self):
        """
        Post-initialization method that calls function calculating results.

        This method is automatically called after the dataclass instance is initialized. 
        It solves dispersion equation for torsional, longitudinal and flexural modes and 
        calculates wavestructure.
        """
        super().__post_init__()
        self.velocites_torsional = self.solve_freq_equations('torsional')
        self.velocites_longitudinal = None
        self.velocites_flexural = None

    def calculate_dispersion_components(self, phase_velocity: float, freq_thickness: float) -> float | float | float:
        """
        Calculates dispersion componetns k, alpha_sq and beta_sq.

        Parameters:
            phase_velocity (float): Current phase velocity value.
            freq_thickness (float): Current frequency and thickness product value.

        Returns:
            float: wavenumber
            float: alpha_sq component
            float: beta_sq component
        """
        omega = 2 * np.pi * (freq_thickness / self.material.thickness) 
        k = omega / phase_velocity
        alpha_sq = omega**2/self.material.longitudinal_wave_velocity**2 - k**2
        beta_sq = omega**2/self.material.shear_wave_velocity**2 - k**2
        return k, alpha_sq, beta_sq

    def get_group_velocity_eq(self, cp: float, fd: float, cp_prime: float) -> float:
        """
        Calculates group velocity from current phase velocity and fd.

        Parameters:
            phase_velocity (float): Current phase velocity value.
            freq_thickness (float): Current frequency and thickness product value.
            cp_prime (float): Current phase velocity derivative.

        Returns:
            float: Group velocity
        """
        cg = cp**2 * (cp - fd * cp_prime(fd))**-1
        return cg

    def calculate_bessel_functions(self, n: int, r: float, alpha_sq: float, beta_sq: float) -> dict:
        """
        Calculates bessel or modified bessel functions depending on the input.

        Parameters:
            n (int) : Current circumfential order
            r (float) : Radius of the cylinder
            alpha_sq (float) : Dispersion alpha_sq component
            beta_sq (float) : Dispersion beta_sq component

        Returns:
            dict: Dictionary of Bessel functions
        """
        alpha_r = np.sqrt(np.abs(alpha_sq)) * r
        beta_r = np.sqrt(np.abs(beta_sq)) * r

        if alpha_sq > 0 and beta_sq > 0:
            return {'alpha': (jv(n, alpha_r), yv(n, alpha_r), 1), 'beta': (jv(n, beta_r), yv(n, beta_r), 1)}
        elif alpha_sq < 0 and beta_sq > 0:
            return {'alpha': (iv(n, alpha_r), kv(n, alpha_r), -1), 'beta': (jv(n, beta_r), yv(n, beta_r), 1)}
        elif alpha_sq < 0 and beta_sq < 0:
            return {'alpha': (iv(n, alpha_r), kv(n, alpha_r), -1), 'beta': (iv(n, beta_r), kv(n, beta_r), -1)}
        else: 
            return {'alpha': (np.nan, np.nan, np.nan), 'beta': (np.nan, np.nan, np.nan)}

    def calculate_characteristic_equation(self, phase_velocity: float, freq_thickness: float, n: int) -> np.ndarray:
        """
        Calculates bessel or modified bessel functions depending on the input.

        Defines the characteristic equation matrix C and solves it.
        D matrix is extracted from C depending on the type of wave.
        Finally determinant is calculated and returned.

        Parameters:
            phase_velocity (float): Current phase velocity value.
            freq_thickness (float): Current frequency and thickness product value.
            n (int) : Current circumfential order

        Returns:
            np.ndarray: Matrix determinant.
        """
        k, alpha_sq, beta_sq = self.calculate_dispersion_components(phase_velocity, freq_thickness)

        alpha = np.sqrt(alpha_sq)
        beta = np.sqrt(beta_sq)

        # Create an 2x2 matrix filled with zeros
        D = np.zeros((2, 2), dtype=complex)

        # Coefficients of the matrix C
        for i, r in enumerate((self.material.inner_radius, self.material.outer_radius)):
            current_bessel = self.calculate_bessel_functions(n, r, alpha_sq, beta_sq)
            next_bessel = self.calculate_bessel_functions(n + 1, r, alpha_sq, beta_sq)

            Zn_a, Wn_a, lam_1 = current_bessel['alpha']
            Zn_b, Wn_b, lam_2 = current_bessel['beta']

            Zn_1_a, Wn_1_a, _ = next_bessel['alpha']
            Zn_1_b, Wn_1_b, _ = next_bessel['beta']

            D[i, :] = [
                    -(2 * n * (n - 1) - beta**2 * r**2) * Zn_b - 2 * lam_2 * beta * r * Zn_1_b, 
                    -(2 * n * (n - 1) - beta**2 * r**2) * Wn_b - 2 * beta * r * Wn_1_b 
            ]

        return np.linalg.det(D)

    def solve_freq_equations(self, mode_type: str) -> dict:
        """
        Solves frequency equations and returns a function of cp with respect to fd.

        The algorithm is based on the one presented in J.L. Rose's "Ultrasonic Guided Waves in Solid Media," 
        Chapter 6. The algorithm proceeds as follows:
        
        1. Choose a frequency-thickness product fd_0.
        2. Make an initial estimate of the phase velocity cp_0.
        3. Evaluate the signs of the left-hand sides of the frequency equations.
        4. Choose another phase velocity cp_1 > cp_0 and re-evaluate the signs of the frequency equations.
        5. Repeat steps (3) and (4) until a sign change occurs, assuming this happens between cp_n and cp_n+1.
        6. Use bisection to precisely locate the phase velocity in the interval cp_n < cp < cp_n+1 
           where the left-hand side of the equation is close enough to zero.
        7. After finding the root, continue searching at this fd for other roots according to steps (2) through (6).
        8. Choose another fd product and repeat steps (2) through (7).

        Parameters:
            mode_type (str): The type of modes shown on the plot: 'sym', 'anti', or 'both'.

        Returns:
            dict: The results of the frequency equation solutions, indexed by mode type.
        """
        fd_range = np.linspace(0, self.freq_thickness_max, self.freq_thickness_points)
        result = {}
        prefix = 'T_'
        initial_fd = []
        result[self.increasing_mode] = []

        for n in range(0, self.modes_nums['circumfential_order']):
            for fd in fd_range:
                mode = 2
                init_cp = 0
                init_cp2 = init_cp + self.cp_step

                while init_cp2 < self.cp_max:
                    cp_0 = self.calculate_characteristic_equation(init_cp, fd, n)
                    cp_1 = self.calculate_characteristic_equation(init_cp2, fd, n)

                    if mode < self.modes_nums['wavenumber']:
                        if not (np.isnan(cp_0) or np.isnan(cp_1)) and np.sign(cp_0) != np.sign(cp_1):
                            cp_root = bisect(f=self.calculate_characteristic_equation, a=init_cp, b=init_cp2, args=(fd, n))
                            if np.abs(self.calculate_characteristic_equation(cp_root, fd, n)) < self.tolerance:
                                key = f'{prefix}({n},{mode})'

                                if key not in result:
                                    result[key] = []                       

                                while result[key] and cp_root > result[key][-1][1]:
                                    mode += 1
                                    key = f'{prefix}({n},{mode})'
                                    if key not in result and mode:
                                        result[key] = []
                            
                                if not result[key]:
                                    if fd in initial_fd:
                                        break
                                    initial_fd.append(fd)

                                if mode < self.modes_nums['wavenumber']:
                                    result[key].append([fd, cp_root])
                                    mode += 1
                                else:
                                    result.pop(key, None)

                    init_cp += self.cp_step
                    init_cp2 += self.cp_step

        result[self.increasing_mode] = [[freq, self.material.shear_wave_velocity] for freq in fd_range]      
        result = self.interpolate_result(result)
        result = self.calculate_group_wavenumber(result)

        return result