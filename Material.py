import numpy as np
from dataclasses import dataclass, field

@dataclass
class Material:
    density : float            
    youngs_modulus : float  # Young's modulus in Pascal
    poissons_ratio : float  # Poisson's ratio (unitless)
    thickness : float # plate thickness, in m
    longitudinal_wave_velocity : float = None # Longitudinal wave velocity of the material, in m/s 
    shear_wave_velocity : float = None # Shear wave velocity of the material, in m/s 
    rayleigh_wave_velocity : float = None # Rayleigh wave velocity of the material, in m/s 
    half_thickness : float = field(init=False) # Half thickness, in m
    name : str = field(default="No material name") # name of the material 

    def __post_init__(self):
        # temporary 
        self.half_thickness = self.thickness / 2
        self.shear_wave_velocity = np.sqrt(self.youngs_modulus / (2*self.density*(1+self.poissons_ratio)))
        self.longitudinal_wave_velocity=np.sqrt(self.youngs_modulus*(1-self.poissons_ratio) / (self.density*(1+self.poissons_ratio)*(1-2*self.poissons_ratio)))
        if self.rayleigh_wave_velocity is None:
            self.rayleigh_wave_velocity = self.shear_wave_velocity * ((0.862+1.14*self.poissons_ratio) / (1+self.poissons_ratio))