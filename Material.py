from dataclasses import dataclass, field

@dataclass
class Material:
    density : float            
    youngs_modulus : float  # Young's modulus in Pascal
    poissons_ratio : float  # Poisson's ratio (unitless)
    thickness : float # plate thickness, in m
    longitudinal_wave_velocity : float # Longitudinal wave velocity of the material, in m/s 
    shear_wave_velocity : float # Shear wave velocity of the material, in m/s 
    rayleigh_wave_velocity : float # Rayleigh wave velocity of the material, in m/s 
    half_thickness : float = field(init=False) # Half thickness, in m
    name : str = field(default="no_material") # name of the material 

    def __post_init__(self):
        self.half_thickness = self.thickness / 2