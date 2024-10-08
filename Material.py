from dataclasses import dataclass, field

@dataclass
class Material:
    thickness : float # plate thickness, in mm
    longitudinal_wave_velocity : float # Longitudinal wave velocity of the material, in m/s 
    shear_wave_velocity : float # Shear wave velocity of the material, in m/s 
    rayleigh_wave_velocity : float = field(default=None) # Rayleigh wave velocity of the material, in m/s 
    half_thickness : float = field(init=False) # Half thickness, in mm
    name : str = field(default="no_material") # name of the material 

    def __post_init__(self):
        self.thickness = self.thickness/1e3
        self.half_thickness = self.thickness / 2