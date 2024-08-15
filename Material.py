from dataclasses import dataclass, field

@dataclass()
class Material:
    thickness : float # plate thickness, in mm
    longitudinal_wave_velocity : float # Longitudinal wave velocity of the material, in m/s 
    shear_wave_velocity : float # Shear wave velocity of the material, in m/s 
    rayleigh_wave_velocity : float = field(default=None, kw_only=True) # Rayleigh wave velocity of the material, in m/s 
    name : str = field(default="no_material", kw_only=True) # name of the material 

@dataclass
class Plate(Material):
    half_thickness : float = field(init=False) # Half thickness, in mm

    def __post_init__(self):
        self.thickness = self.thickness/1e3
        self.half_thickness = self.thickness / 2

@dataclass
class Cylinder(Material):
    inner_radius : float
   
    def __post_init__(self):
       self.inner_radius, self.thickness = self.inner_radius/1e3, self.thickness/1e3 
       self.outer_radius = self.inner_radius + self.thickness