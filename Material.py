import numpy as np

class Material:
    density = 2700            
    youngs_modulus = 68.9e9  # Young's modulus in Pascal
    poissons_ratio = 0.33  # Poisson's ratio (unitless)
    longitudinal_wave_velocity = np.sqrt(youngs_modulus*(1-poissons_ratio) / (density*(1+poissons_ratio)*(1-2*poissons_ratio))) # Longitudinal wave velocity of the material, in m/s r
    shear_wave_velocity = np.sqrt(youngs_modulus / (2*density*(1+poissons_ratio))) # Shear wave velocity of the material, in m/s r
    rayleigh_wave_velocity = shear_wave_velocity * ((0.862+1.14*poissons_ratio) / (1+poissons_ratio)) # Rayleigh wave velocity of the material, in m/s r  
    thickness = 1e-2 # plate thickness, in m
    half_thickness = thickness / 2 # Half thickness, in m