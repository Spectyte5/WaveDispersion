from dataclasses import dataclass, field

@dataclass()
class Material:
    """
    Represents a material used in dispersion simulations.

    This class defines the physical properties of a material, including its wave velocities 
    and optional attributes. It is used to characterize materials in simulations where 
    wave propagation through different media is analyzed.

    Attributes:
        thickness (float): The thickness of the material plate, measured in millimeters (mm).
        longitudinal_wave_velocity (float): The velocity of longitudinal waves in the material, measured in meters per second (m/s).
        shear_wave_velocity (float): The velocity of shear waves in the material, measured in meters per second (m/s).
        rayleigh_wave_velocity (float, optional): The velocity of Rayleigh waves in the material, measured in meters per second (m/s). Defaults to None if not provided.
        name (str, optional): The name of the material. Defaults to "no_material" if not provided.

    Methods:
        None
    """
    thickness : float
    longitudinal_wave_velocity : float
    shear_wave_velocity : float
    rayleigh_wave_velocity : float = field(default=None, kw_only=True)
    name : str = field(default="no_material", kw_only=True)

@dataclass
class Plate(Material):
    """
    Represents a plate used in dispersion simulations.

    This class defines the physical properties of a plate.
    It is inheriting all material properties from the Material class.

    Attributes:
        half_thickness (float) : Plate thickness divided by 2, in m

    Methods:
        __post_init__ : changes units of thickness mm -> mm, initializes half_thickness
    """
    half_thickness : float = field(init=False) # Half thickness, in mm

    def __post_init__(self):
        """
        Post-initialization method that converts the thickness from millimeters to meters and calculates the
        half-thickness of the plate.

        This method is automatically called after the dataclass instance is initialized. It converts the 
        thickness attribute from millimeters to meters and computes the half-thickness of the plate.
        """
        self.thickness = self.thickness / 1e3
        self.half_thickness = self.thickness / 2

@dataclass
class Cylinder(Material):
    """
    Represents a cylinder used in dispersion simulations.

    This class defines the physical properties of a cylinder.
    It is inheriting all material properties from the Material class.

    Attributes:
        inner_radius (float) : Inner radius of the cylinder, in mm.

    Methods:
        __post_init__ : changes units of thickness mm -> m, initializes inner and outer radiuses.
    """
    inner_radius : float
   
    def __post_init__(self):
        """
        Post-initialization method that converts the thickness from millimeters to meters and calculates the
        the outer radius of the cylinder.

        This method is automatically called after the dataclass instance is initialized. It converts the 
        thickness and inner_ radius attributes from millimeters to meters and computes the outer radius.
        """
        self.inner_radius, self.thickness = self.inner_radius / 1e3, self.thickness / 1e3 
        self.outer_radius = self.inner_radius + self.thickness