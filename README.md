# WaveDispersion

## Goal 
Project idea was to create a alternative numerical libary for providing dispersion characteristics for different wave types and boundaries.
Dispersion curves show a relation between relevant wave speed and its frequency.
They are commonly used in different fields of engineering and science for mode identification, 
solution optimization and generally representing the behavior of waves in given medium.

## Description
This application, allows user to visualize dispersion curves for different wave types (Lamb waves, Shear-Horizontal waves, etc.) 
and in different mediums (such traction-free isotropic plates, hollow cyllinders from different materials, etc.).
User has many options and parameters that allow checking wave behavior under desired condition and also changing visual side of plots.

## Prerequisites
```Python
Python >= 3.11
```
## Setup
Install the package using pip - package installer for Python:
```Python
pip3 install ./WaveDispersion
```
## Usage
Import modules:
```Python
from Material import Plate, Cylinder
from Plot import Plot
from Wave import Lambwave, Shearwave, Axialwave
```
Initialize plate or cylinder object:

### Common parameters
- **thickness** (`float`):  
  The thickness of the material plate, measured in millimeters (mm).
- **longitudinal_wave_velocity** (`float`):  
  The velocity of longitudinal waves in the material, measured in meters per second (m/s).
- **shear_wave_velocity** (`float`):  
  The velocity of shear waves in the material, measured in meters per second (m/s).
- **rayleigh_wave_velocity** (`float`, optional):  
  The velocity of Rayleigh waves in the material, measured in meters per second (m/s). Defaults to `None` if not provided.
- **name** (`str`, optional):  
  The name of the material. Defaults to `"no_material"` if not provided.

### Plate
```Python
plate = Plate(10, 6130, 3130, rayleigh_wave_velocity = 2881.6, name="Aluminium")
```

### Cylinder
```Python
cylinder = Cylinder(1.22, 6290, 3230, 15.24, name="Steel")
```
#### Parameters
- **inner_radius** (`float`):
  Inner radius of the cylinder, in mm.

Create wave object:
```Python
# Lambwave
lamb = Lambwave(plate, (5, 5), 15000, 12000, \
           'A_0', [500, 1000, 1500, 2000, 2500, 3000], 3, 2)
# Shearwave
shear = Shearwave(plate, (5, 5), 10000, 12000, \
           'SH_1', [2000, 3500, 5000, 7500, 9000, 10000], 3, 2)
# Axialwave
axial = Axialwave(cylinder, (8, 1), 15000, 10000)
```
### Parameters
- **material** (`Material`):  
  Material class object providing the material information.

- **modes_nums** (`tuple`):  
  Number of symmetric and antisymmetric modes (*Plate*), circumferencial order and number of wavenumber modes (*Cylinder*).

- **freq_thickness_max** (`int`):  
  Maximum value of Frequency x Thickness [kHz x mm].

- **cp_max** (`int`):  
  Maximum value of phase velocity [m/s].

- **structure_mode** (`str`):  
  Specifies the mode to be used for the wavestructure plot.

- **structure_freq** (`array`):  
  Array of frequencies at which to check wavestructure.

- **rows** (`int`):  
  Number of rows for the wavestructure plot.

- **columns** (`int`):  
  Number of columns for the wavestructure plot.

- **freq_thickness_points** (`int`):  
  Number of frequency x thickness points to calculate.

- **cp_step** (`int`):  
  Step size between phase velocity points to be checked.

- **increasing_mode** (`str`):  
  Specifies a unique mode that increases rather than decreases with frequency x thickness.

Initialize Plot object
```Python
lamb_plotter = Plot(lamb, 'both', True, True)
```
### Parameters

- **wave** (`Wave`):  
  Wave class object providing the results to be plotted.

- **mode_type** (`str`):  
  Types of the modes shown on the plot: `sym`, `anti`, or `both`.

- **cutoff_frequencies** (`bool`, optional):  
  Show cutoff frequencies on the plot.

- **add_velocities** (`bool`, optional):  
  Show plate velocities on the plot.

- **path** (`str`, optional):  
  Path to save the result files.

- **symmetric_style** (`dict`, optional):  
  Dictionary with style keyword arguments for symmetric modes.

- **antisymmetric_style** (`dict`, optional):  
  Dictionary with style keyword arguments for antisymmetric modes.

- **torsional_style** (`dict`, optional):  
  Dictionary with style keyword arguments for torsional modes.

- **longitudinal_style** (`dict`, optional):  
  Dictionary with style keyword arguments for longitudinal modes.

- **flexural_style** (`dict`, optional):  
  Dictionary with style keyword arguments for flexural modes.

- **dashed_line_style** (`dict`, optional):  
  Dictionary with style keyword arguments for all dashed lines on plots.

- **continuous_line_style** (`dict`, optional):  
  Dictionary with style keyword arguments for all continuous lines on plots.

- **in_plane_style** (`dict`, optional):  
  Dictionary with style keyword arguments for the in-plane component on wavestructure plots.

- **out_of_plane_style** (`dict`, optional):  
  Dictionary with style keyword arguments for the out-of-plane component on wavestructure plots.

- **velocity_style** (`dict`, optional):  
  Dictionary with style keyword arguments for the plate velocity option.

- **padding_factor** (`float`, optional):  
  Padding thickness for plots.

- **axial_factor** (`float`, optional):  
  Axial wave factor for scaling the distance from the y-axis.

Call one of possible plot methods:
```Python
lamb_plotter.add_plot('Phase')
lamb_plotter.add_plot('Group')
lamb_plotter.add_plot('Wavenumber')
lamb_plotter.add_plot('Wavestructure')
lamb_plotter.save_plots()
lamb_plotter.show_plots()
lamb_plotter.save_txt_results()
```
### Methods

- **`switch_backend()`**:  
  Switches Matplotlib's backend to `'agg'`.

  The `'agg'` backend is a non-interactive backend that can be used for generating plot images without displaying them. This is useful for saving plots to files in a script or when working in environments where graphical display is not available.

  **Parameters**:  
  None

  **Returns**:  
  None

- **`close_all_plots()`**:  
  Closes all open Matplotlib plots.

  This method calls `plt.close('all')` to close all figure windows. It is useful for cleaning up after generating plots to free up resources or to start a new plotting session without interference from previous figures.

  **Parameters**:  
  None

  **Returns**:  
  None

- **`add_plot(plot_type)`**:  
  Adds a plot based on the specified type.

  This method is used for adding desired plots to the simulation. It checks the `plot_type` parameter and calls the appropriate method to generate the plot.

  **Parameters**:  
  - `plot_type` (`str`): Type of the plot, e.g., "Phase velocity plot".

  **Returns**:  
  None

- **`save_plots(format, transparent, **kwargs)`**:  
  Saves plots in the specified format.

  **Parameters**:  
  - `format` (`str`): Format to save the plots (e.g., `'png'`, `'pdf'`, `'svg'`, etc.).
  - `transparent` (`bool`): Whether to save the plots with a transparent background.
  - `**kwargs`: Additional keyword arguments to pass to the `savefig` method.

  **Returns**:  
  - `list[str]`: List of paths to the saved plot files.

- **`save_txt_results(date)`**:  
  Saves results to a text file.

  The name of the file is generated automatically.

  **Example**:  
  `filename = "Shearwaves_in_10_mm_Aluminium_plate.txt"`

  **Parameters**:  
  - `date` (`bool`): Whether to include the current date in the filename to prevent overriding.

  **Returns**:  
  - `str`: Path to the saved text file.

- **`show_plots()`**:  
  Displays all currently active Matplotlib plots.

  This method calls the `plt.show()` function to render and display any plots that have been created.

  **Parameters**:  
  None

  **Returns**:  
  None

## References
1. This repository was created as master thesis made by Mechatronic Engineering student at Faculty of Mechanical Engineering and Robotics. 
2. Project was done at AGH University of Science and Technology in Cracow under the supervision of Dr hab. inż. Wiesław Staszewski. 
3. The equations and algorithms were obtained from the great book: Rose, J. L., Ultrasonic Guided Waves in Solid Media, Cambridge University Press, 1999.
4. The algorithm was based on the one from [LambWaveDispersion](https://github.com/franciscorotea/Lamb-Wave-Dispersion/tree/master).
5. The results were validated using [Disperse Software](http://www.disperse.software).
