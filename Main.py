from Material import Plate, Cylinder
from Plot import Plot
from Wave import Lambwave, Shearwave, AxialWave

# Initialize plate
plate = Plate(10, 6130, 3130, rayleigh_wave_velocity = 2881.6, name="Aluminium")

### Lambwave
# Initialize wave
lamb = Lambwave(plate, (5, 5), 15000, 12000, \
           'A_0', [500, 1000, 1500, 2000, 2500, 3000], 3, 2)

# Initialize Plot 
lamb_plotter = Plot(lamb, 'both', True, True)

# Add plots and show
lamb_plotter.add_plot('Phase')
lamb_plotter.add_plot('Group')
lamb_plotter.add_plot('Wavenumber')
lamb_plotter.add_plot('Wavestructure')
lamb_plotter.save_plots()
lamb_plotter.show_plots()
lamb_plotter.save_txt_results()

### Shearwave
# Initialize wave
shear = Shearwave(plate, (5, 5), 10000, 12000, \
           'A_0', [2000, 3500, 5000, 7500, 9000, 10000], 3, 2)

# Initialize Plot 
shear_plotter = Plot(shear, 'both', True, True)

# Add plots and show
shear_plotter.add_plot('Phase')
shear_plotter.add_plot('Group')
shear_plotter.add_plot('Wavenumber')
shear_plotter.add_plot('Wavestructure')
shear_plotter.show_plots()

# Initialize cyllinder
cylinder = Cylinder(1.22, 6290, 3230, 15.24, name="Steel")

# Initialize wave
axial = AxialWave(cylinder, (8, 1), 15000, 10000)

# Initialize Plot 
axial_plotter = Plot(axial, 'torsional', False, False)

# Add plots and show
axial_plotter.add_plot('Phase')
axial_plotter.add_plot('Group')
axial_plotter.add_plot('Wavenumber')
axial_plotter.show_plots()