from Plot import Plot
from Wave import Lambwave, Shearwave

### Lambwave
# Initialize wave
lamb = Lambwave('A_0', [500, 1000, 1500, 2000, 2500, 3000], 3, 2)

# Initialize Plot 
lamb_plotter = Plot(lamb, 'both')

# Add plots and show
lamb_plotter.add_plot('Phase')
lamb_plotter.add_plot('Group')
lamb_plotter.add_plot('Wavenumber')
lamb_plotter.add_plot('Wavestructure')
lamb_plotter.show_plots()

### Shearwave
# Initialize wave
shear = Shearwave('S_0', [2000, 3500, 5000, 7500, 9000, 10000], 3, 2)

# Initialize Plot 
shear_plotter = Plot(shear, 'both')

# Add plots and show
shear_plotter.add_plot('Phase')
shear_plotter.add_plot('Group')
shear_plotter.add_plot('Wavenumber')
shear_plotter.add_plot('Wavestructure')
shear_plotter.show_plots()
