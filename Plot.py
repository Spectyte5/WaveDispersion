from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import numpy as np
from Wave import Wave, Lambwave, Shearwave
from io import BytesIO
import base64, os

@dataclass
class Plot:
    wave : Wave
    mode_type : str
    cutoff_frequencies : bool = field(default=True)
    add_velocities : bool = field(default=True)
    path : str = field(default= "results")
    symmetric_style : dict = field(default_factory=lambda: {'color': 'green', 'linestyle': '-'})
    antisymmetric_style : dict = field(default_factory=lambda: {'color': 'purple', 'linestyle': '--'})
    dashed_line_style : dict = field(default_factory=lambda: {'color': 'black', 'linestyle': '--', 'linewidth': 0.5})
    continuous_line_style : dict = field(default_factory=lambda: {'color': 'black', 'linestyle': '-', 'linewidth': 0.75})
    in_plane_style : dict = field(default_factory=lambda: {'color': 'green', 'linestyle': '-', 'label': 'In plane'})
    out_of_plane_style : dict = field(default_factory=lambda: {'color': 'purple', 'linestyle': '--', 'label': 'Out of plane'})
    velocity_style : dict = field(default_factory=lambda: {'color': 'black', 'va': 'center'})
    padding_factor : dict = field(default_factory=lambda: {'x' : 1.00, 'y' : 1.05})
    get_figures = lambda _: [plt.figure(n) for n in plt.get_fignums()]

    def generate_latex(self, string):
        return r'$\mathregular{' + string + '}$'

    def find_max_value(self, index):
        max_value_sym = max(max([point[index] for point in values]) for values in self.wave.velocites_symmetric.values())
        max_value_anti = max(max([point[index] for point in values]) for values in self.wave.velocites_antisymmetric.values())
        return max(max_value_sym, max_value_anti)

    def draw_arrow(self, arrow):
        plt.axvline(x=arrow['x'], **self.dashed_line_style)
        plt.text(x=arrow['x'],y=arrow['y'], s=arrow['s'], va=arrow['dir'], ha='center', clip_on=True)

    def add_cutoff_frequencies(self, mode, max_value, plot_type : str):

        arrow_y, arrow_dir, arrow_s = (max_value, 'top', r'$\downarrow$') \
            if plot_type == 'Phase' else (0, 'bottom', r'$\uparrow$')

        if isinstance(self.wave, Shearwave):
            n = self.wave.get_converted_mode(mode)
            arrow_x = n*self.wave.velocities_dict['C_S']/2          
            self.draw_arrow({'x' : arrow_x, 'y' : arrow_y, 'dir': arrow_dir, 's' : arrow_s})

        else: 
            n = int(mode[2:]) + 1
            arrow_x = n*self.wave.velocities_dict['C_S'] if mode.startswith('S') else n*self.wave.velocities_dict['C_L']
            self.draw_arrow({'x' : arrow_x, 'y' : arrow_y, 'dir': arrow_dir, 's' : arrow_s})
        
            if n % 2 != 0:
                arrow_x = n*self.wave.velocities_dict['C_L']/2 if mode.startswith('S') else n*self.wave.velocities_dict['C_S']/2
                self.draw_arrow({'x' : arrow_x, 'y' : arrow_y, 'dir': arrow_dir, 's' : arrow_s})  


    def add_plate_velocities(self, plot_type):
        if plot_type == 'Phase' and self.add_velocities:
            for name, value in self.wave.velocities_dict.items():
                plt.axhline(value, **self.dashed_line_style)
                if name.endswith('R'):
                    ha = 'right'
                    cord = 0.5
                else: 
                    ha = 'left'
                    cord = self.wave.freq_thickness_max
                plt.text(cord, value, self.generate_latex(name), ha=ha, **self.velocity_style)  

    def plot_velocity(self, plot_type : str):
        title = f'{plot_type} Velocity' if plot_type != 'Wavenumber' else f'{plot_type[:4] + " " + plot_type[4:]}'
        plt.figure(num=title, figsize=(10, 6))
        plt.title(title)

        symmetric_lines, antisymmetric_lines = [], []
        index_map = {'Phase': 1,'Group': 2,'Wavenumber': 3 }
        max_value = self.find_max_value(index_map.get(plot_type))

        mode_mapping = {
            'symmetric': [(self.wave.velocites_symmetric, symmetric_lines, self.symmetric_style)],
            'antisymmetric': [(self.wave.velocites_antisymmetric, antisymmetric_lines, self.antisymmetric_style)],
            'both': [(self.wave.velocites_symmetric, symmetric_lines, self.symmetric_style),
                     (self.wave.velocites_antisymmetric, antisymmetric_lines, self.antisymmetric_style)]
            }

        selected_modes = mode_mapping[self.mode_type] 

        for data, lines, style in selected_modes:
            for mode, values in data.items():
                    x = [point[0] for point in values]
                    y = [point[index_map.get(plot_type)] for point in values]
                    if plot_type in ('Phase', 'Group') and self.cutoff_frequencies:
                        self.add_cutoff_frequencies(mode, max_value, plot_type)
                    line, = plt.plot(x, y, **style)
                    if isinstance(self.wave, Shearwave):
                        mode = 'SH_' + str(self.wave.get_converted_mode(mode))
                    plt.text(x[0], y[0], self.generate_latex(mode), ha='right', va='bottom', color=style['color']) 
                    lines.append(line) 
        

        self.add_plate_velocities(plot_type)
        
        # Create custom legend entries
        values, labels = ([mode_mapping[self.mode_type][0][1][0], mode_mapping[self.mode_type][1][1][0]], ['Symmetric', 'Antisymmetric']) \
            if self.mode_type == 'both' else ([mode_mapping[self.mode_type][0][1][0]], [self.mode_type])
        plt.legend(values, labels, loc='lower right')

        plt.xlim(0, self.wave.freq_thickness_max * self.padding_factor['x'])
        plt.ylim(0, max_value * self.padding_factor['y'])
        plt.xlabel('$\mathregular{f_d}$ (KHz x mm)')
        plt.ylabel('$\mathregular{c_p}$ (m/sec)') if plot_type != 'Wavenumber' else plt.ylabel('Wavenumber (1/m)') 
        
    def plot_wave_structure(self, title):
        # Create a figure and an array of subplots
        if self.wave.structure_result:
            fig, axes = plt.subplots(self.wave.rows, self.wave.columns, figsize=(10, 6),num=f"{title[:4] + ' ' + title[4:]}")

            # Flatten the axes array to easily access individual subplots
            axes = axes.flatten()

            # Loop through each key in the dictionary
            for i, key in enumerate(self.wave.structure_result.keys()):
                # Get the data associated with the current key
                data_list = self.wave.structure_result[key]

                # Determine which subplot to plot on
                ax = axes[i]

                # Loop through each element in the data list
                for entry in data_list:
                    u, w, x = entry  # Extracting the individual components from each element

                    # Plot the data on the current subplot
                    if np.all(np.iscomplex(u)):
                        ax.plot(np.imag(u), x, **self.in_plane_style)
                    else:
                        ax.plot(np.real(u), x, **self.in_plane_style)

                    if isinstance(self.wave, Lambwave):
                        if np.all(np.isreal(w)):
                            ax.plot(np.real(w), x, **self.out_of_plane_style)
                        else:
                            ax.plot(np.imag(w), x, **self.out_of_plane_style)

                # Disable the frame
                ax.set(frame_on=False) 

                # Moving x and y axis to 0,0
                ax.axhline(0, **self.continuous_line_style)  # Horizontal line at y=0
                ax.axvline(0, **self.continuous_line_style, ymin=0 + 0.05, ymax=1 - 0.05)  # Vertical line at x=0       
                ax.spines['left'].set_position(('data', 0))
                ax.spines['bottom'].set_position(('data', 0))

                # Plot ticks
                ax.set_yticks([-self.wave.material.half_thickness, 0, self.wave.material.half_thickness])
                ax.set_yticklabels(['-d/2', '0', 'd/2'])
                ax.text(ax.get_xlim()[1] / 2, self.wave.material.half_thickness / 5, 'u, w' 
                        if isinstance(self.wave, Lambwave) else 'u', ha='center', va='center')

                # Set title
                ax.set_title('$\mathregular{f_d}$' + f'={key} (kHz x mm)')

            # Adjust layout to prevent overlap of subplots add title
            plt.tight_layout
            fig.suptitle('Wave structure for mode ' + self.generate_latex(self.wave.structure_mode))

            # Get handles and labels for the first two legend entries
            handles, labels = ax.get_legend_handles_labels()
            fig.legend(handles, labels, loc='lower center', ncol=2)

    def add_plot(self, plot_type):
        if plot_type == 'Wavestructure':
            self.plot_wave_structure(plot_type)
        elif plot_type in ['Phase', 'Group', 'Wavenumber']:
            self.plot_velocity(plot_type)

    def save_plots(self, format='png', transparent=False, **kwargs):
        """
        Save plots in specified format.

        Parameters:
            format (str): Format to save the plots (e.g., 'png', 'pdf', 'svg', etc.).
            transparent (bool): Whether to save the plots with a transparent background.
            **kwargs: Additional keyword arguments to pass to the `savefig` method.

        Returns:
            plots (list[str])
        """
        if not os.path.exists(self.path):
            os.makedirs(self.path)

        plots = []

        for i, fig in enumerate(self.get_figures()):
            img = f"figure_{i+1}.{format}"
            img_path = os.path.join(self.path, img)
            fig.savefig(img_path, format=format, transparent=transparent, **kwargs)
            plots.append(img)

        return plots

    def save_txt_results(self):
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        
        wave_type = "Lambwaves" if isinstance(self.wave, Lambwave) else "Shearwaves"

        filename = f"{wave_type}_in_{self.wave.material.thickness}_mm_{self.wave.material.name}_plate.txt" 

        filepath = os.path.join(self.path, filename)

        mode_mapping = {
        'symmetric': [self.wave.velocites_symmetric],
        'antisymmetric': [self.wave.velocites_antisymmetric],
        'both': [self.wave.velocites_symmetric,self.wave.velocites_antisymmetric]
        }

        selected_modes = mode_mapping[self.mode_type] 

        with open(filepath, 'w') as file:
            for data in selected_modes:
                for mode, values in data.items():
                    header = "\t\t".join(["fd", "cp", "cg", "k"])
                    file.write(f"{mode}:\n")
                    file.write(f"\t{header}\t\n")
                    for set_values in values:
                        file.write("\t")
                        file.write("\t\t".join(f"{value:.2f}" for value in set_values))
                        file.write("\n")
                    file.write("\n")

        return filepath

    def show_plots(self):
        plt.show()