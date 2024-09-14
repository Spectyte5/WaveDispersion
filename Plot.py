from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import numpy as np
from Wave import Wave, Lambwave, Shearwave
import os
from datetime import datetime

@dataclass
class Plot:
    """
    A class to manage and manipulate Matplotlib plots.

    This class is used for managing result display in both text and plot forms. 
    Supports creating, displaying, and managing plots in a variety of scenarios.
    Provides methods for generating LaTeX formatted strings for mathematical
    expressions, switching the Matplotlib backend to a non-interactive mode and saving, 
    showing and closing all open plot windows.

    Attributes:
        wave (Wave) : Wave class object providing the results to be plotted.
        mode_type (str) : Types of the modes shown on the plot: sym, anti or both.
        cutoff_frequencies (bool, optional) : Show cutoff frequencies on the plot.
        add_velocities (bool, optional) : Show plate velocities on the plot.
        path (str, optional) : Path to save the result files.
        symmetric_style (dict, optional) : Dict with style kwargs for symmetric modes.
        antisymmetric_style (dict, optional) : Dict with style kwargs for antisymmetric modes.
        torsional_style (dict, optional) : Dict with style kwargs for torsional modes.
        longitudinal_style (dict, optional) : Dict with style kwargs for longitudinal modes.
        flexural_style (dict, optional) : Dict with style kwargs for flexural modes.
        dashed_line_style (dict, optional) : Dict with style kwargs for all dashed lines on plots.
        continuous_line_style (dict, optional) : Dict with style kwargs for all continous lines on plots.
        in_plane_style (dict, optional) : Dict with style kwargs for the in_plane component on wavestructure plots.
        out_of_plane_style (dict, optional) : Dict with style kwargs for the out_of_plane component on wavestructure plots.
        velocity_style (dict, optional) : Dict with style kwargs for plate velocity option.
        padding_factor (float, optional) : Padding thickness for plots.
        get_figures (lambda) : Gets all plots in an array

    Methods:
        _generate_latex(string):
            Converts a string into LaTeX format for displaying mathematical expressions.
        switch_backend():
            Switches Matplotlib's backend to 'agg' for non-interactive use.
        close_all_plots():
            Closes all open Matplotlib plot windows.
        _find_max_value(index): 
            Find max value for the given wave parameter.
        _draw_arrow(arrow): 
            Draws the arrow on the plot.
        _add_cutoff_frequencies(mode, max_value, plot_type): 
            Add cutoff frequencies to the plot.
        _add_plate_velocities(plot_type): 
            Add plate velocities to the plot. 
        _plot_velocity(plot_type):
            Plot the given velocity or wavenumber.
        _plot_wave_structure(title): 
            Plot the given velocity or wavenumber.
        add_plot(plot_type): 
            Add plot velocity.
        save_plots(format, transparent, **kwargs):
            Save plots in specified format.
        save_txt_results(date):
            Save results in text file.
        show_plots(): 
            Displays all currently active Matplotlib plots.
    """
    wave : Wave
    mode_type : str
    cutoff_frequencies : bool = field(default=True)
    add_velocities : bool = field(default=True)
    path : str = field(default= "results")
    symmetric_style : dict = field(default_factory=lambda: {'color': 'green', 'linestyle': '-'})
    antisymmetric_style : dict = field(default_factory=lambda: {'color': 'purple', 'linestyle': '--'})
    torsional_style : dict = field(default_factory=lambda: {'color': 'green', 'linestyle': '-'})
    longitudinal_style : dict = field(default_factory=lambda: {'color': 'purple', 'linestyle': '--'})
    flexural_style : dict = field(default_factory=lambda: {'color': 'orange', 'linestyle': '-'})
    dashed_line_style : dict = field(default_factory=lambda: {'color': 'black', 'linestyle': '--', 'linewidth': 0.5})
    continuous_line_style : dict = field(default_factory=lambda: {'color': 'black', 'linestyle': '-', 'linewidth': 0.75})
    in_plane_style : dict = field(default_factory=lambda: {'color': 'green', 'linestyle': '-', 'label': 'In plane'})
    out_of_plane_style : dict = field(default_factory=lambda: {'color': 'purple', 'linestyle': '--', 'label': 'Out of plane'})
    velocity_style : dict = field(default_factory=lambda: {'color': 'black', 'va': 'center'})
    padding_factor : dict = field(default_factory=lambda: {'x' : 1.00, 'y' : 1.05})
    get_figures = lambda _: [plt.figure(n) for n in plt.get_fignums()]

    def _generate_latex(self, string: str) -> str:
        """
        Converts a string into LaTeX format for displaying mathematical expressions.

        This method handles strings with and without subscripts. If the input string
        contains an underscore ('_'), it is interpreted as a subscript, and the string
        is formatted accordingly. If there is no underscore, the string is formatted 
        as a regular mathematical expression.

        Example:
            'x_2' -> r'$\mathregular{x_{2}}$'
            'y' -> r'$\mathregular{y}$'

        Parameters:
            string (str): The input string to be converted into LaTeX format.

        Returns:
            str: The LaTeX formatted string.
        """
        if "_" in string:
            base, sub = string.split("_")
            return r'$\mathregular{' + base + '_{' + sub + '}}$'
        else:
            return r'$\mathregular{' + string + '}$'

    def switch_backend(self):
        """
        Switches Matplotlib's backend to 'agg'.

        The 'agg' backend is a non-interactive backend that can be used for 
        generating plot images without displaying them. This is useful for 
        saving plots to files in a script or when working in environments where 
        graphical display is not available.

        Parameters:
            None

        Returns:
            None
        """
        plt.switch_backend('agg')

    def close_all_plots(self):
        """
        Closes all open Matplotlib plots.

        This method calls `plt.close('all')` to close all figure windows. It is useful 
        for cleaning up after generating plots to free up resources or to start a new 
        plotting session without interference from previous figures.

        Parameters:
            None

        Returns:
            None
        """
        plt.close('all')

    def _find_max_value(self, index: int) -> float:
        """
        Find max value for the given wave parameter.

        This method find the max value for the given wave parameter.
        The parameter can be phase, group velocity or wavenumber. 
        It returns single highest value found in both symmetric and antisymmetric modes.

        Parameters:
            None

        Returns:
            max_value (float)
        """
        max_value_sym = max(max([point[index] for point in values]) for values in self.wave.velocites_symmetric.values())
        max_value_anti = max(max([point[index] for point in values]) for values in self.wave.velocites_antisymmetric.values())
        return max(max_value_sym, max_value_anti)

    def _draw_arrow(self, arrow: dict):
        """
        Closes all open Matplotlib plots.

        This method draws arrow dictionary using combination of plt.axvline 
        for the arrow and text for the arrow sign. This function is used 
        mainly for marking cutoff frequencies using arrows.


        Parameters:
            arrow (dict): Arrow to be drown with given parameters.

        Returns:
            None
        """
        plt.axvline(x=arrow['x'], **self.dashed_line_style)
        plt.text(x=arrow['x'],y=arrow['y'], s=arrow['s'], va=arrow['dir'], ha='center', clip_on=True)

    def _add_cutoff_frequencies(self, mode: str, max_value: float, plot_type : str):
        """
        Add cutoff frequencies to the plot.

        This method is used for marking cutoff frequencies on the plot.
        Calculates cutoff frequencies from the wave class parameters.

        Parameters:
            mode (str): Which mode the cutoff frequency is calcuated for.
            max_value (float): Max value, used for positioning of the arrow.
            plot_type (str): Type of the plot, ex. Phase velocity plot.

        Returns:
            None
        """
        arrow_y, arrow_dir, arrow_s = (max_value, 'top', r'$\downarrow$') \
            if plot_type == 'Phase' else (0, 'bottom', r'$\uparrow$')

        if isinstance(self.wave, Shearwave):
            n = self.wave.get_converted_mode(mode)
            arrow_x = n*self.wave.velocities_dict['C_S']/2          
            self._draw_arrow({'x' : arrow_x, 'y' : arrow_y, 'dir': arrow_dir, 's' : arrow_s})

        else: 
            n = int(mode[2:]) + 1
            arrow_x = n*self.wave.velocities_dict['C_S'] if mode.startswith('S') else n*self.wave.velocities_dict['C_L']
            self._draw_arrow({'x' : arrow_x, 'y' : arrow_y, 'dir': arrow_dir, 's' : arrow_s})
        
            if n % 2 != 0:
                arrow_x = n*self.wave.velocities_dict['C_L']/2 if mode.startswith('S') else n*self.wave.velocities_dict['C_S']/2
                self._draw_arrow({'x' : arrow_x, 'y' : arrow_y, 'dir': arrow_dir, 's' : arrow_s})  

    def _add_plate_velocities(self, plot_type: str):
        """
        Add plate velocities to the plot.

        This method is used for marking bulk plate velocities on the plot.
        It marks the velocities when the user toggles the option on.

        Parameters:
            plot_type (str): Type of the plot, ex. Phase velocity plot.

        Returns:
            None
        """
        if plot_type == 'Phase' and self.add_velocities:
            for name, value in self.wave.velocities_dict.items():
                if not value:
                    continue
                plt.axhline(value, **self.dashed_line_style)
                if name.endswith('R'):
                    ha = 'right'
                    cord = 0.5
                else: 
                    ha = 'left'
                    cord = self.wave.freq_thickness_max
                plt.text(cord, value, self._generate_latex(name), ha=ha, **self.velocity_style)  

    def _plot_velocity(self, plot_type : str):
        """
        Plot the given velocity or wavenumber.

        This method is used for generation of the velocity plots.
        It gets parameters from the self.wave plots them and adds legend, 
        labels, plot limits and the format specified by user in kwargs.

        Parameters:
            plot_type (str): Type of the plot, ex. Phase velocity plot.

        Returns:
            None
        """
        title = f'{plot_type} Velocity' if plot_type != 'Wavenumber' else f'{plot_type[:4] + " " + plot_type[4:]}'
        plt.figure(num=title, figsize=(10, 6))
        plt.title(title)

        symmetric_lines, antisymmetric_lines = [], []
        torsional_lines, longitudinal_lines, flexural_lines = [], [], []
        index_map = {'Phase': 1,'Group': 2,'Wavenumber': 3 }
        max_value = self._find_max_value(index_map.get(plot_type))

        mode_mapping =  {
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
                        self._add_cutoff_frequencies(mode, max_value, plot_type)
                    line, = plt.plot(x, y, **style)
                    if isinstance(self.wave, Shearwave):
                        mode = 'SH_' + str(self.wave.get_converted_mode(mode))
                    plt.text(x[0], y[0], self._generate_latex(mode), ha='right', va='bottom', color=style['color']) 
                    lines.append(line) 
        

        self._add_plate_velocities(plot_type)
        
        # Create custom legend entries
        values, labels = ([mode_mapping[self.mode_type][0][1][0], mode_mapping[self.mode_type][1][1][0]], ['Symmetric', 'Antisymmetric']) \
            if self.mode_type == 'both' else ([mode_mapping[self.mode_type][0][1][0]], [self.mode_type])
        plt.legend(values, labels, loc='upper left')

        plt.xlim(0, self.wave.freq_thickness_max * self.padding_factor['x'])
        plt.ylim(0, max_value * self.padding_factor['y'])
        plt.xlabel('$\mathregular{f_d}$ (KHz x mm)')
        plt.ylabel('$\mathregular{c_p}$ (m/sec)') if plot_type != 'Wavenumber' else plt.ylabel('Wavenumber (1/m)') 
        
    def _plot_wave_structure(self, title : str):
        """
        Plot the given velocity or wavenumber.

        This method is used for generation of the wave structure plots.
        It gets parameters from the self.wave, plots them and adds legend, 
        labels, plot limits and the format specified by user in kwargs.

        Parameters:
            title (str): Type of the plot, ex. Phase velocity plot.

        Returns:
            None
        """
        if self.wave.structure_result:
            if self.wave.rows == 1 and self.wave.columns == 1:
                fig, axes = plt.subplots(1, 1, figsize=(10, 6), num=f"{title[:4] + ' ' + title[4:]}")
                axes = [axes]
            else:
                fig, axes = plt.subplots(self.wave.rows, self.wave.columns, figsize=(10, 6), num=f"{title[:4] + ' ' + title[4:]}")
                axes = axes.flatten()

            for i, key in enumerate(self.wave.structure_result.keys()):
                data_list = self.wave.structure_result[key]

                ax = axes[i]

                for entry in data_list:
                    u, w, x = entry  
                    if np.all(np.iscomplex(u)):
                        ax.plot(np.imag(u), x, **self.in_plane_style)
                    else:
                        ax.plot(np.real(u), x, **self.in_plane_style)

                    if isinstance(self.wave, Lambwave):
                        if np.all(np.isreal(w)):
                            ax.plot(np.real(w), x, **self.out_of_plane_style)
                        else:
                            ax.plot(np.imag(w), x, **self.out_of_plane_style)

                ax.set(frame_on=False) 
                ax.axhline(0, **self.continuous_line_style)
                ax.axvline(0, **self.continuous_line_style, ymin=0 + 0.05, ymax=1 - 0.05)      
                ax.spines['left'].set_position(('data', 0))
                ax.spines['bottom'].set_position(('data', 0))

                ax.set_yticks([-self.wave.material.half_thickness, 0, self.wave.material.half_thickness])
                ax.set_yticklabels(['-d/2', '0', 'd/2'])
                ax.text(ax.get_xlim()[1] / 2, self.wave.material.half_thickness / 5, 'u, w' 
                        if isinstance(self.wave, Lambwave) else 'u', ha='center', va='center')

                ax.set_title('$\mathregular{f_d}$' + f'={key} (kHz x mm)')

            plt.tight_layout
            mode = 'SH_' + str(self.wave.get_converted_mode(self.wave.structure_mode)) if isinstance(self.wave, Shearwave) \
               else self.wave.structure_mode
            fig.suptitle('Wave structure for mode ' + self._generate_latex(mode))

            handles, labels = ax.get_legend_handles_labels()
            fig.legend(handles, labels, loc='lower center', ncol=2)

    def add_plot(self, plot_type: str):
        """
        Add plot velocity.

        This method is used for adding desired plots to simulation.
        It checks the plot_type parameter and calls the appropriate method that 
        generates the plot.

        Parameters:
            plot_type (str): Type of the plot, ex. Phase velocity plot.

        Returns:
            None
        """
        if plot_type == 'Wavestructure':
            self._plot_wave_structure(plot_type)
        elif plot_type in ['Phase', 'Group', 'Wavenumber']:
            self._plot_velocity(plot_type)

    def save_plots(self, format: str='png', transparent: bool=False, **kwargs) -> list[str]:
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

    def save_txt_results(self, date=False) -> str:
        """
        Save results in text file.

        The name of the file is generated automatically.
        
        Example:
             filename = f"Shearwaves_in_10_mm_Aluminium_plate.txt" 

        Parameters:
            date (bool): Whether to save the plots with a date to prevent overriding.

        Returns:
            filepath (str)
        """
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        
        wave_type_map = {
            Lambwave: "Lambwaves",
            Shearwave: "Shearwaves"
        }

        wave_type = wave_type_map.get(type(self.wave), "Unknown wave type")

        filename = f"{wave_type}_in_{self.wave.material.thickness}_mm_{self.wave.material.name}_plate"
        
        if date:
            now = datetime.now()
            file_name += f"_{now}"
        
        filename += ".txt"

        filepath = os.path.join(self.path, filename)

        mode_mapping = { 
            'symmetric': [self.wave.velocites_symmetric],
            'antisymmetric': [self.wave.velocites_antisymmetric],
            'both': [self.wave.velocites_symmetric,self.wave.velocites_antisymmetric],
        }

        selected_modes = mode_mapping[self.mode_type] 

        with open(filepath, 'w') as file:
            for data in selected_modes:
                for mode, values in data.items():
                    if isinstance(self.wave, Shearwave):
                        mode = 'SH_' + str(self.wave.get_converted_mode(mode))
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
        """
        Displays all currently active Matplotlib plots.

        This method calls the `plt.show()` function to render and display
        any plots that have been created.

        Parameters:
            None

        Returns:
            None
        """
        plt.show()