from dataclasses import dataclass
import matplotlib.pyplot as plt

@dataclass
class Plot:
    max_fd : float
    symmetric_data : dict = None
    antisymmetric_data : dict = None
    horizontal_args : dict = None
    symmetric_color : str = 'green'
    antisymmetric_color : str = 'purple'

    def generate_latex(self, string):
        return r'$\mathregular{' + string + '}$'

    def show_plots(self):
        self.plot_velocity('Phase')
        self.plot_velocity('Group')
        self.plot_velocity('Wavenumber')
        plt.show()

    def find_max_value(self, index):
        max_value_sym = max(max([point[index] for point in values]) for values in self.symmetric_data.values())
        max_value_anti = max(max([point[index] for point in values]) for values in self.antisymmetric_data.values())
        return max(max_value_sym, max_value_anti)

    # rework
    def add_cutoff_frequencies(self, mode, max_value, plot_type : str):
        n = int(mode[2:]) + 1
        if plot_type == 'Phase':
            if mode[0] == 'S':
                arrow_x = n*self.horizontal_args['C_L']/2 
            else: 
                arrow_x = n*self.horizontal_args['C_S']/2
            arrow_y = max_value
            plt.axvline(x=arrow_x, color='black', linestyle='--', linewidth=0.5)
            plt.text(x = arrow_x,y=arrow_y, s=r'$\downarrow$', va='top', ha='center', clip_on=True)
            if n % 2 != 0:
                if mode[0] == 'S':
                    arrow_x = n*self.horizontal_args['C_S'] 
                else: 
                    arrow_x = n*self.horizontal_args['C_L']
                plt.axvline(x=arrow_x, color='black', linestyle='--', linewidth=0.5)
                plt.text(x=arrow_x,y=arrow_y, s=r'$\downarrow$', va='top', ha='center', clip_on=True)   
        elif plot_type == 'Group':
            if mode[0] == 'S':
                arrow_x = n*self.horizontal_args['C_L'] 
            else: 
                arrow_x = n*self.horizontal_args['C_S']
            arrow_y = 0
            plt.axvline(x=arrow_x, color='black', linestyle='--', linewidth=0.5)
            plt.text(x = arrow_x,y=arrow_y, s=r'$\uparrow$', va='bottom', ha='center', clip_on=True)
            if n % 2 != 0:
                if mode[0] == 'S':
                    arrow_x = n*self.horizontal_args['C_S']/2 
                else: 
                    arrow_x = n*self.horizontal_args['C_L']/2
                plt.axvline(x=arrow_x, color='black', linestyle='--', linewidth=0.5)
                plt.text(x = arrow_x,y=arrow_y, s=r'$\uparrow$', va='bottom', ha='center', clip_on=True)  

    def add_plate_velocities(self, plot_type):
        if plot_type == 'Phase' and self.horizontal_args:
            for name, value in self.horizontal_args.items():
                plt.axhline(value, color='black', linestyle='--', linewidth=0.5)
                if name.endswith('R'):
                    ha = 'right'
                    cord = 0.5
                else: 
                    ha = 'left'
                    cord = self.max_fd
                plt.text(cord, value, self.generate_latex(name), color='black', ha=ha, va='center')  

    def plot_velocity(self, plot_type : str):
        if plot_type == 'Wavenumber':
             title = f"{plot_type[:4] + ' ' + plot_type[4:]}"
        else:
            title = f"{plot_type} Velocity"
        plt.figure(num=title, figsize=(10, 6))
        symmetric_lines = []
        antisymmetric_lines = []
        index_map = {'Phase': 1,'Group': 2,'Wavenumber': 3 }

        max_value = self.find_max_value(index_map.get(plot_type))

        for data, color, lines in [(self.symmetric_data, self.symmetric_color, symmetric_lines), 
                           (self.antisymmetric_data, self.antisymmetric_color, antisymmetric_lines)]:
            for mode, values in data.items():
                x = [point[0] for point in values]
                y = [point[index_map.get(plot_type)] for point in values]
                self.add_cutoff_frequencies(mode, max_value, plot_type)
                line_style = '-' if data == self.symmetric_data else '--'
                line, = plt.plot(x, y, color=color, linestyle=line_style)
                lines.append(line)
                plt.text(x[0], y[0], self.generate_latex(mode), ha='right', va='bottom', color=color)  

        self.add_plate_velocities(plot_type)
        # Create custom legend entries
        plt.legend([symmetric_lines[0], antisymmetric_lines[0]], ['Symmetric', 'Antisymmetric'], loc='lower right')
        plt.xlim(0, self.max_fd)
        plt.ylim(0, max_value)
        plt.xlabel('Frequency (KHz x mm)')
        if plot_type == 'Wavenumber':
            plt.ylabel('Wavenumber (1/m)')
        else:
            plt.ylabel('$\mathregular{c_p}$ (m/sec)')
        plt.title(title)
        