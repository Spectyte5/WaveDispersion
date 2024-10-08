import matplotlib.pyplot as plt
from Material import Material
from Wave import Lambwave, Shearwave

def load_data(file_path, thickness):
    data = {}
    current_mode = None
    frequencies = []
    velocities = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if any(prefix in line for prefix in ['S_', 'A_']):
                if current_mode is not None:
                    data[current_mode] = (frequencies, velocities)
                current_mode = line.split(':')[0]
                frequencies = []
                velocities = []
            else:
                frequency, velocity = map(float, line.split(', '))
                frequencies.append(frequency * 1000 * thickness)  # Convert MHz to kHz * mm
                velocities.append(velocity * 1000)  # Convert m/ms to m/s
        if current_mode:
            data[current_mode] = (frequencies, velocities)  # Ensure last mode's data is stored
    return data

def setup_shear_wave(file_path):
    thickness = 10
    data_shear = load_data(file_path, thickness)
    plate_shear = Material(thickness, 6060, 3230, None, "Titanium")
    shear = Shearwave(plate_shear, (3, 3), 10000, 10000, cp_step=50)
    return data_shear, shear

def setup_lamb_wave(file_path):
    thickness = 1
    data_lamb = load_data(file_path, thickness)
    plate_lamb = Material(thickness, 5770, 3050, None, "Magnesium")
    lamb = Lambwave(plate_lamb, (5, 5), 10000, 10000, cp_step=50)
    return data_lamb, lamb

def plot_data(data, wave_model, type, title, save_path, background=False):
    plt.figure(figsize=(10, 8))
    if background:
        plt.switch_backend('agg') 
    color_my_software = 'purple'
    color_wave_dispersion_software = 'green'
    legend_added = {'wave_dispersion_software': False, 'disperse_software': False}

    for mode, (freqs, vels) in data.items():
        if not legend_added['disperse_software']:
            plt.plot(freqs, vels, label='Disperse Software', color=color_wave_dispersion_software)
            legend_added['disperse_software'] = True
        else:
            plt.plot(freqs, vels, color=color_wave_dispersion_software)

        if mode in wave_model.velocites_symmetric:
            sim_freqs, sim_vels = zip(*[(fd, vel) for fd, vel, *_ in wave_model.velocites_symmetric[mode]]) if type == 'Phase' \
                else zip(*[(fd, vel) for fd, _, vel, *_ in wave_model.velocites_symmetric[mode]])

            if not legend_added['wave_dispersion_software']:
                plt.plot(sim_freqs, sim_vels, label='WaveDispersion Software', color=color_my_software)
                legend_added['wave_dispersion_software'] = True
            else:
                plt.plot(sim_freqs, sim_vels, color=color_my_software)

        if mode in wave_model.velocites_antisymmetric:
            sim_freqs, sim_vels = zip(*[(fd, vel) for fd, vel, *_ in wave_model.velocites_antisymmetric[mode]]) if type == 'Phase' \
                else zip(*[(fd, vel) for fd, _, vel, *_ in wave_model.velocites_antisymmetric[mode]])
            if not legend_added['wave_dispersion_software']:
                plt.plot(sim_freqs, sim_vels, label='WaveDispersion Software', color=color_my_software)
                legend_added['wave_dispersion_software'] = True
            else:
                plt.plot(sim_freqs, sim_vels, color=color_my_software)

    plt.xlabel('Frequency x Thickness (KHz x mm)')
    plt.ylabel('Velocity (m/s)')
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.savefig(save_path)

def plot_show():
    plt.show()

def plot_close_all():
    plt.close('all')

if __name__ == "__main__":
    # Phase velocity
    data_shear, shear_wave = setup_shear_wave('validation/Titanium_Shear_Phase.txt')
    plot_data(data_shear, shear_wave, 'Phase', 'Shear Wave Phase Velocity Test', 'validation/Titanium_Shear_Phase.png')
    plot_show()

    data_lamb, lamb_wave = setup_lamb_wave('validation/Magnesium_Lamb_Phase.txt')
    plot_data(data_lamb, lamb_wave, 'Phase', 'Lamb Wave Phase Velocity Test', 'validation/Magnesium_Lamb_Phase.png')
    plot_show()

    # Group velocity
    data_shear, shear_wave = setup_shear_wave('validation/Titanium_Shear_Group.txt')
    plot_data(data_shear, shear_wave, 'Group', 'Shear Wave Group Velocity Test', 'validation/Titanium_Shear_Group.png')
    plot_show()

    data_lamb, lamb_wave = setup_lamb_wave('validation/Magnesium_Lamb_Group.txt')
    plot_data(data_lamb, lamb_wave, 'Group', 'Lamb Wave Group Velocity Test', 'validation/Magnesium_Lamb_Group.png')
    plot_show()