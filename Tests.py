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
                frequencies.append(frequency * 1000 * thickness)  # Convert MHz to kHz
                velocities.append(velocity * 1000)  # Convert m/ms to m/s
        if current_mode:
            data[current_mode] = (frequencies, velocities)  # Ensure last mode's data is stored
    return data

def setup_shear_wave(file_path):
    thickness = 10
    data_shear = load_data(file_path, thickness)
    plate_shear = Material(4460, 121e9, 0.3, thickness, 6060, 3230, None, "Titanium")
    shear = Shearwave(plate_shear, (3, 3), 10000, 10000, cp_step=50)
    return data_shear, shear

def setup_lamb_wave(file_path):
    thickness = 1
    data_lamb = load_data(file_path, thickness)
    plate_lamb = Material(1700, 42e9, 0.28, thickness, 5770, 3050, None, "Magnesium")
    lamb = Lambwave(plate_lamb, (5, 5), 10000, 10000, cp_step=50)
    return data_lamb, lamb

def plot_data(data, wave_model, title, save_path):
    plt.figure(figsize=(10, 8))
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
            sim_freqs, sim_vels = zip(*[(fd, cp) for fd, cp, *_ in wave_model.velocites_symmetric[mode]])
            if not legend_added['wave_dispersion_software']:
                plt.plot(sim_freqs, sim_vels, label='WaveDispersion Software', color=color_my_software)
                legend_added['wave_dispersion_software'] = True
            else:
                plt.plot(sim_freqs, sim_vels, color=color_my_software)

        if mode in wave_model.velocites_antisymmetric:
            sim_freqs, sim_vels = zip(*[(fd, cp) for fd, cp, *_ in wave_model.velocites_antisymmetric[mode]])
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

if __name__ == "__main__":
    data_shear, shear_wave = setup_shear_wave('validation/Titanium_Shear.txt')
    plot_data(data_shear, shear_wave, 'Shear Wave Test', 'validation/Titanium_Shear.png')
    plot_show()

    data_lamb, lamb_wave = setup_lamb_wave('validation/Magnesium_Lamb.txt')
    plot_data(data_lamb, lamb_wave, 'Lamb Wave Test', 'validation/Magnesium_Lamb.png')
    plot_show()
