import unittest
import matplotlib.pyplot as plt
from Material import Material
from Wave import Lambwave, Shearwave

class TestWaveVisualization(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Set up common properties
        cls.setupShearWave()
        cls.setupLambWave()

    @staticmethod
    def setupShearWave():
        # Initialize ShearWave data and model
        thickness = 10
        TestWaveVisualization.data_shear = TestWaveVisualization.load_data('validation/Titanium_Shear.txt', thickness)
        TestWaveVisualization.plate_shear = Material(4460, 121e9, 0.3, thickness/1000, 6060, 3230, None, "Titanium")
        TestWaveVisualization.shear = Shearwave(TestWaveVisualization.plate_shear, (3, 3), 10000, 10000, 'S_0', [2000, 3500, 5000, 7500, 9000, 10000], 3, 2)

    @staticmethod
    def setupLambWave():
        thickness = 1
        # Initialize LambWave data and model
        TestWaveVisualization.data_lamb = TestWaveVisualization.load_data('validation/Magnesium_Lamb.txt', thickness)
        TestWaveVisualization.plate_lamb = Material(1700, 42e9, 0.28, thickness/1000, 5770, 3050, None, "Magnesium")
        TestWaveVisualization.lamb = Lambwave(TestWaveVisualization.plate_lamb, (5, 5), 10000, 11000, 'S_0', [2000, 3500, 5000, 7500, 9000, 10000], 3, 2)

    @staticmethod
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
            data[current_mode] = (frequencies, velocities)  # Store mode's data
        return data

    def test_shear_wave_visual_overlap(self):
        save_path = "validation/Titanium_Shear.png"
        self.plot_data(self.data_shear, self.shear, 'Shear Wave Test', save_path)

    def test_lamb_wave_visual_overlap(self):
        save_path = "validation/Magnesium_Lamb.png"
        self.plot_data(self.data_lamb, self.lamb, 'Lamb Wave Test', save_path)

    def plot_data(self, data, wave_model, title, save_path):
        plt.figure(figsize=(10, 8))
        for mode, (freqs, vels) in data.items():
            plt.plot(freqs, vels, label=f'Experimental {mode}')
            if mode in wave_model.velocites_symmetric:
                sim_freqs, sim_vels = zip(*[(fd, cp) for fd, cp, *_ in wave_model.velocites_symmetric[mode]])
                plt.plot(sim_freqs, sim_vels, label=f'Simulated {mode} (Symmetric)')
            if mode in wave_model.velocites_antisymmetric:
                sim_freqs, sim_vels = zip(*[(fd, cp) for fd, cp, *_ in wave_model.velocites_antisymmetric[mode]])
                plt.plot(sim_freqs, sim_vels, label=f'Simulated {mode} (Antisymmetric)')
        plt.xlabel('Frequency x Thickness (KHz x mm)')
        plt.ylabel('Velocity (m/s)')
        plt.title(title)
        plt.grid(True)
        plt.savefig(save_path)

if __name__ == '__main__':
    unittest.main()
