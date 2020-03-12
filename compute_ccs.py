from scipy.stats import linregress
import numpy as np

class SteppedFieldCCS:
    """compute the ccs for the multi-fields (stepped field method)
    """
    def __init__(self, params):
        """
        params: mass, temperatures, pressures, voltages, arrival_time, drift_tube_length=90.33, neutral_mass=28.013
        """
        self.mass = params['mass']
        self.temperatures = params['temp']
        self.pressures = params['pressures']
        self.voltages = params['voltages']
        self.arrival_time = params['arrival_time']
        self.drift_tube_length = params['drift_tube_length']
        self.neutral_mass = params['neutral_mass']

    def compute(self):
        """compute the ccs values based on the multi-field parameters
        """
        # ========================
        # given parameters
        # ========================
        # mass: scalar
        # drift_tube_length (cm): scalar
        # temperatures, T(C): array --> T(K) = T(C)+273.15
        T_K = np.array(self.temperatures) + 273.15
        # pressures, P(torr): array --> P(Pa) = P(torr)/760*101325
        P_torr = np.array(self.pressures)
        P_Pa = P_torr / 760 * 101325
        # voltage_cell, Vcell: array --> E = Vcell / drift_tube_length
        Vcell = np.array(self.voltages)
        E = Vcell / self.drift_tube_length
        inv_E = 1.0 / (E * 100.0)
        # arrival_time (ms): array
        arrival_sec = np.array(self.arrival_time) / 1000
        # neutral_mass = 28.013 (N2 by default)
        # ========================
        # constant parameters
        # ========================
        # 1.60217657E-19 or 1.6021766208E-19
        e = 1.6021766208E-19
        charge_state = 1
        boltzmann_constant = 1.38064852E-23
        N0 = 101325/boltzmann_constant/273.15 # N0_(m-3)
        # ========================
        # computed parameters by given
        # ========================
        # P/V = P(torr) / Vcell
        P_V = P_torr / Vcell
        # E/N (Td) = E / P(torr) / 0.3535
        E_N = (E / P_torr) / 0.3535
        mass_in_kg = self.mass * 1.66054E-27
        neutral_mass_in_kg = self.neutral_mass * 1.66054E-27
        reduced_mass_in_kg = (mass_in_kg * neutral_mass_in_kg / (mass_in_kg + neutral_mass_in_kg))
        # ========================

        slope, intercept, r_value, p_value, std_err = linregress(P_V, arrival_sec)
        # drift_time (sec) = arrival_sec - intercept
        drift_time = arrival_sec - intercept

        # compute CCS by Mason-Schamp Equation
        # ccs = 3 * e / 16 / N0 * np.sqrt(2 * np.pi / reduced_mass_in_kg / boltzmann_constant / T_K) \
        # * drift_time * 760 * T_K * Vcell / (drift_tube_length / 100)**2 / P_torr / 273.15 * 1E20
        K0 = self.drift_tube_length * self.drift_tube_length / slope * 273.15 / 760 / np.mean(T_K)
        ccs = 3 * e / 16 / N0 / K0 / 0.0001 * np.sqrt(2 * np.pi / (boltzmann_constant * reduced_mass_in_kg * np.mean(T_K))) * 1e20
        properties = {'mass':self.mass, 'slope': slope, 'intercept': intercept, 'r_value': r_value, 'p_value':p_value, 'k0':K0, 'ccs':ccs}
        return ccs, properties