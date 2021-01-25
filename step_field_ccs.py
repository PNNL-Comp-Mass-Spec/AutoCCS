from scipy.stats import linregress
import numpy as np

def mass_ppm_error(x, mass):
    return abs(x - mass) / mass * 1e6

class SteppedFieldCCS:
    """compute the ccs for the multi-fields (stepped field method)
    """
    def __init__(self, meta_df, adduct_mass, old_drift_tube_length, charge_state=1):
        """
        metadata: a dictionary for
            {mass, temperatures, pressures, voltages, arrival_time}
        """
        self._metadata = {}
        # self.mz = params['mass']
        # self.temperatures = params['temp']
        # self.pressures = params['pressures']
        # self.voltages = params['voltages']
        # self.arrival_time = params['arrival_time']
        # self.drift_tube_length = params['drift_tube_length']
        # self.neutral_mass = params['neutral_mass']
        self._metadata['adduct_mz'] = adduct_mass
        self._metadata['num_features'] = len(list(meta_df.frame.drop_duplicates()))
        
        self._mppid = []
        self._dt = []
        self._num_isotopes = []
        self._intensity_org = []
        self._intensity_z = []
        self._intensity = []
        self._mass_ppm_error = []
        for feature in meta_df.itertuples():
            self._metadata['dt_' + str(feature.frame)] = feature.dt
            self._metadata['intensity_org_' + str(feature.frame)] = feature.intensity_org
            self._metadata['intensity_z_' + str(feature.frame)] = feature.intensity_z
            self._metadata['intensity_' + str(feature.frame)] = feature.intensity
            self._metadata['mass_error_' + str(feature.frame)] = mass_ppm_error(feature.mz, adduct_mass)
            
            self._mppid.append(feature.mppid)
            self._dt.append(feature.dt)
            self._intensity_org.append(feature.intensity_org)
            self._intensity.append(feature.intensity)
            self._intensity_z.append(feature.intensity_z)
            self._mass_ppm_error.append(self._metadata['mass_error_' + str(feature.frame)])

            if 'num_isotopes' in meta_df.columns: 
                self._metadata['num_isotopes_' + str(feature.frame)] = feature.num_isotopes
                self._num_isotopes.append(feature.num_isotopes)
        self.temperatures = meta_df.ImsTemperature.tolist()
        self._pressures = meta_df.ImsPressure.tolist()
        self._fields = meta_df.ImsField.tolist()
        self.voltages = (meta_df.ImsField*old_drift_tube_length).tolist()
        self._arrival_time = meta_df.dt.tolist()
        self.mz = adduct_mass
        self.charge_state = charge_state

        # params['temp'] = df.ImsTemperature.tolist()
        # params['pressures'] = df.ImsPressure.tolist()
        # params['voltages'] = (df.ImsField*config_params['old_drift_tube_length']).tolist()  ## 10.869 * (78.12 / 78.236) = 10.853 for correction
        # params['arrival_time'] = df.dt.tolist()
        # params['neutral_mass'] = config_params['neutral_mass']
        # params['drift_tube_length'] = config_params['drift_tube_length']
        # params['mz'] = ion_mz

    @property
    def r2(self):
        return self._metadata['r2']
    
    @property
    def ccs(self):
        return self._metadata['ccs']

    @property
    def p_v(self):
        return self._p_v

    @property
    def arrival_time(self):
        return self._arrival_time

    @property
    def fields(self):
        return self._fields
    @property
    def pressures(self):
        return self._pressures
    @property
    def intensity_org(self):
        return self._intensity_org
    @property
    def intensity(self):
        return self._intensity
    @property
    def intensity_z(self):
        return self._intensity_z
    @property
    def mass_ppm_error(self):
        return self._mass_ppm_error
    @property
    def dt(self):
        return self._dt
    @property
    def num_isotopes(self):
        return self._num_isotopes
    @property
    def mppid(self):
        return self._mppid

    def compute(self,
                drift_tube_length=90.33,
                neutral_mass=28.013):
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
        E = Vcell / drift_tube_length
        inv_E = 1.0 / (E * 100.0)
        # arrival_time (ms): array
        arrival_sec = np.array(self.arrival_time) / 1000
        # neutral_mass = 28.013 (N2 by default)
        # ========================
        # constant parameters
        # ========================
        # 1.60217657E-19 or 1.6021766208E-19
        e = 1.6021766208E-19
        charge_state = self.charge_state
        boltzmann_constant = 1.38064852E-23
        N0 = 101325/boltzmann_constant/273.15 # N0_(m-3)
        # ========================
        # computed parameters by given
        # ========================
        # P/V = P(torr) / Vcell
        self._p_v = P_torr / Vcell
        # E/N (Td) = E / P(torr) / 0.3535
        E_N = (E / P_torr) / 0.3535
        mass_in_kg = self.mz * self.charge_state * 1.66054E-27
        neutral_mass_in_kg = neutral_mass * 1.66054E-27
        reduced_mass_in_kg = (mass_in_kg * neutral_mass_in_kg / (mass_in_kg + neutral_mass_in_kg))
        # ========================

        slope, intercept, r_value, p_value, std_err = linregress(self._p_v, arrival_sec)
        # drift_time (sec) = arrival_sec - intercept
        drift_time = arrival_sec - intercept

        # compute CCS by Mason-Schamp Equation
        # ccs = 3 * e / 16 / N0 * np.sqrt(2 * np.pi / reduced_mass_in_kg / boltzmann_constant / T_K) \
        # * drift_time * 760 * T_K * Vcell / (drift_tube_length / 100)**2 / P_torr / 273.15 * 1E20
        K0 = drift_tube_length * drift_tube_length / slope * 273.15 / 760 / np.mean(T_K)
        ccs = 3 * charge_state * e / 16 / N0 / K0 / 0.0001 * np.sqrt(2 * np.pi / (boltzmann_constant * reduced_mass_in_kg * np.mean(T_K))) * 1e20
        properties = {'slope': slope, 'intercept': intercept, 'r2': r_value**2, 'p_value':p_value, 'k0':K0, 'ccs':ccs}
        for p in properties: self._metadata[p] = properties[p]
    
    def to_dict(self):
        return self._metadata
