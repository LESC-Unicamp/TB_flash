import numpy as np


class Parameters:

    def __init__(self):
        pass

    @staticmethod
    def R32_R125(T):
        """
        Reference:
        :return:
        """
        Components = np.array(['R32', 'R125a'])  # ----------------------------> Components
        R          = 8.314  # ------------------------------------------------> Gas constant (J*mol-1*K-1)
        Tc         = np.array([351.25, 339.33])  # ---------------------------> Critical temperature (K)
        Pc         = np.array([5782000, 3629000])  # -------------------------> Critical pressure (Pascal)
        w          = np.array([0.2769, 0.3035])  # ---------------------------> acentric factor (dimensionless)
        mi         = np.array([2.76400, 3.1030])  # ---------------------------> Number of segments
        sigma      = np.array([2.75400, 3.1270])  # ---------------------------> Segment diameter (A°)
        Ek         = np.array([173.158, 156.891])  # --------------------------> Segment energy/Boltzmann's constant (K)
        Kij        = np.array([1.78571434e-06*T**2-1.00464289e-03*T+1.41316077e-01])

        Tr        = T/Tc
        Pa_R32    = -7.46398401811*(1-Tr[0])
        Pb_R32    = 2.43444074405*(1-Tr[0])**1.55
        Pc_R32    = -1.94666444684*(1-Tr[0])**2
        Pd_R32    = -2.49913094754*(1-Tr[0])**4
        a         = np.array([Pa_R32, Pb_R32, Pc_R32, Pd_R32])
        b         = np.array([-8.2627, 23.446, -15.19])
        Psat_R32  = Pc[0]*np.exp((1/Tr[0])*(a[0]+a[1]+a[2]+a[3]))
        Psat_R125 = Pc[1]*np.exp((b[0]*Tr[1]**2)+(b[1]*Tr[1])+b[2])
        Psat      = np.array([Psat_R32, Psat_R125])

        return Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, Psat

    @staticmethod
    def R32_R152(T):

        Components = np.array(['R32', 'R152a'])  # ---------------> Components
        R = 8.314  # --------------------------------------------> Gas constant (J*mol-1*K-1)
        Tc = np.array([351.25, 386.41])  # -----------------------> Critical temperature (K)
        Pc = np.array([5782000, 4520000])  # ---------------------> Critical pressure (Pascal)
        w = np.array([0.2769, 0.2752])  # -----------------------> acentric factor (dimensionless)
        mi = np.array([2.76400, 2.7180])  # -----------------------> Number of segments
        sigma = np.array([2.75400, 3.1010])  # -----------------------> Segment diameter (A°)
        Ek = np.array([173.158, 192.117])  # ----------------------> Segment energy/Boltzmann's constant (K)
        Kij = [0.0329]

        a = np.array([4.26224, 821.092, -28.554])
        b = np.array([4.23406, 896.171, -34.714])
        Psat_R32 = a[0]-(a[1]/(T+a[2]))
        Psat_R152 = b[0]-(b[1]/(T+b[2]))
        Psat = 1e5*np.array([np.exp(Psat_R32), np.exp(Psat_R152)])

        return Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, Psat

    @staticmethod
    def Benzene_Toluene(T):

        Components = np.array(['Benzene', 'Toluene'])  # ----------> Components
        R = 8.314  # --------------------------------------------> Gas constant (J*mol-1*K-1)
        Tc = np.array([562, 591.75])  # --------------------------> Critical temperature (K)
        Pc = np.array([4895000, 4108000])  # --------------------> Critical pressure (Pascal)
        w = np.array([0.2103, 0.264])  # ----------------------> acentric factor (dimensionless)
        mi = np.array([3.0576, 2.5911])  # -----------------------> Number of segments
        sigma = np.array([3.7983, 3.7171])  # --------------------> Segment diameter (A°)
        Ek = np.array([236.77, 275.6736])  # ---------------------> Segment energy/Boltzmann's constant (K)
        Kij = np.array([0.00245])

        a = np.array([4.72183, 1660.652, -1.461])
        b = np.array([4.07827, 1343.943, -53.773])
        Psat = []
        for i, t in enumerate(T):
            Psat_C3H8 = a[0] - (a[1] / (t + a[2]))
            Psat_C4H10 = b[0] - (b[1] / (t + b[2]))
            Psat.append(1e5 * np.array([np.exp(Psat_C3H8), np.exp(Psat_C4H10)]))

        return Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, Psat

    @staticmethod
    def Propane_Butane(T):

        Components = np.array(['Propane', 'Butane'])  # ---------> Components
        R = 8.314  # ----------------------------------------> Gas constant (J*mol-1*K-1)
        Tc = np.array([369.8, 425.1])  # ---------------------> Critical temperature (K)
        Pc = np.array([4248000, 3796000])  # -----------------> Critical pressure (Pascal)
        w = np.array([0.152, 0.2])  # -----------------------> acentric factor (dimensionless)
        mi = np.array([2.002, 2.3316])  # ---------------------> Number of segments
        sigma = np.array([3.6184, 3.7086])  # --------------------> Segment diameter (A°)
        Ek = np.array([208.11, 222.88])  # --------------------> Segment energy/Boltzmann's constant (K)
        Kij = np.array([0.0008])


        a = np.array([3.98292, 819.296, -24.417])
        b = np.array([4.35576, 1175.581, -2.071])
        Psat_C3H8 = a[0] - (a[1] / (T + a[2]))
        Psat_C4H10 = b[0] - (b[1] / (T + b[2]))
        Psat = 1e5 * np.array([np.exp(Psat_C3H8), np.exp(Psat_C4H10)])

        return Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, Psat

    @staticmethod
    def CO2_Butane(T):

        Components = np.array(['CO\u2082', 'Butano'])  # ---------> Components
        R = 8.314  # --------------------------------------------> Gas constant (J*mol-1*K-1)
        Tc = np.array([304.13, 425.1])  # ----------------------> Critical temperature (K)
        Pc = np.array([7377300, 3796000])  # --------------------> Critical pressure (Pascal)
        w = np.array([0.2239, 0.2])  # -----------------------> acentric factor (dimensionless)
        mi = np.array([2.0729, 2.3316])  # -----------------------------> Number of segments
        sigma = np.array([2.7852, 3.7086])  # --------------------------> Temperature-independent segment diameter (A°)
        Ek = np.array([169.21, 222.88])  # -----------------------------> Depth of pair potential/Boltzmann constant (K)
        Kij = np.array([0.15])

        a = np.array([6.81228, 1301.679, -3.494])
        b = np.array([3.85002, 909.65, -36.146])
        Psat_CO2 = a[0] - (a[1] / (T + a[2]))
        Psat_C4H10 = b[0] - (b[1] / (T + b[2]))
        Psat = 1e5*np.array([np.exp(Psat_CO2),np.exp(Psat_C4H10)])

        return Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, Psat

    @staticmethod
    def Ethylene_CO2(T):

        Components = np.array(['Ethylene', 'CO\u2082'])  # ---------> Components
        R = 8.314  # --------------------------------------------> Gas constant (J*mol-1*K-1)
        Tc = np.array([282.34, 304.13])  # ----------------------> Critical temperature (K)
        Pc = np.array([50.41e5, 7377300])  # --------------------> Critical pressure (Pascal)
        w = np.array([0.087, 0.2239])  # -----------------------> acentric factor (dimensionless)
        mi = np.array([1.5930, 2.0729])  # -----------------------------> Number of segments
        sigma = np.array([3.4450, 2.7852])  # --------------------------> Temperature-independent segment diameter (A°)
        Ek = np.array([176.47, 169.21])  # -----------------------------> Depth of pair potential/Boltzmann constant (K)
        Kij = [0.05420150963997733]


        a = np.array([3.87261, 584.146, -18.307])
        b = np.array([6.81228, 1301.679, -3.494])
        Psat_C3H8 = a[0] - (a[1] / (T + a[2]))
        Psat_C4H10 = b[0] - (b[1] / (T + b[2]))
        Psat = 1e5 * np.array([np.exp(Psat_C3H8), np.exp(Psat_C4H10)])

        return Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, Psat

    @staticmethod
    def Methane_CO2(T):
        Components = np.array(['Methane', 'CO\u2082'])  # ---------> Components
        R = 8.314  # --------------------------------------------> Gas constant (J*mol-1*K-1)
        Tc = np.array([190.6, 304.13])  # ----------------------> Critical temperature (K)
        Pc = np.array([4610000, 7377300])  # --------------------> Critical pressure (Pascal)
        w = np.array([0.01142, 0.2239])  # -----------------------> acentric factor (dimensionless)
        mi = np.array([1.0, 2.0729])  # -----------------------------> Number of segments
        sigma = np.array([3.7039, 2.7852])  # --------------------------> Temperature-independent segment diameter (A°)
        Ek = np.array([150.03, 169.21])  # -----------------------------> Depth of pair potential/Boltzmann constant (K)
        if T == 230:
            Kij = np.array([0.095])
        elif T == 270:
            Kij = np.array([0.096])

        a = np.array([4.22061, 516.689, 11.223])
        b = np.array([6.81228, 1301.679, -3.494])
        Psat_C3H8 = a[0] - (a[1] / (T + a[2]))
        Psat_C4H10 = b[0] - (b[1] / (T + b[2]))
        Psat = 1e5 * np.array([np.exp(Psat_C3H8), np.exp(Psat_C4H10)])

        return Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, Psat

    @staticmethod
    def Methane_Ethane_CO2(T):
        Components = np.array(['Methane', 'Ethane', 'Carbon dioxide'])  # ---------> Components
        R = 8.314  # --------------------------------------------> Gas constant (J*mol-1*K-1)
        Tc = np.array([190.6, 305.3, 304.13])  # ----------------------> Critical temperature (K)
        Pc = np.array([4610000, 4900000, 7377300])  # --------------------> Critical pressure (Pascal)
        w = np.array([0.01142, 0.0995, 0.2239])  # -----------------------> acentric factor (dimensionless)
        mi    = np.array([])  # -----------------------------> Number of segments
        sigma = np.array([])  # --------------------------> Temperature-independent segment diameter (A°)
        Ek    = np.array([])  # -----------------------------> Depth of pair potential/Boltzmann constant (K)
        Kij   = [-0.0026, 0.0919, 0.1322]  # -----------------------------> Binary interaction parameters (dimensionless)

        a         = np.array([4.22061, 516.689, 11.223])
        b         = np.array([3.93835, 659.739, -16.719])
        c         = np.array([6.81228, 1301.679, -3.494])
        Psat_CH4  = a[0]-(a[1]/(T+a[2]))
        Psat_C2H8 = b[0]-(b[1]/(T+b[2]))
        Psat_CO2  = c[0]-(c[1]/(T+c[2]))
        P_sat     = 1e5 * np.array([np.exp(Psat_CH4), np.exp(Psat_C2H8), np.exp(Psat_CO2)])
        return Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, P_sat