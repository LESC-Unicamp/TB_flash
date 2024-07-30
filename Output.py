from Functionalities import Vol_fac
from Experimental_points import Points
from Parameters import Parameters


class Output:

    def __init__(self):
        pass

    @staticmethod
    def system(sys, **kwargs):

        if sys == 0:
            T = kwargs['T']
            x, y, P = Points().R32_R125(T)
            Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, P_sat = Parameters().R32_R125(T)

        elif sys == 1:
            T = kwargs['T']
            x, y, P = Points().R32_R152(T)
            Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, P_sat = Parameters().R32_R152(T)

        elif sys == 2:
            T = kwargs['T']
            x, y, P = Points().Propane_Butane(T)
            Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, P_sat = Parameters().Propane_Butane(T)

        elif sys == 3:
            T = kwargs['T']
            x, y, P = Points().CO2_Butane(T)
            Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, P_sat = Parameters().CO2_Butane(T)

        elif sys == 4:
            T = kwargs['T']
            x, y, P = Points().Ethylene_CO2(T)
            Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, P_sat = Parameters().Ethylene_CO2(T)

        elif sys == 5:
            T = kwargs['T']
            x, y, P = Points().Methane_CO2(T)
            Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, P_sat = Parameters().Methane_CO2(T)

        elif sys == 6:
            T = kwargs['T']
            x, y, P = Points().Methane_Ethane_CO2(kwargs['T'], kwargs['P'])
            Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, P_sat = Parameters().Methane_Ethane_CO2(kwargs['T'])

        else:
            raise TypeError("The selected system doesn't exist.")

        return x, y, P, T, Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, P_sat

    @classmethod
    def results(cls, **kwargs):

        system = kwargs['system']
        tol    = kwargs['tol']
        n_ite  = kwargs['n_ite']
        save   = kwargs['save']
        beta_v = kwargs['beta_v']
        case   = kwargs['case']
        ftol   = kwargs['ftol']
        complete = kwargs['complete']

        x, y, P, T, Components, R, Tc, Pc, w, sigma, Ek, mi, Kij, P_sat = cls.system(system, **kwargs)

        if case == 0:
            Vol_fac(ftol).graph(P, R, T, x, y, Tc, Pc, w, Kij, P_sat, beta_v,
                                Components, tol, n_ite, complete, save)

        elif case == 1:
            Vol_fac(ftol).efficiency_plot(Components, x, y, P, R, Tc, Pc, w, Kij,
                                          P_sat, beta_v, T, tol, n_ite, save)

        elif case == 2:
            Vol_fac(ftol).efficiency_report(Components, x, R, Tc, Pc, w, Kij, P_sat,
                                            beta_v, T, tol, n_ite, save)

        elif case == 3:
            Vol_fac(ftol).AARD(Components, R, Tc, Pc, w, Kij, P_sat, x, y, P,
                               beta_v, T, tol, n_ite)
