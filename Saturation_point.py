import numpy as np
from EoS import PengRobinson


class Saturation:
    """
    * Calculation of saturation points of binary mixture of R32+R125 by different equations of state.
    In this class the following methods are present:

        - Ki: Equilibrium constant
        - bubble_P: Bubble point pressure
        - dew_P: Dew point pressure
        - bubble_and_dew: Bubble and dew points pressure
        - graph: P-x-y diagram
    """

    def __int__(self):
        pass

    @staticmethod
    def Ki(P, T, R, Tc, Pc, w, Kij, x, y):

        Eos    = PengRobinson(P, T, R, Tc, Pc, w, Kij)
        Vol_l  = Eos.Volume(x, phase=1)
        Vol_v  = Eos.Volume(y, phase=0)
        phi_L  = Eos.phi(x, Vol_l)
        phi_V  = Eos.phi(y, Vol_v)
        Ki     = (np.array(phi_L) / np.array(phi_V))
        return Ki

    @classmethod
    def bubble_P(cls, T, R, Tc, Pc, w, Kij, P_sat, x, tol, n_ite):

        #Raoult law (Initial guess)
        P   = sum(P_sat*x)
        y   = (P_sat/P)*x
        Ki  = cls.Ki(P, T, R, Tc, Pc, w, Kij, x, y)
        ite = 0

        while(True):
            ite   += 1
            Ki_old = Ki
            P      = P*sum(Ki*x)
            y      = Ki*x
            y      = y/sum(y)
            Ki     = cls.Ki(P, T, R, Tc, Pc, w, Kij, x, y)
            error  = sum(abs(Ki-Ki_old))

            if (error<tol) or (ite>n_ite):
                return P,y

    @classmethod
    def dew_P(cls, T, R, Tc, Pc, w, Kij, P_sat, y, tol, n_ite):

        #Raoult law (Initial guess)
        P  = (sum(y/P_sat))**-1
        x  = y*P/P_sat
        Ki = cls.Ki(P, T, R, Tc, Pc, w, Kij, x, y)

        ite = 0
        while (True):
            ite   += 1
            Ki_old = Ki
            P      = P*(sum(y/Ki))**-1
            x      = y/Ki
            x      = x/sum(x)
            Ki     = cls.Ki(P, T, R, Tc, Pc, w, Kij, x, y)

            error = sum(abs(Ki - Ki_old))
            if (error<tol) or (ite>n_ite):
                return P, x

    @classmethod
    def bubble_and_dew(cls, T, R, Tc, Pc, w, Kij, P_sat, x, tol, n_ite):

        #Saturation pressure
        P_bubble, y_calc = cls.bubble_P(T, R, Tc, Pc, w, Kij, P_sat, x, tol, n_ite)
        P_dew, x_calc    = cls.dew_P(T, R, Tc, Pc, w, Kij, P_sat, x, tol, n_ite)
        return P_bubble, P_dew, x_calc, y_calc