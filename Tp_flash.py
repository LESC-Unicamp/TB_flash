import numpy as np
from EoS import PengRobinson


class TpFlash:
    """
    PT flash calculation for two-phase systems, using Rachford-Rice equation.
    In this class the following methods are present:
      - Ki: Equilibrium constant
      - sub_liq: Sub-cooled liquid condition
      - super_vapour: superheated steam condition
      - root: Volume fraction in the vapor phase(V) calculation by Newton-Raphson method.
      - flash: PT flash calculation using Rachford-Rice plus successive substitution.
    """

    def __init__(self):
        pass

    @staticmethod
    def Ki(P, T, R, Tc, Pc, w, Kij, x, y):

        Eos   = PengRobinson(P, T, R, Tc, Pc, w, Kij)
        Vol_l = Eos.Volume(x, phase=1)
        Vol_v = Eos.Volume(y, phase=0)
        phi_L = Eos.phi(x, Vol_l)
        phi_V = Eos.phi(y, Vol_v)
        Ki    = (np.array(phi_L) / np.array(phi_V))

        return Ki, phi_L, phi_V

    @classmethod
    def sub_liq(cls, P, T, R, Tc, Pc, w, Kij, z, Ki, ite):

        x = z
        error = 1
        tol = 1e-7
        while error > tol:
            ite += 1
            y = (Ki*z)/sum(Ki*z)
            Ki_new, phi_L, phi_V = cls.Ki(P, T, R, Tc, Pc, w, Kij, x, y)
            error = sum(abs(Ki-Ki_new))
            Ki = Ki_new

            if error <= tol:
                return Ki, x, y, ite

    @classmethod
    def super_vapour(cls, P, T, R, Tc, Pc, w, Kij, z, Ki, ite):

        y     = z
        error = 1
        tol   = 1e-7
        while error > tol:
            ite += 1
            x = (z/Ki)/(sum(z/Ki))
            Ki_new, phi_L, phi_V = cls.Ki(P, T, R, Tc, Pc, w, Kij, x, y)
            error = sum(abs(Ki-Ki_new))
            Ki = Ki_new

            if error <= tol:
                return Ki, x, y, ite

    @staticmethod
    def g(z, Ki, V):
        num = z*(Ki-1)
        den = 1-V+(V*Ki)
        return sum(num/den)

    @classmethod
    def root(cls, z, Ki):

        #Initial variables
        error = 1.0
        tol = 1e-10
        Vmin = 0
        Vmax = 1.0
        cont = 0

        #Step II of Michelsen
        ki = np.array(list(filter(lambda x: x > 1, Ki)))
        beta_upper = []
        beta_lower = []
        for j, k in enumerate(ki):
            i = np.where(Ki==k)
            lower = ((k*z[i])-1)/(k-1)
            upper = (1-z[i])/(1-k)
            if lower[0]>0:
                beta_lower.append(lower)
            if upper[0]>0:
                beta_upper.append(upper)

        if len(beta_lower)>0:
            Vmin = max(beta_lower)

        if len(beta_upper)>0:
            Vmax = min(beta_upper)

        #Step III of Michelsen
        V0 = 0.5 * (Vmin + Vmax)
        f = cls.g(z, Ki, V0)
        if f > 0:
            Vmin = V0
        else:
            Vmax = V0

        #Loop among steps IV and V of Michelsen
        while abs(error)>tol:
            cont += 1

            #Step IV of Michelsen
            Fnum = z*(Ki-1)
            Fden = 1-V0+(V0*Ki)
            Fv = sum(Fnum/Fden)
            dFnum = z*((Ki-1)**2)
            dFden = (1-V0+(V0*Ki))**2
            dFv = -sum(dFnum/dFden)
            f = cls.g(z, Ki, V0)
            if f > 0:
                Vmin = V0
            if f < 0:
                Vmax = V0
            V = V0 - Fv/dFv

            #Step V of Michelsen
            error = abs((V-V0)/V)
            if V>Vmin and V<Vmax:
                V0 = V
            else:
                V0 = 0.5*(Vmin+Vmax)

            if cont > 600:
                break

        return V0

    @classmethod
    def flash(cls, P, T, R, Tc, Pc, w, Ki, Kij, x):

        #Initial variables
        z = x
        error  = 1.0
        tol    = 1e-10
        V      = 0
        ite    = 0

        #Flash calculation
        while error > tol:
            ite += 1
            g0 = cls.g(z,Ki,0)
            if g0 <= 0.0:
                Ki, x, y, ite = cls.sub_liq(P, T, R, Tc, Pc, w, Kij, z, Ki, ite)
                g0       = cls.g(z, Ki, 0)

            g1 = cls.g(z, Ki, 1)
            if g1 >= 0.0:
                Ki, x, y, ite = cls.super_vapour(P, T, R, Tc, Pc, w, Kij, z, Ki, ite)
                g1       = cls.g(z, Ki, 1)

            if g0>0.0 and g1<0.0:
                V = cls.root(z, Ki)
                x = z/(1-V+(V*Ki))
                y = (Ki*z)/(1-V+(V*Ki))

                if error <= tol:
                    F = 1 - V
                    return V

                else:
                    #Choosing the EoS for equilibrium factor calculation
                    Ki_new, phi_L, phi_V = cls.Ki(P, T, R, Tc, Pc, w, Kij, x, y)

                error = sum(abs(Ki-Ki_new))

                Ki    = Ki_new

                if ite > 60:
                    tol = 1e-6

            else:
                V = cls.root(z, Ki)
                break

        return V,y,x,ite,Ki