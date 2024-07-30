import numpy as np

class PengRobinson:
    """
    *Calculation of properties using Peng-Robinson EoS.
    In this class the following methods are present:

        - kij: Binary interaction parameter
        - mixture_rule: Molar volume calculation for a specific phase
        - Volume: Molar volume calculation for a specific phase in a multi-component system
        - phi: Coefficient fugacity (phi) for a multi-component system

    """

    def __init__(self, P, T, R, Tc, Pc, w, Kij):
        self.P = float(P)  # ---------------------------------> Pressure (Pa)
        self.T = T  # ----------------------------------------> Temperature (K)
        self.R = R
        self.Tc = Tc
        self.Pc = Pc
        self.w = w

        # Global parameters of Peng-Robinson equation
        self.sigma = 1 + np.sqrt(2)
        self.epsilon = 1 - np.sqrt(2)
        self.Kij = Kij

    def kij(self):
        """
        Organization of binary interaction parameter
        """

        cont = 0
        row = len(self.w)
        column = len(self.w)
        matrix = np.zeros((row, column))
        Kij = []

        for i in range(row):
            for j in range(column):
                if j > i:
                    matrix[i, j] = self.Kij[cont]
                    cont += 1
                else:
                    matrix[i, j] = 0

        matrix_T = np.transpose(matrix)
        matrix_kij = matrix + matrix_T

        for i in range(row):
            for j in range(column):
                Kij.append(matrix_kij[i, j])

        return Kij

    def mixture_rule(self, x):
        """
        * Mixture attraction and repulsion coefficients by VdW rule mixture
        """
        cont = 0
        a = 0
        b = 0
        ai, bi, aij = [], [], []
        Kij = self.kij()

        for i in range(len(self.w)):
            Tr = self.T / self.Tc[i]
            Kpri = 0.37464 + (1.54226 * self.w[i]) - (0.26992 * (self.w[i] ** 2))
            alfai = (1 + Kpri * (1 - np.sqrt(Tr))) ** 2

            aci = 0.45724 * (self.R ** 2) * (self.Tc[i] ** 2) / self.Pc[i]
            ai.append(aci * alfai)
            bi.append(0.07780 * self.R * (self.Tc[i] / self.Pc[i]))
        for j in range(len(self.w)):
            b += x[j] * bi[j]
            for l in range(len(self.w)):
                aij.append(np.sqrt(ai[j] * ai[l]) * (1 - Kij[cont]))
                a += x[j] * x[l] * aij[cont]
                cont += 1

        return a, b, aij, bi

    def Volume(self, x, phase):

        a, b, _, _ = self.mixture_rule(x)
        c3 = self.P
        c2 = (self.P * self.epsilon * b) + (self.P * self.sigma * b) - (self.P * b) - (self.R * self.T)
        c1 = (self.P * self.epsilon * self.sigma * (b ** 2)) - (self.P * self.epsilon * (b ** 2)) - (
                self.P * self.sigma * (b ** 2)) - (self.R * self.T * self.epsilon * b) - (
                     self.R * self.T * self.sigma * b) + a
        c0 = -(self.P * self.epsilon * self.sigma * (b ** 3)) - (
                self.R * self.T * self.epsilon * self.sigma * (b ** 2)) - (a * b)

        Volumes = np.roots([c3, c2, c1, c0])

        if (phase == 0) and (Volumes[0].real > b) and (Volumes[0].real > 0):
            resp = Volumes[0].real
        elif (phase == 0) and (Volumes[2].real > b) and (Volumes[2].real > 0):
            resp = Volumes[2].real

        if (phase == 1) and (Volumes[2].real > b) and (Volumes[2].real > 0):
            resp = Volumes[2].real
        elif (phase == 1) and (Volumes[0].real > b) and (Volumes[0].real > 0):
            resp = Volumes[0].real

        if (Volumes[0].real < 0) and (Volumes[2].real < 0) and Volumes[1].real > 0:
            resp = Volumes[1].real

        # print(f"Z:{(self.P*resp)/(self.R*self.T)}")

        return resp

    def phi(self, x, Vol):
        """
        * Coefficient fugacity (phi) for a multi-component system
        """

        phi = np.zeros(len(x))
        da, q, beta = [], [], []
        cont = 0
        termo = 0

        a, b, aij, bi = self.mixture_rule(x)
        Z = (self.P * Vol) / (self.R * self.T)
        num = Vol + (self.sigma * b)
        den = Vol + (self.epsilon * b)
        I = (1 / (self.sigma - self.epsilon)) * np.log(num / den)

        for k in range(len(x)):
            for j in range(len(x)):
                termo += aij[cont] * x[j]
                cont += 1

            da.append(2 * termo - a)
            q.append(a / (b * self.R * self.T))
            beta.append((b * self.P) / (self.R * self.T))
            termo = 0

        db = bi

        for l in range(len(x)):
            dq_l = q[l] * (1 + (da[l] / a) - (db[l] / b))
            Pa = (db[l] / b) * (Z - 1)
            Pb = np.log(Z - beta[l])
            Pc = dq_l * I
            phi[l] = (np.exp(Pa - Pb - Pc))

        return phi