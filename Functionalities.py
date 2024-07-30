import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import time
import mpltern
import psutil
from matplotlib.ticker import FormatStrFormatter
from scipy import optimize
from Tp_flash import TpFlash
from Saturation_point import Saturation


class Vol_fac:

    def __init__(self, ftol):
        self.ite_flash = 0
        self.ite_opt   = 0
        self.ftol      = ftol
        self.eps       = 1e-20


    def f_obj(self, P, T, R, Tc, Pc, w, Ki, Kij, x, V_exp):

        V_calc,_,_,ite,_ = TpFlash().flash(P, T, R, Tc, Pc, w, Ki, Kij, x)
        self.ite_flash += ite
        f_obj  = abs(V_calc-V_exp)/V_exp
        return f_obj


    def P_adjustment(self, x0, T, R, Tc, Pc, w, Ki, Kij, P_sat, x, V_exp, tol, n_ite):

        args = (T, R, Tc, Pc, w, Ki, Kij, x, V_exp)
        P_bubble, P_dew, _, _ = Saturation().bubble_and_dew(T, R, Tc, Pc, w, Kij, P_sat, x, tol, n_ite)
        bounds = [(P_dew, P_bubble)]
        res      = optimize.minimize(self.f_obj, x0, args, method='SLSQP', bounds=tuple(bounds), options={'ftol':self.ftol, 'eps':self.eps})
        Pressure = res.x

        if x[0] != 0 and x[0] != 1.0:
            self.ite_opt = res.nit

        return Pressure, P_bubble, P_dew


    def output(self, x0, z, R, Tc, Pc, w, Ki, Kij, P_sat, beta_v, T, tol, num_ite):

        Pressure, P_bubble, P_dew = self.P_adjustment(x0, T, R, Tc, Pc, w, Ki, Kij, P_sat, z, beta_v,
                                                     tol, num_ite)

        print(f'* Pressure at volumetric fraction of {beta_v} and global molar fraction of {z} is {Pressure[0]:.2f} Pa')
        return Pressure


    def binary_Pxy(self, Components, R, Tc, Pc, w, Kij, P_sat, x_exp, y_exp, P_exp, beta_v, T, tol, n_ite, complete, save):

        legend_model = 'T\u03B2-flash (PR EoS)'

        if complete == 1:
            z1  = np.arange(0.0, 1.0+0.01, 0.01)
        else:
            z1 = np.arange(0.0, x_exp[-1]+0.01, 0.01)

        z2       = 1-z1
        z        = np.array([z1[0],z2[0]])
        P_b, P_d, _, _ = Saturation().bubble_and_dew(T, R, Tc, Pc, w, Kij, P_sat, z, tol, n_ite)
        x0       = (1-beta_v)*(P_b-P_d)+P_d
        Ki       = (Pc/x0)*np.exp(5.37*(1+w)*(1-(Tc/T)))
        P_calc   = np.zeros(len(z1))
        y_calc   = np.zeros(len(z1))
        x_calc   = np.zeros(len(z1))
        P_bubble = np.zeros(len(z1))
        P_dew    = np.zeros(len(z1))
        z        = []

        for i in range(len(z1)):
            z.append(np.array([z1[i], z2[i]]))
            P_bubble[i], P_dew[i], _, _ = Saturation().bubble_and_dew(T, R, Tc, Pc, w, Kij, P_sat, z[i], tol, n_ite)
            x0              = (1-beta_v)*(P_bubble[i]-P_dew[i])+P_dew[i]
            P_calc[i]       = self.output(x0, z[i], R, Tc, Pc, w, Ki, Kij, P_sat, beta_v, T, tol, n_ite)
            _,y,x,_, Ki_new = TpFlash().flash(P_calc[i], T, R, Tc, Pc, w, Ki, Kij, z[i])
            Ki              = Ki_new
            y_calc[i]       = y[0]
            x_calc[i]       = x[0]



        NAME = f"{Components[0] + '_' + Components[1] + '_' + str(T)}"
        fig, ax = plt.subplots(1, 1, figsize=(12.8, 8.8))

        ax.plot(x_exp, P_exp/1e6, 'ok', label='Experimental data', markersize=15)
        ax.plot(y_exp, P_exp/1e6, 'ok', markersize=15)
        ax.plot(x_calc, P_calc/1e6, 'k', label=legend_model, linewidth=3.0)
        ax.plot(y_calc, P_calc/1e6, 'k', linewidth=3.0)
        ax.set_xlabel(f"x\u2081,y\u2081", fontsize=30)
        ax.set_ylabel('Pressure [MPa]', fontsize=30)
        ax.tick_params(labelsize=30)
        plt.legend(loc='best', fontsize=22, frameon=False)
        plt.gca().spines['bottom'].set_linewidth(4)
        plt.gca().spines['left'].set_linewidth(4)
        plt.gca().spines['top'].set_linewidth(4)
        plt.gca().spines['right'].set_linewidth(4)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        if save == 1:
            plt.savefig(f"Output/{NAME}.png", dpi=300)

        plt.show()


    def ternary_plot(self, P, R, T, x_exp, y_exp, Tc, Pc, w, Kij, P_sat, beta_v, Components, tol, n_ite, save):

            data_storage_x_EoS = []
            data_storage_y_EoS = []
            data_storage_x_exp = []
            data_storage_y_exp = []

            for j, pressure in enumerate(P):
                P_calc         = np.zeros(len(x_exp[pressure]))
                P_b, P_d, _, _ = Saturation().bubble_and_dew(T, R, Tc, Pc, w, Kij, P_sat, x_exp[pressure][0], tol, n_ite)
                x0             = (1-beta_v)*(P_b-P_d)+P_d
                Ki             = (Pc/x0)*np.exp(5.37*(1+w)*(1-(Tc/T)))
                x1_output, x2_output, x3_output = [], [], []
                y1_output, y2_output, y3_output = [], [], []
                x1_exp, x2_exp, x3_exp          = [], [], []
                y1_exp, y2_exp, y3_exp          = [], [], []

                for i in range(len(x_exp[pressure])):
                    P_bubble, P_dew, _, _  = Saturation().bubble_and_dew(T, R, Tc, Pc, w, Kij, P_sat, x_exp[pressure][i], tol, n_ite)
                    x0                     = (1-beta_v)*(P_bubble-P_dew)+P_dew
                    P_calc[i]              = self.output(x0, x_exp[pressure][i], R, Tc, Pc, w, Ki, Kij, P_sat, beta_v, T, tol, n_ite)
                    _, y, x, _, Ki_new     = TpFlash().flash(P_calc[i], T, R, Tc, Pc, w, Ki, Kij, x_exp[pressure][i])
                    Ki                     = Ki_new

                    x1_output.append(x[0])
                    x2_output.append(x[1])
                    x3_output.append(x[2])

                    y1_output.append(y[0])
                    y2_output.append(y[1])
                    y3_output.append(y[2])

                # Organize the experimental points
                for i in range(len(x_exp[pressure])):
                    x1_exp.append(x_exp[pressure][i][0])
                    x2_exp.append(x_exp[pressure][i][1])
                    x3_exp.append(x_exp[pressure][i][2])

                    y1_exp.append(y_exp[pressure][i][0])
                    y2_exp.append(y_exp[pressure][i][1])
                    y3_exp.append(y_exp[pressure][i][2])

                data_x_EoS = {Components[0]: x1_output, Components[1]: x2_output, Components[2]: x3_output}
                data_x_exp = {Components[0]: x1_exp, Components[1]: x2_exp, Components[2]: x3_exp}
                data_y_EoS = {Components[0]: y1_output, Components[1]: y2_output, Components[2]: y3_output}
                data_y_exp = {Components[0]: y1_exp, Components[1]: y2_exp, Components[2]: y3_exp}

                data_storage_x_EoS.append(data_x_EoS)
                data_storage_y_EoS.append(data_y_EoS)
                data_storage_x_exp.append(data_x_exp)
                data_storage_y_exp.append(data_y_exp)

            # ----------------------------- Ternary diagram ---------------------------#
            fig = plt.figure(figsize=(12.8, 8.8))
            ax = fig.add_subplot(projection="ternary")
            for spine in ax.spines.values():
                spine.set_linewidth(2.5)

            ax.set_tlabel(Components[1], fontsize=28)
            ax.set_llabel(Components[0], fontsize=28)
            ax.set_rlabel(Components[2], fontsize=28)
            ax.taxis.set_label_position('tick1')
            ax.laxis.set_label_position('tick1')
            ax.raxis.set_label_position('tick1')
            ax.tick_params(labelsize=22)

            markers = 'o'
            colors  = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'gray', 'purple', 'cyan', 'orange']
            n = 0
            for i in range(len(data_storage_x_EoS)):
                data_x_EoS = data_storage_x_EoS[i]
                data_x_exp = data_storage_x_exp[i]
                data_y_EoS = data_storage_y_EoS[i]
                data_y_exp = data_storage_y_exp[i]
                ax.plot(data_x_EoS[Components[1]], data_x_EoS[Components[0]], data_x_EoS[Components[2]], colors[i], linewidth=3.0, label=f'{round(P[i]/10e5,2)} MPa')
                ax.plot(data_x_exp[Components[1]], data_x_exp[Components[0]], data_x_exp[Components[2]], colors[i] + markers, markerfacecolor='white', markersize=10)
                ax.plot(data_y_EoS[Components[1]], data_y_EoS[Components[0]], data_y_EoS[Components[2]], colors[i], linewidth=3.0)
                ax.plot(data_y_exp[Components[1]], data_y_exp[Components[0]], data_y_exp[Components[2]], colors[i] + markers, markerfacecolor='white', markersize=10)
                n += 1

            ax.grid()

            NAME = f"{Components[0] + '_' + Components[1] + '_' + Components[1] + '_' + '_' + str(T) + '_' + 'K'}"

            if save == 1:
                plt.savefig(f'Output/{NAME}.png', dpi=300)

            plt.show()


    def graph(self, P, R, T, x_exp, y_exp, Tc, Pc, w, Kij, P_sat, beta_v, Components, tol, n_ite, complete, save):

        if len(Components) == 3:
            self.ternary_plot(P, R, T, x_exp, y_exp, Tc, Pc, w, Kij, P_sat,
                              beta_v, Components, tol, n_ite, save)

        elif len(Components) == 2:
            self.binary_Pxy(Components, R, Tc, Pc, w, Kij, P_sat, x_exp, y_exp,
                                P, beta_v, T, tol, n_ite, complete, save)

        else:
            raise TypeError("This code only can plot diagrams for 2 or 3 components.")


    def AARD(self, Components, R, Tc, Pc, w, Kij, P_sat, x_exp, y_exp, P_exp, beta_v, T, tol, n_ite):

        if len(Components) == 2:
            if type(P_exp) is np.ndarray:
                z1 = x_exp
                z2 = 1 - z1
                z = np.array([z1[0], z2[0]])
                P_b, P_d, _, _ = Saturation().bubble_and_dew(T, R, Tc, Pc, w, Kij, P_sat, z, tol, n_ite)
                x0 = (1 - beta_v) * (P_b - P_d) + P_d
                Ki = (Pc / x0) * np.exp(5.37 * (1 + w) * (1 - (Tc / T)))
                P_calc = np.zeros(len(z1))
                y_calc = np.zeros(len(z1))
                x_calc = np.zeros(len(z1))
                P_bubble = np.zeros(len(z1))
                P_dew = np.zeros(len(z1))
                z = []

                for i in range(len(z1)):
                    z.append(np.array([z1[i], z2[i]]))
                    P_bubble[i], P_dew[i], _, _ = Saturation().bubble_and_dew(T, R, Tc, Pc, w, Kij, P_sat, z[i],
                                                                              tol, n_ite)
                    x0 = (1 - beta_v) * (P_bubble[i] - P_dew[i]) + P_dew[i]
                    P_calc[i] = self.output(x0, z[i], R, Tc, Pc, w, Ki, Kij, P_sat, beta_v, T, tol, n_ite)
                    _, y, x, _, Ki_new = TpFlash().flash(P_calc[i], T, R, Tc, Pc, w, Ki, Kij, z[i])
                    Ki = Ki_new
                    y_calc[i] = y[0]
                    x_calc[i] = x[0]

            elif type(T) is np.ndarray:
                T_exp = T
                z1 = x_exp
                z2 = 1 - z1
                z = np.array([z1[0], z2[0]])
                P_b, P_d, _, _ = Saturation().bubble_and_dew(T[0], R, Tc, Pc, w, Kij, P_sat[0], z, tol, n_ite)
                x0 = (1 - beta_v) * (P_b - P_d) + P_d
                Ki = (Pc / x0) * np.exp(5.37 * (1 + w) * (1 - (Tc / T[0])))
                P_calc = np.zeros(len(z1))
                y_calc = np.zeros(len(z1))
                x_calc = np.zeros(len(z1))
                P_bubble = np.zeros(len(z1))
                P_dew = np.zeros(len(z1))
                z = []
                for i in range(len(z1)):
                    z.append(np.array([z1[i], z2[i]]))
                    P_bubble[i], P_dew[i], _, _ = Saturation().bubble_and_dew(T[i], R, Tc, Pc, w, Kij, P_sat[i], z[i], tol, n_ite)
                    x0 = (1 - beta_v) * (P_bubble[i] - P_dew[i]) + P_dew[i]
                    P_calc[i] = self.output(x0, z[i], R, Tc, Pc, w, Ki, Kij, P_sat[i], beta_v, T_exp[i], tol, n_ite)
                    _, y, x, _, Ki_new = TpFlash().flash(P_calc[i], T_exp[i], R, Tc, Pc, w, Ki, Kij, z[i])
                    Ki = Ki_new
                    y_calc[i] = y[0]
                    x_calc[i] = x[0]

            AARD = (100/len(x_exp)) * sum(abs((P_calc - P_exp) / P_exp))

            print(f"AARD (%): {AARD:.3f}")

        elif len(Components) == 3:

            for j, pressure in enumerate(P_exp):
                P_calc = np.zeros(len(x_exp[pressure]))
                P_b, P_d, _, _ = Saturation().bubble_and_dew(T, R, Tc, Pc, w, Kij, P_sat,
                                                             x_exp[pressure][0], tol, n_ite)
                x0 = (1 - beta_v) * (P_b - P_d) + P_d
                Ki = (Pc / x0) * np.exp(5.37 * (1 + w) * (1 - (Tc / T)))

                for i in range(len(x_exp[pressure])):
                    P_bubble, P_dew, _, _ = Saturation().bubble_and_dew(T, R, Tc, Pc, w, Kij, P_sat,
                                                                        x_exp[pressure][i], tol, n_ite)
                    x0 = (1 - beta_v) * (P_bubble - P_dew) + P_dew
                    P_calc[i] = self.output(x0, x_exp[pressure][i], R, Tc, Pc, w, Ki, Kij, P_sat, beta_v,
                                            T, tol, n_ite)
                    _, y, x, _, Ki_new = TpFlash().flash(P_calc[i], T, R, Tc, Pc, w, Ki, Kij,
                                                         x_exp[pressure][i])
                    Ki = Ki_new

                AARD = (100/len(x_exp[pressure]))*sum(abs((P_calc-pressure)/pressure))
                print(f"AARD (%) for pressure {pressure/1e6:.2f} MPa: {AARD:.3f}")

        return 1


    def efficiency_report(self, Components, x_exp, R, Tc, Pc, w, Kij, P_sat, beta_v, T, tol, num_ite, save):

        z1          = np.arange(0.0, x_exp[-1]+0.01, 0.01)
        z2          = 1-z1
        Z           = np.array([z1[0], z2[0]])
        z           = []
        V_calc      = np.zeros(len(z1))
        P_calc      = np.zeros(len(z1))
        P_bubble    = np.zeros(len(z1))
        P_dew       = np.zeros(len(z1))
        n_ite_flash = np.zeros(len(z1))
        n_ite_opt   = np.zeros(len(z1))
        Time        = np.zeros(len(z1))
        x           = np.zeros(len(z1))
        y           = np.zeros(len(z1))
        cpu         = np.zeros(len(z1))
        memory      = np.zeros(len(z1))

        P_b, P_d, _, _  = Saturation().bubble_and_dew(T, R, Tc, Pc, w, Kij, P_sat, Z, tol, num_ite)
        x0 = (1-beta_v)*(P_b-P_d)+P_d
        Ki = (Pc/x0)*np.exp(5.37*(1+w)*(1-(Tc/T)))

        for i in range(len(z1)):
            z.append(np.array([z1[i], z2[i]]))
            P_bubble[i], P_dew[i], _, _ = Saturation().bubble_and_dew(T, R, Tc, Pc, w, Kij, P_sat, z[i], tol, num_ite)
            x0 = (1-beta_v)*(P_bubble[i]-P_dew[i])+P_dew[i]

            time_init = time.time()
            start_cpu = psutil.cpu_percent()
            start_mem = psutil.virtual_memory().percent

            P_calc[i] = self.output(x0, z[i], R, Tc, Pc, w, Ki, Kij, P_sat, beta_v, T, tol, num_ite)
            end_cpu = psutil.cpu_percent()
            end_mem = psutil.virtual_memory().percent
            time_end  = time.time()
            n_ite_flash[i]  = self.ite_flash
            n_ite_opt[i]    = self.ite_opt
            self.ite_flash  = 0
            self.ite_opt    = 0
            Time[i]   = (time_end-time_init)

            cpu[i]    = abs(end_cpu-start_cpu)
            memory[i] = abs(end_mem-start_mem)

            V_calc[i], Y, X, _, Ki_new = TpFlash().flash(P_calc[i], T, R, Tc, Pc, w, Ki, Kij, z[i])
            Ki = Ki_new
            y[i] = Y[0]
            x[i] = X[0]

        NAME   = f"{Components[0]+'_'+Components[1]+'_'+'Beta_specified'+'_'+str(beta_v)+'_'+str(T)}"
        output = {'x': z1,
                  'Bubble_point': P_bubble,
                  'Dew_point': P_dew,
                  'Calculated_pressure': P_calc,
                  'Molar fraction liquid (x)': x,
                  'Molar fraction vapor (y)': y,
                  'Specified_beta': beta_v,
                  'Calculated_beta': V_calc,
                  'Number of iterations flash': n_ite_flash,
                  'Number of iterations optimization': n_ite_opt,
                  'Processing time': Time,
                  'CPU percentage': cpu,
                  'Memory percentage': memory
                 }

        df = pd.DataFrame(output)
        print("")
        print(df.to_markdown())

        if save == 1:
            df.to_csv(f'Output/{NAME}.csv', index=False)

        return df


    def efficiency_plot(self, Components, x_exp, y_exp, P_exp, R, Tc, Pc, w, Kij, P_sat, beta_v, T, tol, n_ite, save):

        df     = self.efficiency_report(Components, x_exp, R, Tc, Pc, w, Kij, P_sat, beta_v, T, tol, n_ite, save=0)
        fig1, ax1 = plt.subplots(1, 1, figsize=(12.8, 8.8))
        NAME1 = f"{'Efficency' + '_' + Components[0] + '_' + Components[1] + '_'+'_'+str(T)}"

        colors = df['Processing time']  # Usar os valores de y como intensidade de cor
        norm   = plt.Normalize(min(df['Processing time']), max(df['Processing time']))  # Normalizar os valores de y para o intervalo [0, 1]
        cmap   = 'turbo'
        ax1.scatter(df['Molar fraction liquid (x)'], df['Calculated_pressure']/1e6, marker='o', c=colors, cmap=cmap, norm=norm, s=200)
        ax1.scatter(df['Molar fraction vapor (y)'], df['Calculated_pressure']/1e6, marker='o', c=colors, cmap=cmap, norm=norm, s=200)
        ax1.plot(x_exp, P_exp / 1e6, 'ok', markerfacecolor='white', label='Experimental data', markersize=15)
        ax1.plot(y_exp, P_exp / 1e6, 'ok', markerfacecolor='white', markersize=15)
        ax1.set_xlabel(f"x\u2081,y\u2081", fontsize=28)
        ax1.set_ylabel('Pressure [MPa]', fontsize=28)
        ax1.tick_params(labelsize=30)
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: '{:.2f}'.format(val)))
        plt.legend(loc='upper left', fontsize=22, frameon=False)
        plt.gca().spines['bottom'].set_linewidth(4)
        plt.gca().spines['left'].set_linewidth(4)
        plt.gca().spines['top'].set_linewidth(4)
        plt.gca().spines['right'].set_linewidth(4)
        cbar = fig1.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax1)
        cbar.ax.tick_params(labelsize=30)
        cbar.set_label(label='Computational effort [s]', size=30)
        if save == 1:
            plt.savefig(f"Output/{NAME1}.png", dpi=300)

        fig2, ax2 = plt.subplots(1, 1, figsize=(12.8, 8.8))
        NAME2 = f"{'Number of iterations flash' + '_' + Components[0] + '_' + Components[1] + '_' +'_'+str(T)}"
        ax2.bar(df['x'],df['Number of iterations flash'], width=0.005)
        ax2.set_xlabel('x\u2081', fontsize=18)
        ax2.set_ylabel("Iterations inside the TP-flash", fontsize=18)
        ax2.tick_params(labelsize=16)
        if save == 1:
            plt.savefig(f"Output/{NAME2}.png", dpi=300)

        fig3, ax3 = plt.subplots(1, 1, figsize=(12.8, 8.8))
        NAME3 = f"{'Number of iterations optimization' + '_' + Components[0] + '_' + Components[1] + '_' + '_' + str(T)}"
        ax3.bar(df['x'], df['Number of iterations optimization'], width=0.005)
        ax3.set_xlabel('x\u2081', fontsize=18)
        ax3.set_ylabel("Iterations of optimization", fontsize=18)
        ax3.tick_params(labelsize=16)
        if save == 1:
            plt.savefig(f"Output/{NAME3}.png", dpi=300)

        plt.show()