from Fluids import steam
import numpy as np
import matplotlib.pyplot as plt


class Otto():
    def __init__(self, ratio=8, v_1=0.02, t_1=540, p_1=1, t_3=3600):
        """
        Otto Cycle info goes here!
        """
        self.ratio = ratio
        self.v_1 = v_1
        self.v_2 = 0
        self.v_3 = 0
        self.v_4 = 0
        self.p_1 = p_1
        self.p_2 = 0
        self.p_3 = 0
        self.p_4 = 0
        self.t_1 = t_1
        self.t_2 = 0
        self.t_3 = t_3
        self.t_4 = 0
        self.s_1 = 0
        self.s_2 = 0
        self.s_3 = 0
        self.s_4 = 0
        self.state1 = None
        self.state2 = None
        self.state3 = None
        self.state4 = None
        self.u1 = None
        self.u2 = None
        self.u3 = None
        self.u4 = None
        self.compressionstrokework = 0
        self.heataddition = 0
        self.powerstrokework = 0
        self.heatrejection = 0
        self.cycleefficiency = 0
        self.k = 1.4  # standard cold air

    def states(self):
        # calculate the 4 states

        self.t_2 = self.t_1 * self.ratio ** (self.k - 1)
        self.p_2 = self.p_1 * self.ratio * (self.t_2 / self.t_1)
        self.v_2 = self.v_1 / self.ratio
        self.p_3 = self.p_2 * (self.t_3 / self.t_2)
        self.t_4 = self.t_3 * (1 / self.ratio) ** (self.k - 1)
        self.p_4 = self.p_3 * (1 / self.ratio) * (self.t_4 / self.t_3)
        self.v_3 = self.v_2
        self.v_4 = self.v_1
        self.state1 = steam(T=540, Pressure=self.p_1, v=0.02)
        self.state2 = steam(T=self.t_2, Pressure=self.p_2, v=self.v_2)
        self.state3 = steam(T=self.t_3, Pressure=self.p_3, v=self.v_3)
        self.state4 = steam(T=self.t_4, Pressure=self.p_4, v=self.v_4)
        self.compressionstrokework = self.state2.u - self.state1.u
        self.powerstrokework = self.state3.u - self.state4.u
        self.heataddition = self.state3.u - self.state2.u
        self.heatrejection = self.state3.u - self.state1.u
        self.cycleefficiency = (self.powerstrokework - self.compressionstrokework) / (
                    self.heataddition - self.heatrejection)
        print(
            'Point 1 (pressure, temp, volume): {:.2f} atm, {:.2f} R, {:.2f} ft^3'.format(self.p_1, self.t_1, self.v_1))
        print(
            'Point 2 (pressure, temp, volume): {:.2f} atm, {:.2f} R, {:.6f} ft^3'.format(self.p_2, self.t_2, self.v_2))
        print(
            'point 3 (pressure, temp, volume): {:.2f} atm, {:.2f} R, {:.6f} ft^3'.format(self.p_3, self.t_3, self.v_3))
        print(
            'point 4 (pressure, temp, volume): {:.2f} atm, {:.2f} R, {:.2f} ft^3'.format(self.p_4, self.t_4, self.v_4))
        print('Internal Energies: {:.2f}BTU/lb∙°F, {:.2f}BTU/lb∙°F, {:.2f}BTU/lb∙°F, {:.2f}BTU/lb∙°F'.format(
            self.state1.u, self.state2.u,
            self.state3.u, self.state4.u))
        print('compression stroke work: {:.2f} ft*lbf'.format(self.compressionstrokework))
        print('power stroke work: {:.2f} ft*lbf'.format(self.powerstrokework))
        print('heat addition: {:.2f}'.format(self.heataddition))
        print('heat rejection: {:.2f}'.format(self.heataddition))
        print('Cycle efficiency: {:.2f}% '.format(abs(self.cycleefficiency * 100)))
        print('Entropies: {:.2f} BTU/lb∙°F, {:.2f} BTU/lb∙°F, {:.2f} BTU/lb∙°F, {:.2f} BTU/lb∙°F'.format(self.state1.s,
                                                                                                         self.state3.s,
                                                                                                         self.state3.s,
                                                                                                         self.state1.s))

    def plot_cycle_PV(self):
        """
        This function will graph the rankine cycle from HW6 part 3 on a T-S diagram.
        The two halves of the main curve,SaturatedLiquidLine and SaturatedVaporLine, are from sat_water_table.txt file.
        Graph also includes isobars for plow and phigh (constructed from state objects above and steam objects in Steam_work.py).
        :return: none, just the graph
        """
        plt.xlim(0, .0205)  # sets limits on x
        plt.ylim(0, 60)  # sets-limits on y
        plt.plot(self.v_1, self.p_1, marker='.', markerfacecolor='blue', markeredgecolor='black',
                 markersize=15)
        plt.plot(self.v_2, self.p_2, marker='.', markerfacecolor='blue', markeredgecolor='black', markersize=15)
        plt.plot(self.v_3, self.p_3, marker='.', markerfacecolor='blue', markeredgecolor='black', markersize=15)
        plt.plot(self.v_4, self.p_4, marker='.', markerfacecolor='blue', markeredgecolor='black', markersize=15)
        plt.plot()
        x1vals = [self.v_1, self.v_2, self.v_3, self.v_4, self.v_1]
        y1vals = [self.p_1, self.p_2, self.p_3, self.p_4, self.p_1]
        plt.plot(x1vals, y1vals, color='black')
        # print(self.state1.s)
        plt.fill_between(x1vals, y1vals, facecolor='gray',
                         alpha=0.5)  # fills in the graph between phigh and plow isobars
        plt.xlabel(r'V $ (ft^3) $', fontsize=12)  # x-label
        plt.ylabel(r'P $ (bar) $', fontsize=12)  # y-label
        txt = 'Summary:'
        txt += '\n$\eta_{cycle} = $' + '{:0.2f}%'.format(self.cycleefficiency)
        txt += '\n${power stroke} = $' + '{:0.2f}'.format(self.powerstrokework) + r'  ${ft*lbf}$'
        txt += '\n$W_{compression stroke} = $' + '{:0.2f}'.format(self.compressionstrokework) + r'  ${ft*lbf}$'
        txt += '\n$W_{heat addition} = $' + '{:0.2f}'.format(self.heataddition) + r'  $\frac{BTU}{lb*F}$'
        txt += '\n$Q_{heat rejected} = $' + '{:0.2f}'.format(self.heatrejection) + r'  $\frac{BTU}{lb*F}$'
        plt.text(.01, 58, txt, ha='left', va='top')
        plt.title('P-v Diagram')
        plt.plot()
        plt.show()

        ts, ps, hfs, hgs, sfs, sgs, vfs, vgs = np.loadtxt('sat_water_table.txt', skiprows=1,
                                                          unpack=True)  # reads-values from water table txt file
        # tcol, hcol, scol, pcol = np.loadtxt('superheated_water_table.txt', skip-rows=1, unpack=True)
        plt.xlim(0, 2.2)  # sets limits on x
        plt.ylim(0, 3500)  # sets-limits on y
        plt.plot(sfs * 0.23885, ts * (5 / 9))  # plots the sat fluid entropy vs. Tsat
        plt.plot(sgs * 0.23885, ts * (5 / 9),  # kelvin
                 color='red')  # plots the sat vapor entropy vs. Tsat #resposible for the right half of the curve
        TC1 = self.t_1*(5/9)  # converted to Kelvin
        TC2 = self.t_2*(5/9)
        TC3 = self.t_3*(5/9)
        TC4 = self.t_4*(5/9)
        x1vals = [self.state1.s, self.state3.s, self.state3.s, self.state1.s, self.state1.s]
        y1vals = [TC1, TC2, TC3, TC4, TC1]
        plt.plot(x1vals, y1vals, color='black')  # plots the plow isobar
        plt.fill_between(x1vals, y1vals, facecolor='gray', alpha=0.5)
        plt.plot(self.state1.s, TC1, marker='o', markeredgecolor='k', markerfacecolor='w')
        plt.plot(self.state3.s, TC2, marker='o', markeredgecolor='k', markerfacecolor='w')
        plt.plot(self.state3.s, TC3, marker='o', markeredgecolor='k', markerfacecolor='w')
        plt.plot(self.state1.s, TC4, marker='o', markeredgecolor='k', markerfacecolor='w')
        plt.xlabel(r'S $BTU/lb∙°F$', fontsize=12)  # x-label
        plt.ylabel(r'T $\left(K\right)$', fontsize=12)  # y-label
        txt = 'Summary:'
        txt += '\n$\eta_{cycle} = $' + '{:0.2f}%'.format(self.cycleefficiency)
        txt += '\n${power stroke} = $' + '{:0.2f}'.format(self.powerstrokework) + r'  ${ft*lbf}$'
        txt += '\n$W_{compression stroke} = $' + '{:0.2f}'.format(self.compressionstrokework) + r'  ${ft*lbf}$'
        txt += '\n$W_{heat addition} = $' + '{:0.2f}'.format(self.heataddition) + r'  $\frac{BTU}{lb*F}$'
        txt += '\n$Q_{heat rejected} = $' + '{:0.2f}'.format(self.heatrejection) + r'  $\frac{BTU}{lb*F}$'
        plt.text(.05, 3450, txt, ha='left', va='top')
        plt.title('T-s diagram')
        plt.show()


def main():
    Otto1 = Otto(ratio=8, v_1=0.02, t_1=540, p_1=1, t_3=3600)  # instantiate Otto object
    energies = Otto1.states()
    y = Otto1.plot_cycle_PV()


if __name__ == "__main__":
    main()
