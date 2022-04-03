import numpy as np
from scipy.interpolate import griddata


class steam():
    """
    The steam class is used to find thermodynamic properties of steam along an isobar.
    The Gibbs phase rule tells us we need two independent properties in order to find
    all the other thermodynamic properties.  Hence, the constructor requires press of
    the isobar and one other property.
    """

    def __init__(self, Pressure, T=None, h=None, u=None, v=None, s=None):
        """
        constructor for steam
        :param pressure: pressure in kPa
        :param T: Temperature in degrees C
        :param x: quality of steam x=1 is saturated vapor, x=0 is saturated liquid
        :param v: specific volume in m^3/kg
        :param h: specific enthalpy in kJ/kg
        :param s: specific entropy in kJ/(kg*K)
        :param name: a convenient identifier
        """
        self.T = T
        self.h = h
        self.p = Pressure
        self.u = 0
        self.v = v
        self.u1 = 0
        self.u2 = 0
        self.u3 = 0
        self.u4 = 0
        self.x = 0
        self.s = s  # entropy - kj/(kg*K)
        if T == None and h == None and Pressure == None and u == None and v == None:
            return
        else:
            self.calc()

    def calc(self):
        " Interpolating stuff. I like Grapes"
        ts, ps, hfs, hgs, sfs, sgs, vfs, vgs = np.loadtxt(fname='sat_water_table.txt', skiprows=1,
                                                          unpack=True)  # np.loadtxt to read the saturated properties
        tcol, hcol, scol, pcol = np.loadtxt(fname='superheated_water_table.txt', skiprows=1,
                                            unpack=True)  # use np.loadtxt to read the superheated properties
        tR, tF, h, p, u = np.loadtxt(fname='properties.txt', skiprows=1, unpack=True)

        Pbar = self.p*1.01325
        TC = (self.T - 491.67) / 1.8
        v = self.v / 35.3146667
        self.u = float(griddata(p, tR, Pbar))
        Tsat = float(griddata(ps, ts, Pbar))
        hf = float(griddata(ps, hfs, Pbar))
        hg = float(griddata(ps, hgs, Pbar))
        sf = float(griddata(ps, sfs, Pbar))
        sg = float(griddata(ps, sgs, Pbar))
        vf = float(griddata(ps, vfs, Pbar))
        vg = float(griddata(ps, vgs, Pbar))
        if self.x != None:  # manual interpolation
            self.region = 'Saturated'
            self.s = sf + self.x * (sg - sf)
        elif self.h != None:
                self.region = 'Saturated'
                self.T = Tsat
                self.s = sf + self.x * (sg - sf)
        self.s = self.s*0.23885
        return self.u, self.s


def main():
    Otto = steam(Pressure=1, T=540, v=0.02)  # not enough information to calculate
    Otto.calc()


if __name__ == "__main__":
    main()
