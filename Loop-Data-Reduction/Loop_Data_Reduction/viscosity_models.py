import numpy as np

"""
Several viscosity models are available. See - Ye, X.; Kandlikar, S. G.; Li, C. Viscosity of Nanofluids Containing
Anisotropic Particles: A Critical Review and a Comprehensive Model. Eur. Phys. J. E 2019, 42 (12),
60–65. https://doi.org/10.1140/epje/i2019-11923-7.

Typical values of nu (intrinsic viscoirt), phi_max are found page 25 of the above reference and Table 4 page 24
type = ["Einstein","Simha rods","Simha ellipsoids","Simha disks", "Mooney","Khun","Brenner Cardiff","Jeffrey Acrivos",
"Batchelor Sphere","Batchelor Higgins","K-D rods","K-D discs","Modified K-D", "Modified M-P","Kandlikar Dilute"]
:return:  visc viscosity enhancement


"""

class viscosity(object):
    def Einstein(self,phi,phi_e,d_p,p):
        phi *=phi_e
        return 1 + 2.5 * phi

    def SimhaRods(self,phi,d_p,p=10):
        # valid for p>>1
        lam = 1.5
        nu = 8/5+ p**2 * (15*np.ln(2*p)-lam)**-1 +p**2 * (5*np.ln(2*p)-lam+1)**-1
        visc = 1+nu*phi
        return visc
    def SimhaEllipsoids(self,phi,d_p,p):
        # valid for p>>1
        lam = 1.8
        nu = 8/5+ p**2 * (15*np.ln(2*p)-lam)**-1 +p**2 * (5*np.ln(2*p)-lam+1)**-1
        visc = 1+nu*phi
        return visc
    def Simhadisks(self,phi,d_p,p):
        # valid for p>>1
        lam = 1.5
        visc = 16/5 *np.tanh(p**-1)*p**-1
        return visc
    def Mooney(self,phi,nu,phi_max,d_p,p):
        # phi_max is the maximum particle packing fraction - obtained by fitting
        visc = np.exp(nu*phi*((1-phi)*phi_max**-1)**-1)
        return visc
    def Khun(self,phi,d_p,p):
        if p < 1:
            nu = 2.5 +32*np.divide(1-p,15*np.pi*p)-.628*np.divide(1-p,1-0.075*p)
        elif p >= 1 and p <15:
            nu = 2.5 +0.4075*(p-1)**1.508
        elif p >=15  and p <40:
            nu = 8 / 5 + p ** 2 * (15 * np.ln(2 * p) - 1.5) ** -1 + p ** 2 * (
                        5 * np.ln(2 * p) - .5) ** -1
        elif p >=40 :
            nu = 2+ np.divide(.312*p-.5,np.log(2*p)-1.5)-1.872*p**-1
        visc = 1+nu*phi
        return visc
    def BrennerCardiff(self,phi,d_p,p):
        nu = 2+ np.divide(.312*p-.5,np.log(2*p)-1.5)-1.872*p**-1
        visc = 1 + nu* phi
        return visc
    def JeffreyAcrivos(self,phi,d_p,p):
        visc = 3+4/3* np.multiply(phi , np.divide(p,np.log(np.pi/phi)))
        return visc
    def BatchelorSphere(self,phi,phi_e,d_p,p):
        # good for phi <=.1
        phi *=phi_e
        visc = 1 + 2.5 * phi + 6.2 * phi ** 2.0
        return visc
    def BatchelorHiggins(self,phi,nu,k_h,d_p,p):
        # k_h is higgins constnt, 1 for sphreres, .4 for rods
        # nu is the intrinsic viscosity
        visc = 1 + np.multiply(nu, phi) + np.multiply(k_h,(np.multiply(nu, phi)) ** 2.0)
        return visc
    def KDrods(self,phi,p,d_p,phi_max):
        nu = .07*p**(5/3)
        visc = np.power(1-np.multiply(phi,phi_max**-1), -np.multiply(nu,phi_max))
        return visc
    def KDdiscs(self,phi,p,d_p,phi_max):
        nu = .3*p**(-1)
        visc = np.power(1-np.multiply(phi,phi_max**-1), -np.multiply(nu,phi_max))
        return visc
    def ModifiedKD(self,phi,p,d_a,d_p,phi_max,D):
        # D is fraction index of aggregates reasonable values of 1.6 to 2.5 for spherical  and 1.5 tl 2.45 for rod-like
        # d_a is effective aggregate diameter and d_p is single particle
        # phi a is the effective volume fraction of aggregates (but does this not depend on the preparation)
        nu = 2+ np.divide(.312*p-.5,np.log(2*p)-1.5)-1.872*p**-1
        phi_a = np.multiply(phi,np.power(d_a*d_p**-1,3-D))
        visc = np.power(1-np.multiply(phi_a,phi_max**-1), -np.multiply(nu,phi_max))
        return visc
    def ModifiedMP(self,phi,p,d_a,d_p,D):
        phi_a = np.multiply(phi,np.power(d_a*d_p**-1,3-D))
        phi_max = 2*(.321*p+3.02)**-1
        visc = np.power(1-np.multiply(phi_a,phi_max**-1),-2)
        return visc
    def ModifiedMP2(self,phi,p,d_a,d_p,D,phi_max):
        phi_a = np.multiply(phi,np.power(d_a*d_p**-1,3-D))
        visc = np.power(1-np.multiply(phi_a,phi_max**-1),-2)
        return visc
    def KandlikarDilute(self,phi,p,d_a,d_p,phi_max,D):
        if p < 1:
            nu = 2.5 +32*np.divide(1-p,15*np.pi*p)-.628*np.divide(1-p,1-0.075*p)
        elif p >= 1 and p <15:
            nu = 2.5 +0.4075*(p-1)**1.508
        elif p >=15  and p <40:
            nu = 8 / 5 + p ** 2 * (15 * np.ln(2 * p) - 1.5) ** -1 + p ** 2 * (
                        5 * np.ln(2 * p) - .5) ** -1
        elif p >=40 :
            nu = 2+ np.divide(.312*p-.5,np.log(2*p)-1.5)-1.872*p**-1
        phi_a = np.multiply(phi,np.power(d_a*d_p**-1,3-D))
        visc = np.power(1-np.multiply(phi_a,phi_max**-1), -np.multiply(nu,phi_max))
        return visc


def VogelFletcherTammann(A,n,T0,T):

    return np.multiply(n,np.exp(np.divide(A,T-T0)))

