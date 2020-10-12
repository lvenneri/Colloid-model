import numpy as np
import matplotlib.pyplot as plt




def v_nf(phi, type, par={}, effective_phi_factor=1):
    """
    Several viscosity models are available. See - Ye, X.; Kandlikar, S. G.; Li, C. Viscosity of Nanofluids Containing
    Anisotropic Particles: A Critical Review and a Comprehensive Model. Eur. Phys. J. E 2019, 42 (12),
    60–65. https://doi.org/10.1140/epje/i2019-11923-7.

    Typical values of nu (intrinsic viscoirt), phi_max are found page 25 of the above reference and Table 4 page 24
    type = ["Einstein","Simha rods","Simha ellipsoids","Simha disks", "Mooney","Khun","Brenner Cardiff","Jeffrey Acrivos",
    "Batchelor Sphere","Batchelor Higgins","K-D rods","K-D discs","Modified K-D", "Modified M-P","Kandlikar Dilute"]
    :return:  visc viscosity enhancement


    """
    phi = phi* effective_phi_factor

    if type == "Einstein":
        visc = 1 + 2.5 * phi

    if type == "Simha rods":
        # valid for p>>1
        lam = 1.5
        par['nu'] = 8/5+ par['p']**2 * (15*np.ln(2*par['p'])-lam)**-1 +par['p']**2 * (5*np.ln(2*par['p'])-lam+1)**-1
        visc = 1+par['nu']*phi

    if type == "Simha ellipsoids":
        # valid for p>>1
        lam = 1.8
        par['nu'] = 8/5+ par['p']**2 * (15*np.ln(2*par['p'])-lam)**-1 +par['p']**2 * (5*np.ln(2*par['p'])-lam+1)**-1
        visc = 1+par['nu']*phi

    if type == "Simha disks":
        # valid for p>>1
        lam = 1.5
        visc = 16/5 *np.tanh(par['p']**-1)*p**-1

    if type == "Mooney":
        # phi_max is the maximum particle packing fraction - obtained by fitting
        visc = np.exp(par['nu']*phi*((1-phi)*par['phi_max']**-1)**-1)

    if type == "Khun":
        if par['p'] < 1:
            par['nu'] = 2.5 +32*np.divide(1-par['p'],15*np.pi*par['p'])-.628*np.divide(1-par['p'],1-0.075*par['p'])
        elif par['p'] >= 1 and par['p'] <15:
            par['nu'] = 2.5 +0.4075*(par['p']-1)**1.508
        elif par['p'] >=15  and par['p'] <40:
            par['nu'] = 8 / 5 + par['p'] ** 2 * (15 * np.ln(2 * par['p']) - 1.5) ** -1 + par['p'] ** 2 * (
                        5 * np.ln(2 * par['p']) - .5) ** -1
        elif par['p'] >=40 :
            par['nu'] = 2+ np.divide(.312*par['p']-.5,np.log(2*par['p'])-1.5)-1.872*par['p']**-1
        visc = 1+par['nu']*phi

    if type == "Brenner Cardiff":
        par['nu'] = 2+ np.divide(.312*par['p']-.5,np.log(2*par['p'])-1.5)-1.872*par['p']**-1
        visc = 1 + par['nu']* phi

    if type == "Jeffrey Acrivos":
        visc = 3+4/3* np.multiply(phi , np.divide(par['p'],np.log(np.pi/phi)))

    if type == "Batchelor Sphere":
        # good for phi <=.1
        visc = 1 + 2.5 * phi + 6.2 * phi ** 2.0

    if type == "Batchelor Higgins": # disregards third order terms
        # k_h is higgins constnt, 1 for sphreres, .4 for rods
        # nu is the intrinsic viscosity
        visc = 1 + np.multiply(par['nu'], phi) + np.multiply(par['k_h'],(np.multiply(par['nu'], phi)) ** 2.0)

    if type == "K-D rods":
        par['nu'] = .07*par['p']**(5/3)
        visc = np.power(1-np.multiply(phi,par['phi_max']**-1), -np.multiply(par['nu'],par['phi_max']))

    if type == "K-D discs":
        par['nu'] = .3*par['p']**(-1)
        visc = np.power(1-np.multiply(phi,par['phi_max']**-1), -np.multiply(par['nu'],par['phi_max']))

    if type == "Modified K-D":
        # D is fraction index of aggregates reasonable values of 1.6 to 2.5 for spherical  and 1.5 tl 2.45 for rod-like
        # d_a is effective aggregate diameter and d_p is single particle
        # phi a is the effective volume fraction of aggregates (but does this not depend on the preparation)
        par['nu'] = 2+ np.divide(.312*par['p']-.5,np.log(2*par['p'])-1.5)-1.872*par['p']**-1
        par['phi_a'] = np.multiply(phi,np.power(par['d_a']*par['d_p']**-1,3-par['D']))
        visc = np.power(1-np.multiply(par['phi_a'],par['phi_max']**-1), -np.multiply(par['nu'],par['phi_max']))

    if type == "Modified M-P":
        par['phi_a'] = np.multiply(phi,np.power(par['d_a']*par['d_p']**-1,3-par['D']))
        par['phi_max'] = 2*(.321*par['p']+3.02)**-1
        visc = np.power(1-np.multiply(par['phi_a'],par['phi_max']**-1),-2)


    if type == "Modified M-P 2":
        par['phi_a'] = np.multiply(phi,np.power(par['d_a']*par['d_p']**-1,3-par['D']))
        visc = np.power(1-np.multiply(par['phi_a'],par['phi_max']**-1),-2)


    if type == "Kandlikar Dilute":
        if par['p'] < 1:
            par['nu'] = 2.5 +32*np.divide(1-par['p'],15*np.pi*par['p'])-.628*np.divide(1-par['p'],1-0.075*par['p'])
        elif par['p'] >= 1 and par['p'] <15:
            par['nu'] = 2.5 +0.4075*(par['p']-1)**1.508
        elif par['p'] >=15  and par['p'] <40:
            par['nu'] = 8 / 5 + par['p'] ** 2 * (15 * np.ln(2 * par['p']) - 1.5) ** -1 + par['p'] ** 2 * (
                        5 * np.ln(2 * par['p']) - .5) ** -1
        elif par['p'] >=40 :
            par['nu'] = 2+ np.divide(.312*par['p']-.5,np.log(2*par['p'])-1.5)-1.872*par['p']**-1
        par['phi_a'] = np.multiply(phi,np.power(par['d_a']*par['d_p']**-1,3-par['D']))
        visc = np.power(1-np.multiply(par['phi_a'],par['phi_max']**-1), -np.multiply(par['nu'],par['phi_max']))



    return visc



def viscosity_tempdep():

    return

if __name__ == "__main__":
    print(v_nf(.1,"Einstein"))
    print(v_nf(.1,"Batchelor Sphere"))
    print(v_nf(.1,"Batchelor Higgins",{"k_h":.4,"nu":3.8},))

    phi = np.linspace(0,.1,50)
    nu_local = 10
    phi_max_local = .2
    k_h_local = .4
    p_local = 10
    d_a_local = 30
    d_p_local = 10
    D_local = 2
    phi_a_local = np.power(d_a_local*d_p_local**-1,3-D_local)

    if 1==1:
        plt.plot(phi,v_nf(phi,"Einstein"),label="Einstein")
        plt.plot(phi,v_nf(phi,"Mooney",{"phi_max":phi_max_local,"nu":nu_local}),label="Mooney")
        plt.plot(phi,v_nf(phi,"Batchelor Sphere"),label="Batchelor Sphere")
        plt.plot(phi,v_nf(phi,"Batchelor Sphere",effective_phi_factor=phi_a_local),label="Batchelor Sphere $\phi_e = $"+str(phi_a_local))
        plt.plot(phi,v_nf(phi,"Batchelor Higgins",{"k_h":k_h_local,"nu":nu_local}),label="Batchelor Higgins Rod")
        plt.plot(phi,v_nf(phi,"K-D rods",{"phi_max":phi_max_local,"p":p_local}),label="K-D rods")
        plt.plot(phi,v_nf(phi,"K-D discs",{"phi_max":phi_max_local,"p":p_local}),label="K-D discs")
        plt.plot(phi,v_nf(phi,"Modified K-D",{"phi_max":phi_max_local,"p":p_local,"d_a":d_a_local,"d_p":d_p_local,"D":D_local}),label="M-K-D")
        plt.plot(phi,v_nf(phi,"Modified M-P",{"p":p_local,"d_a":d_a_local,"d_p":d_p_local,"D":D_local}),label="Modified M-P")
        plt.plot(phi,v_nf(phi,"Kandlikar Dilute",{"phi_max":phi_max_local,"p":p_local,"d_a":d_a_local,"d_p":d_p_local,"D":D_local}),label="Kandlikar Dilute")
        plt.xlabel('$\phi$')
        plt.ylabel('$\mu_{R}$')
        # plt.plot(phi,v_nf(phi,"Einstein"),label="Einstein")
        plt.yscale('log')
        plt.legend()
    if 1==1:
        plt.figure()
        D_locals = np.linspace(1,2.5,5)
        plt.subplot(121)
        for D_local_i in D_locals:
            plt.plot(phi,v_nf(phi,"Modified M-P 2",{"p":p_local,"phi_max":phi_max_local,"d_a":d_a_local,"d_p":d_p_local,"D":D_local_i}),label="D="+str(D_local_i))
        plt.title("Vary D")
        plt.yscale('log')
        plt.legend()
        plt.xlabel('$\phi$')
        plt.ylabel('$\mu_{R}$')
        phi_max_locals = np.linspace(.01,1,5)
        plt.subplot(122)
        for phi_max_local_i in phi_max_locals:
            plt.plot(phi,v_nf(phi,"Modified M-P 2",{"p":p_local,"phi_max":phi_max_local_i,"d_a":d_a_local,"d_p":d_p_local,"D":D_local}),label="$\phi_{max}$="+str(phi_max_local_i))
        plt.title("Vary $\phi_{max}$")
        plt.legend()
        plt.yscale('log')
        plt.xlabel('$\phi$')
        plt.ylabel('$\mu_{R}$')

    if 1 == 1:
        # essentially, when does Batchelor with effective volume fraction match Maron Pierce

        batcherlor_set_lubrizol_P1 = 1+2.5*(phi*9) + 6.5*(phi*9)**2
        batcherlor_set_lubrizol_P2 = 1+2.5*(phi*4.3) + 6.5*(phi*4.3)**2

        D_measured_P1 = 1.85
        D_measured_P2 = 2.45

        plt.figure()
        # plt.plot(phi,batcherlor_set_lubrizol_P1,c='b',linestyle='--',label='batcherlor_set_lubrizol_P1')
        plt.plot(phi,batcherlor_set_lubrizol_P2,c='r',linestyle='--',label='batcherlor_set_lubrizol_P2')
        plt.plot(phi, v_nf(phi, "Modified M-P 2",
                           {"p": 1, "phi_max": .785, "d_a": 600, "d_p": 50,
                            "D": D_measured_P2}),c='r',linestyle='-', label="M-M-P")
        plt.yscale('log')
        plt.xlabel('$\phi$')
        plt.ylabel('$\mu_{R}$')
        plt.legend()
        # D_a_locals = np.linspace(50,200,10)
        # plt.subplot(221)
        # for D_local_i in D_locals:
        #     plt.plot(phi, v_nf(phi, "Modified M-P 2",
        #                        {"p": p_local, "phi_max": phi_max_local, "d_a": d_a_local, "d_p": d_p_local,
        #                         "D": D_local_i}), label="D=" + str(D_local_i))
        # plt.title("Vary D")
        plt.legend()

    plt.show()