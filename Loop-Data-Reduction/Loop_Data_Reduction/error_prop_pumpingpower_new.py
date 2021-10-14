
import numpy as np
import matplotlib.pyplot as plt
import itertools
import pandas as pd
L_t = 6.510 / 100
D_t = 1.825 / 100
S_t = 2 / 1000
A_c = np.sqrt(3) / 4 * (D_t + S_t) ** 2 - 0.5 * np.pi * (D_t / 2) ** 2
Wetted_p = np.pi * D_t / 2
d_hp = 4 * A_c / Wetted_p

def setup_text_plots(fontsize=8, usetex=True, style='default'):
    """
    e.g. setup_text_plots(fontsize=14, usetex=True,style='default')
    :param fontsize:
    :param usetex:
    :param style:
    :return: setup for nice plotting
    """
    import matplotlib
    matplotlib.style.use(style)
    matplotlib.rcParams['savefig.dpi'] = 200
    matplotlib.rc('legend', fontsize=fontsize, handlelength=3)
    matplotlib.rc('axes', titlesize=fontsize)
    matplotlib.rc('axes', titlesize=fontsize)
    matplotlib.rc('axes', labelsize=fontsize)
    matplotlib.rc('xtick', labelsize=fontsize)
    matplotlib.rc('ytick', labelsize=fontsize)
    matplotlib.rc('text', usetex=usetex)
    matplotlib.rc('font', size=fontsize, family='serif',
                  style='normal', variant='normal',
                  stretch='normal', weight='normal')
    matplotlib.rcParams['grid.color'] = 'k'
    matplotlib.rcParams['axes.grid'] = True
    matplotlib.rcParams['grid.linestyle'] = ':'
    matplotlib.rcParams['grid.linewidth'] = 0.5
    matplotlib.rcParams['figure.figsize'] = [6.0, 6.0]
    matplotlib.rcParams['errorbar.capsize'] = 3
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'
    matplotlib.rcParams['xtick.top'] = True
    matplotlib.rcParams['ytick.right'] = True
    matplotlib.rcParams['image.interpolation'] = 'nearest'
    matplotlib.rcParams['image.resample'] = False
    matplotlib.rcParams['figure.dpi'] = 100
    matplotlib.rcParams['mathtext.fontset'] = 'cm'
    matplotlib.rcParams['mathtext.rm'] = 'serif'
    matplotlib.rc('lines', linewidth=1.0, color='k')
    matplotlib.rc('legend', fancybox=False, shadow=False, framealpha=1, borderpad=0, frameon=True)


setup_text_plots()


def function_to_exp(xs):
    return np.divide(np.multiply(xs[:,0],xs[:,2]),xs[:,1])

def single_exp(mean = [40,1000,0.00089000],errors_frac = [.002,.01,.05],population= 1000):
    """

    :param mean:  mass flow, density, viscosity
    :param errors_frac:
    :param population:
    :return:
    """

    errors = np.multiply(mean,errors_frac)

    made_up_data = np.ones((population,3))*np.tile(mean,[population,1])

    made_up_data = made_up_data + np.multiply(np.random.standard_normal(size=(np.shape(made_up_data)))
                                              ,np.tile(errors,[population,1]))
    newdata = function_to_exp(made_up_data)

    if 1==1:
        plt.subplot(121)
        plt.hist(made_up_data[:,0], bins=35,  histtype='step', color='k',  lw=.5)
        plt.hist(made_up_data[:,1], bins=35,  histtype='step', color='k',  lw=.5)
        plt.hist(made_up_data[:,2], bins=35,  histtype='step', color='k', lw=.5)
        plt.subplot(122)
        plt.hist(newdata, bins=35,  histtype='step', color='k',  lw=.5)
        plt.show()
        print('MonteCarlo estimate: ', np.nanstd(newdata)/np.nanmean(newdata))
        norm_errors = np.divide(errors,mean)

        alpha = 26*L_t/A_c/d_hp
        print('Dumbass estimate: ', np.power(2*(errors_frac[0]**2+errors_frac[1]**2)+errors_frac[2],.5))


    return np.nanstd(newdata)/np.nanmean(newdata),mean[0],mean[1],mean[2],errors_frac[0],errors_frac[1],errors_frac[2]

single_exp(mean = [40,1000,10],errors_frac = [.002,.01,.05],population= 1000)

if __name__ == "__main__":
    if 1==0:
        water_prop = {"name":"water","rho":997.04740, "k":0.60650, "c":4181.40000000,"mu":0.00089000}
        isopar_ZnOp1_prop = {"name":"isopar_ZnOp1_prop","rho":1272.40000, "k": 0.19950, "c": 1356.02562087,"mu":0.00428558}
        ipa_ZnOp1_prop = {"name":"ipa_ZnOp1_prop","rho":1267, "k": 0.18623, "c": 1648.68981847,"mu":0.00314880}
        if 1==0:
            pars = [['mass_flow', np.linspace(25, 25, 1)],
                    ['rho', np.linspace(500,2000, 1)],
                    ['mu', np.linspace(0.00428558, 0.00428558*5, 1)],
                    ['mu error', np.linspace(0.01, 0.1, 1)]]
        else:
            pars = [['mass_flow', np.linspace(5, 100, 10)],
                    ['rho', np.linspace(500, 2000, 5)],
                    ['mu', np.linspace(0.00089000, 0.00428558 * 2, 5)],
                    ['mu error', np.linspace(0.01, 0.1, 20)]]
        parameter = []
        possible_vals = []

        # gets us the entire set of measurement space. We also know the accuracy on each measurement device, so we can put uncertainties
        # we'll see magnitude of the outputs measured at the two extremes, and we'll know the precision of our measurements
        # good plot is the fluid volume (bad) vs measured quantity delta (good)
        cols = ["Pumpin Power error","mass flow", "density", "mu","mass flow error", "density error", "mu error"]
        if 1==1:

            for par in pars:
                parameter.append(par[0])
                possible_vals.append(par[1])
            par_tests = list(itertools.product(*possible_vals))
            print("Number of tests: ",len(par_tests))
            # print(par_tests)

            data = pd.DataFrame(columns=cols)
            for p,par in enumerate(par_tests):
                print(p)
                test_bank = []
                results = single_exp(mean = [par[0],par[1],par[2]],errors_frac = [.002,.01,par[3]],population= 1000)
                data.loc[len(data.index)] = results

                # print(results)
            # data.to_csv('experimental_design_space.csv',index=True)
            # print(data)
            plt.scatter(data["mu error"],data["Pumpin Power error"],s=3,alpha=.3,c='k')
            plt.xlabel('$\sigma_{mu}$')
            plt.ylabel('$\sigma_{\dot{W}}$')

            plt.show()

