import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

def CalculateSolution(time, thermal_conductivity, r_outer):
    plt.matplotlib.rc('text', usetex = True)
    plt.matplotlib.rc('grid', linestyle = 'dotted')
    plt.matplotlib.rc('figure', figsize = (6.4,4.8)) # (width,height) inches

    #time_arr = [0.005, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8]
    r = np.linspace(0, r_outer, 1000)

    # Roots of the Bessel function
    roots = np.array(sp.jn_zeros(0, 100))

    J_0_x_arr = np.multiply.outer(r, roots) / r_outer

    J_0_arr = sp.jv(0, J_0_x_arr)
    J_1_arr = sp.jv(1, roots)

    #Calculate the resulting constant
    div = J_1_arr * roots
    res = J_0_arr/div

    T = thermal_conductivity * time / r_outer**2

    pre_sum = np.exp(-roots**2 * T) * res
    sum = np.sum(pre_sum, 1) * 2.0

    result = 1 - sum

    legend_list = []
    legend_list.append("t=" + str(time))
    plt.plot(r, result)

    #plt.xlim((0, 15))
    #plt.ylim((-0.5, 1.1))
    #plt.legend(('$t=0s$', '${J}_1(x)$', '${J}_2(x)$',
    #'${J}_3(x)$', '${J}_4(x)$', '${J}_5(x)$'), loc = 0)
    plt.legend(legend_list, loc=0)
    plt.xlabel('$r/r_{outer}$')
    plt.ylabel('${T}/T_0$')
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    #plt.savefig('example-04-fig.pdf')
    # save the data for later use by pgfplots
    #np.savetxt('example-04.txt',list(zip(x,sp.jv(0,x),sp.jv(1,x),sp.jv(2,x),
    #sp.jv(3,x),sp.jv(4,x),sp.jv(5,x))), fmt="% .10e")

    return result