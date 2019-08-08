#!/usr/bin/env python

"""Plot electrochemical Pourbaix diagram.

Author: Michal Bajdich; Raul A. Flores
"""

#| - Import Modules
import os

from pylab import fill_between
from pylab import plot
from pylab import legend
import matplotlib

import matplotlib.pyplot as plt
import numpy as np
#__|

def plot_Pourbaix(
    surfs=None,

    # Gas References **********************************************************
    # #########################################################################
    h2=-6.77014123,
    zpe_h2=0.27283,  # exp. NIST 0.27283
    ts_h2=0.41,  # exp.
    cv_h2=0.,

    h2o=-14.21744725,
    zpe_h2o=0.5741,  # exp. NIST 0.5584250
    ts_h2o=0.67,  # 1 bar pressure
    cv_h2o=0.,


    # Adsorbates **************************************************************
    # #########################################################################
    zpe_o=0.065,  # 0.113 Max. G IrO2
    ts_o=0.,
    cv_o=0.,

    zpe_oh=0.37,  # 0.397 Max. G IrO2
    ts_oh=0.,
    cv_oh=0.,

    zpe_ooh=0.44,  # stays the same
    ts_ooh=0.,
    cv_ooh=0.,


    Umin=0.0,
    Umax=2.2,
    print_out=True,
    save_dir=None,
    file_name=None,
    close_plt=True,
    ):
    """
    """
    #| - plot_Pourbaix

    #| - Inputs
    if surfs is None:
        surfs = [
            [-178.96067932, 0, 0, 0],  # clean
            [-197.44772242, 0, 4, 0],  # 4O
            [-219.47269277, 0, 0, 4],  # 4OH*
            ]
    #__|

    #| - Setup
    # with golden ration and the whole shebang ...
    # settings size and font for revtex stylesheet
    fig_width_pt = 1.8 * 246.0  # Get from LaTeX using \showthe\columnwidth
    # fig_width_pt *= 300./72 # convert to 300 dpi
    inches_per_pt = 1.0 / 72.27               # Convert pt to inches
    # inches_per_pt = 1.0/300               # Convert pt to inches
    golden_mean = (np.sqrt(5) - 1.0) / 2.0         # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches
    fig_height = fig_width * golden_mean       # height in inches
    fig_size = [fig_width, fig_height]

    font_size = 10
    tick_font_size = 10

    # xlabel_pad = 8
    # ylabel_pad = 18

    matplotlib.rcParams['ps.usedistiller'] = 'xpdf'
    matplotlib.rcParams['font.family'] = 'sans-serif'
    #matplotlib.rcParams['font.family'] = 'serif'
    #matplotlib.rcParams['font.serif'] = 'Arial'
    matplotlib.rcParams['font.sans-serif'] = 'Arial'
    matplotlib.rcParams['font.size'] = 10
    matplotlib.rcParams['axes.labelsize'] = 2 * font_size
    matplotlib.rcParams['legend.fontsize'] = font_size
    matplotlib.rcParams['xtick.labelsize'] = tick_font_size
    matplotlib.rcParams['ytick.labelsize'] = tick_font_size
    matplotlib.rcParams['mathtext.default'] = 'regular'

    matplotlib.rcParams['lines.linewidth'] = 1.
    fig = plt.figure(figsize=fig_size, dpi=300)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])

    # Umin = 0.0
    # Umax = 2.2
    # Values for X,Y
    # pH = np.arange(0, 14, 0.10)
    # U = np.arange(Umin, Umax, 0.01)

    Umax2 = Umax + 0.06 * 14
    U2 = np.arange(Umin, Umax2, 0.01)

    ax.axis([0, 14, Umin, Umax])
    ax.set_xlabel(r'pH')
    ax.set_ylabel(r'U/V')

    #plt.xlabel('$pH$',family = 'Helvetica', size='40')
    #plt.ylabel('$U/V$',family = 'Helvetica',size='40')

    # Constants
    kbt = 0.0256
    const = kbt * np.log(10)

    extraticks = [1.23]
    plt.yticks(list(plt.yticks()[0]) + extraticks)
    #__|

    #| - Methods

    def addO(x, y):
        """
        """
        #| - addO
        out = -(h2o - h2) - 2 * (y + x * const) + ds_o
        return(out)
        #__|

    def addOH(x, y):
        """
        """
        #| - addOH
        out = -(h2o - 0.5 * h2) - (y + x * const) + ds_oh
        return(out)
        #__|

    def addH(x, y):
        """
        """
        #| - addH
        out = -0.5 * h2 + 1 * (y + x * const) + ds_h
        return(out)
        #__|

    def dg(i, x, y):
        """Function to calculate DG"""
        #| - dg
        if(surfs[i][1] + surfs[i][2] + surfs[i][3] > 0 and x == 0 and y == 0):

            #| - Print
            if print_out:
                print(
                    (
                        surfs[i][0] - surfs[0][0] + surfs[i][1] * addH(x, y) +
                        surfs[i][2] * addO(x, y) +
                        surfs[i][3] * addOH(x, y)
                        ) / (surfs[i][1] + surfs[i][2] + surfs[i][3])
                    )
            #__|

        out = 0. + \
            surfs[i][0] - surfs[0][0] + \
            surfs[i][1] * addH(x, y) + \
            surfs[i][2] * addO(x, y) + \
            surfs[i][3] * addOH(x, y) + \
            0.

        return(out)
        #__|

        # return surfs[i][0] - surfs[0][0] + surfs[i][1] * addH(x, y) +
        # surfs[i][2] * addO(x, y) + surfs[i][3] * addOH(x, y)
    #__|

    #| - Free Energy Corrections For Adsorbates

    #| - Previous Version
    # ds_o = zpe_o - (
    #     +(zpe_h2o - ts_h2o) + \
    #     -(zpe_h2 - ts_h2) + \
    #     0.
    #     )
    #
    # ds_oh = zpe_oh - (
    #     +(zpe_h2o - ts_h2o) + \
    #     -0.5 * (zpe_h2 - ts_h2) + \
    #     0.
    #     )
    #
    # sd_ooh = zpe_ooh - (
    #     +2 * (zpe_h2o - ts_h2o) + \
    #     -1.5 * (zpe_h2 - ts_h2) + \
    #     0.
    #     )
    #
    # ds_h = ds_oh - ds_o

    #__|


    #| - TEMP
    corr_h2o = zpe_h2o + cv_h2o - ts_h2o
    corr_h2 = zpe_h2 + cv_h2 - ts_h2

    ds_o = (zpe_o + cv_o - ts_o) - (
        +corr_h2o + \
        -corr_h2 + \
        0.
        )

    ds_oh = (zpe_oh + cv_oh - ts_oh) - (
        +corr_h2o + \
        -corr_h2 + \
        0.
        )

    sd_ooh = (zpe_ooh + cv_ooh - ts_ooh) - (
        +corr_h2o + \
        -corr_h2 + \
        0.
        )

    ds_h = ds_oh - ds_o
    #__|

    if print_out:
        print(ds_oh, ds_o, ds_h)
    #__|

    #| -  Find Intersects
    i = 0  # Ph=0
    lowest_surfaces = []
    for j in U2:
        # print i,j
        values = []

        for k in range(len(surfs)):
            # print k,dg(k,i,j)
            values.append(dg(k, i, j))
        sorted_values = sorted(range(len(values)), key=lambda k: values[k])
        lowest_surfaces.append(sorted_values[0])
        # print j, values
    # print lowest_surfaces
    #__|

    #| - Finding Crossover and Unique surface
    crossover = []
    uniquesurf = []
    uniquesurf.append(lowest_surfaces[0])
    old_value = lowest_surfaces[0]
    crossover.append(Umin)
    for j in range(len(U2)):
        if(lowest_surfaces[j] != old_value):
            uniquesurf.append(lowest_surfaces[j])
            crossover.append(U2[j])
            old_value = lowest_surfaces[j]

    crossover.append(Umax2)

    if print_out:
        print(crossover)
        print(uniquesurf)
    #__|

    #| - Main Plotting Loop
    color = [
        "turquoise",
        "green",
        "red",
        "blue",
        "gray",
        "gold",
        # "gray20",
        "purple",
        "orange",
        "magenta",

        "cyan",
        "yellow"
        ]

    pH2 = np.arange(0, 14, 0.01)
    for i in range(len(uniquesurf)):
        k = uniquesurf[i]
        foo = r"S$_{%i}$(H-%i O-%i OH-%i)" % (
            k,
            surfs[k][1],
            surfs[k][2],
            surfs[k][3],
            )

        #fbk = {'lw':0.0, 'edgecolor':color[i]}
        fbk = {'lw': 0.5, 'edgecolor': 'black'}

        fill_between(
            pH2,
            crossover[i] - pH2 * const,
            crossover[i + 1] - pH2 * const,
            facecolor=color[i],
            alpha=0.3,
            **fbk,
            )

        plot([], [], color=color[i], alpha=0.3, linewidth=5, label=foo)
    #__|

    #| - Additional Ploting and Plot Settings
    Vover = 0.526
    y = 1.23 + Vover - pH2 * const
    llabel = '$\eta$ = ' + repr(Vover) + ' V at S$_{12}$'
    plot(pH2, y, '-', color='black', lw=1, dashes=(3, 1), label=llabel)
    plot(pH2, 1.23 - pH2 * const, '--', color='blue', lw=1, dashes=(3, 1))

    ax.text(
        0.2,
        1.25,
        r'2H$_2$O $\leftrightarrow$ 4H$^+$ +O$_2$+4e$^-$',
        color='blue',
        rotation=-11.6,
        fontsize='x-small',
        )

    # legend()
    legend(
        ncol=2,
        fancybox=True,
        shadow=True,
        fontsize='x-small',
        handlelength=3)

    # TEMP | Raul added this to make background white so that I can view the
    # figure in jupyter
    fig.patch.set_facecolor("white")
    fig.patch.set_alpha(0.9)

    # ax = fig.add_subplot(111)
    # ax.plot(range(10))
    # ax.patch.set_facecolor('red')
    # ax.patch.set_alpha(0.5)

    if file_name is None:
        file_name = "Pourbaix_diagram.pdf"

    if save_dir is None:
        save_dir = "Pourbaix_plots_tmp"

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    save_file = os.path.join(
        save_dir,
        file_name,
        )

    fig.savefig(
        save_file,
        # 'Pourbaix_simple_generic.pdf',
        bbox_inches='tight',
        )

    if close_plt:
        plt.close()
    #__|

    #__|
