#/usr/bin/env python

"""Plot 2D, 3D volcano, and Scaling for OER systems.

Illustrate simple contour plotting, contours on an image with
a colorbar for the contours, and labelled contours.

See also contour_image.py.
"""

#| - Import Modules
import os
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

from pylab import polyfit
from pylab import poly1d
from pylab import plot
#__|

def create_OER_plots(
    data,
    annotation_placement_dict={},
    plots_folder="OER_plots",
    exp_lines=True,
    annotate_data=True,
    axes_ranges={
        "x1": 1.0,
        "x2": 2.0,
        "y1": 1.4,
        "y2": 2.0,
        },
    ):
    """I'm just wrapping Michals whole script in a method.

    Args:
        data:
        plots_folder:
    """
    #| - create_OER_plots
    calc_systems = data

    #| - Styling and Setup
    # settings size and font for revtex stylesheet

    # Get this from LaTeX using \showthe\columnwidth
    fig_width_pt = 1.8 * 246.0
    #fig_width_pt *= 300./72 # convert to 300 dpi
    inches_per_pt = 1.0 / 72.27               # Convert pt to inches
    #inches_per_pt = 1.0/300               # Convert pt to inches
    golden_mean = (np.sqrt(5) - 1.0) / 2.0         # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches
    fig_height = fig_width * golden_mean       # height in inches
    fig_size = [fig_width, fig_height]
    fig = plt.figure(figsize=fig_size, dpi=300)

    font_size = 9
    tick_font_size = 8
    xlabel_pad = 8
    ylabel_pad = 18
    matplotlib.rcParams['ps.usedistiller'] = 'xpdf'

    matplotlib.rcParams['font.size'] = 10
    #matplotlib.rcParams['axes.labelsize'] = 2*font_size
    matplotlib.rcParams['axes.labelsize'] = font_size
    matplotlib.rcParams['legend.fontsize'] = font_size
    matplotlib.rcParams['xtick.labelsize'] = tick_font_size
    matplotlib.rcParams['ytick.labelsize'] = tick_font_size

    font_default = 'helvetica'
    #font_default='cmss'

    def setfont(font=font_default, unicode=True):
        """Set font.

        Set Matplotlibs rcParams to use LaTeX for font rendering.
        Revert all changes by calling rcdefault() from matplotlib.

        Parameters:
        -----------
        font: string
            "Helvetica"
            "Times"
            "Computer Modern"

        usetex: Boolean
            Use unicode. Default: False.

        """
        #| - setfont
        # Use TeX for all figure text!
        plt.rc('text', usetex=True)

        font = font.lower().replace(" ", "")
        if font == 'times':
            # Times
            font = {'family': 'serif', 'serif': ['Times']}
            preamble = r"""
                          \usepackage{color}
                          \usepackage{mathptmx}
                       """
        elif font == 'helvetica':
            # Helvetica
            # set serif, too. Otherwise setting to times and then
            # Helvetica causes an error.
            font = {'family': 'sans-serif', 'sans-serif': ['Helvetica'],
                    'serif': ['cm10']}
            preamble = r"""
                          \usepackage{color}
                          \usepackage[tx]{sfmath}
                          \usepackage{helvet}
                          \usepackage{sansmath}
                       """
        else:
            # Computer modern serif
            font = {'family': 'serif', 'serif': ['cm10']}
            # preamble = r"""
            preamble = r"""
                        \usepackage{color}
                        """

        if font == 'cmss':
            # Computer modern sans serif
            font = {'family': 'sans-serif', 'serif': ['cmss']}
            preamble = r"""
                          \usepackage{color}
                          \usepackage[tx]{sfmath}
                       """

        if unicode:
            # Unicode for Tex
            #preamble =  r"""\usepackage[utf8]{inputenc}""" + preamble
            # inputenc should be set automatically
            plt.rcParams['text.latex.unicode'] = True

        # print font, preamble
        plt.rc('font', **font)
        plt.rcParams['text.latex.preamble'] = preamble
        #__|

    setfont(
        font_default,
        unicode=True,
        # unicode=False,
        )

    matplotlib.rcParams['lines.linewidth'] = 1.

    #matplotlib.rcParams['ytick.direction'] = 'out'
    #matplotlib.rcParams['xtick.direction'] = 'out'

    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])

    zoom = 0.5
    d1 = 3 * zoom
    d2 = 4 * zoom
    xcenter = 1.5  # 0.65
    #ycenter=1.23#2.4
    ycenter = 0.8  # 2.4

    x1 = xcenter - d1  # -0.6
    x2 = xcenter + d1  # 2.2
    y1 = ycenter - d2  # 1#0.5
    y2 = ycenter + d2  # 5
    ax.axis([x1, x2, y1, y2])
    ax.set_xlabel(r'$\Delta$G$_{\sf O}$ - $\Delta$G$_{\sf OH}$ (eV)')
    #ax.set_ylabel(r'$\Delta$G$_{\sf OOH}$ -$\Delta$G$_{\sf O}$ (eV)')
    ax.set_ylabel(r'$\Delta$G$_{\sf OH}$')

    delta = 0.025
    x = np.arange(x1, x2 + delta, delta)
    y = np.arange(y1, y2 + delta, delta)
    X, Y = np.meshgrid(x, y)

    #Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
    #Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
    # difference of Gaussians
    #Z = 10.0 * (Z2 - Z1)
    #__|

    #| - Methods
    #fit=[0.84527288, 3.38026638]
    def ooh_oh_scaling(doh):
        """ooh_oh_scaling equation."""
        #| - ooh_oh_scaling
        #like ambars
        #dooh=0.5*doh  + 3.0		 #O
        #normal one

        dooh = doh + 3.2
        return(dooh)
        #__|

    def overpotential(doh, do):
        """Calculate overpotential.

        Args:
            doh:
            do:
        """
        #| - overpotential
        dooh = ooh_oh_scaling(doh)
        dg14 = [doh, do - doh, dooh - do, - dooh + 4.92]
        m = max(dg14)
        return(m - 1.23)
        #return doh*do
        #__|

    def overpotential2(x, doh):
        """Calculate overpotential (version 2).

        Args:
            x:
            doh:
        """
        #| - overpotential2
        dooh = ooh_oh_scaling(doh)
        dg14 = [doh, x, -x + 2.46, -dooh + 4.92]
        m = max(dg14)
        return(m - 1.23)
        #return doh*do
        #__|

    def overpotential3(x, doh):
        """Calculate overpotential (version 3).

        Args:
            x:
            doh:
        """
        #| - overpotential3
        dooh = ooh_oh_scaling(doh)
        dg14 = [doh, x, dooh - (x + doh), -dooh + 4.92]
        m = max(dg14)
        return(m - 1.23)

        #return doh*do
        #__|

    def overpotential_label(doh, do):
        """Return overpotential label.

        Args:
            doh:
            do:
        """
        #| - overpotential_label
        dooh = ooh_oh_scaling(doh)
        dg14 = [doh, do - doh, dooh - do, -dooh + 4.92]
        m = max(dg14)
        for i in range(len(dg14)):
            if(m == dg14[0]):
                return(r'OH lim.')
            if(m == dg14[1]):
                return(r'OH-O lim.')
            if(m == dg14[2]):
                return(r'O-OOH lim.')
            if(m == dg14[3]):
                return( r'OOH-O$_{\sf 2}$ lim.')
        #return doh*do
        #__|

    #Z=overpotential(X,Y)
    #__|

    # *************************************************************************
    #| - OER_contour_plot *****************************************************
    Z = []
    for j in y:
        tmp = []
        for i in x:
            tmp.append(overpotential3(i, j))
        Z.append(tmp)


    #print overpotential(0.8,2.4)

    Z = np.array(Z)


    #im = plt.imshow(Z, origin='lower',interpolation='bilinear',
    #                cmap=cm.jet_r, extent=(x1,x2,y1,y2), vmin=0, vmax=2)

    origin = 'lower'
    levels = np.arange(0.0, 2, 0.1)
    #levels = np.arange(0.2, 2, 0.1)
    CS = plt.contourf(
        X,
        Y,
        Z,
        levels,
        #20,
        # [-1, -0.1, 0, 0.1],
        #alpha=0.8,
        #cmap=plt.cm.bone,
        cmap=plt.cm.jet_r,
        #extend='both',
        extend='max',
        origin=origin,
        )

    # Note that in the following, we explicitly pass in a subset of
    # the contour levels used for the filled contours.  Alternatively,
    # We could pass in additional levels to provide extra resolution,
    # or leave out the levels kwarg to use all of the original levels.

    CS2 = plt.contour(
        CS,
        levels=CS.levels,
        colors='white',
        linewidths=0.05,
        alpha=0.3,
        origin=origin,
        # hold='on',
        )

    cbar = plt.colorbar(CS)
    #cbar.ax.set_ylabel('Overpotential [V]')
    #cbar.ax.set_ylabel(r'$\eta_{\sf calc.}$')
    cbar.ax.set_ylabel(r'$\eta_{\sf OER}$')

    ax.tick_params(axis='both', direction='out')
    ax.get_xaxis().tick_bottom()   # remove unneeded ticks
    ax.get_yaxis().tick_left()

    offset = [0.0, 0.08]

    #foo=r': %f' % (calc_systems[i][3])
    for i in range(len(calc_systems)):

        # ax.plot(calc_systems[i][1]-calc_systems[i][0], calc_systems[i][0],
        # 'or',color=calc_systems[i][5])

        marker_color = calc_systems[i][10]
        marker_border_color = calc_systems[i][5]
        marker_border_width = 1.

        if marker_border_color is None:
            marker_border_width = 0.

        x_i = calc_systems[i][1] - calc_systems[i][0]
        y_i = calc_systems[i][0]
        size_i = calc_systems[i][9]

        lim_pot_i = calc_systems[i][3] + 1.23
        label_i = calc_systems[i][4] + ' : %.2f V' % (lim_pot_i)

        ax.plot(
            x_i,
            y_i,
            size_i,
            mec=marker_border_color,
            mew=marker_border_width,
            mfc=marker_color,
            zorder=4,
            marker=calc_systems[i][11],
            label=calc_systems[i][4] + ' : %.2f V' % (calc_systems[i][3])
            )

    #| - __old__
    # if i!=0 and 1:
    # ax.text(calc_systems[i][1]-calc_systems[i][0]+calc_systems[i][6],
    # calc_systems[i][0]+calc_systems[i][7],
    # calc_systems[i][4]+'(%.2f)' %(calc_systems[i][3]),color='black',
    # fontsize=6,horizontalalignment='center',rotation=0,zorder=1)
    #  else:
    #      ax.text(calc_systems[i][1]-calc_systems[i][0]+calc_systems[i][6],
    # calc_systems[i][0]+calc_systems[i][7],
    # calc_systems[i][4]+'(%.2f)' %(calc_systems[i][3]),
    # color='white',fontsize=6,horizontalalignment='center',
    # rotation=0,zorder=1)
    #ax.text(calc_systems[i][0],calc_systems[i][1],'%i' %(i+1),
    # color='black',fontsize=4,
    #        horizontalalignment='center',
    #        verticalalignment='center',
    #        rotation=0,zorder=2)
    #__|

    corners = [
        [1.3, 1.0],
        [x1 + (x2 - x2) * 0.2, y1 + (y2 - y1) * 0.9],
        [x1 + (x2 - x2) * 0.8, y1 + (y2 - y1) * 0.1],
        [-2, 0],
        ]

    #for i in range(len(corners)):
    #   ax.text(corners[i][0],corners[i][1], overpotential_label(corners[i][0],
    # corners[i][1]), color='white',fontsize='x-small',
    # horizontalalignment='center',rotation=0,zorder=3)

    ax.legend(
        bbox_to_anchor=(1.25, 1.05),
        loc=2,
        borderaxespad=1,
        ncol=1,
        fancybox=True,
        shadow=True,
        fontsize='x-small',
        handlelength=2,
        )

    fig_path_i = os.path.join(
        plots_folder,
        "OER_2D_Volcano.pdf",
        )

    if not os.path.exists(plots_folder):
        os.makedirs(plots_folder)

    fig.savefig(
        fig_path_i,
        bbox_inches="tight",
        )

    # fig.savefig('OER_contour_plot_v13.pdf', bbox_inches='tight')
    fig.clf()

    #__| **********************************************************************

    # *************************************************************************
    #| - OER_scaling **********************************************************
    # fig = plt.figure(figsize=fig_size, dpi=300)
    # ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    # x1 = -1
    # x2 = 2.5
    # ax.axis([x1, x2, x1, ooh_oh_scaling(x2)])
    #
    # ax.set_xlabel(r'$\Delta$G$_{\sf OH}$ (eV)')
    # ax.set_ylabel(r'$\Delta$G$_{\sf OOH}$,$\Delta$G$_{\sf O}$ (eV)')
    #
    # xdata = []
    # ydata = []
    # y2data = []
    #
    # #for i in range(3):
    # for i in range(len(calc_systems)):
    #     xdata.append(calc_systems[i][0])
    #     ydata.append(calc_systems[i][2])
    #     y2data.append(calc_systems[i][1])
    #
    # # print(xdata)
    # # print(ydata)
    #
    # fit = polyfit(xdata, ydata, 1)
    # fit_fn = poly1d(fit)
    # print(fit_fn)
    # aa = fit_fn[1]
    # bb = fit_fn[0]
    #
    # fit1 = polyfit(xdata, y2data, 1)
    # fit_fn1 = poly1d(fit1)
    # print(fit_fn1)
    #
    # #print fit_fn[0], fit_fn[1]
    # # #how bad is scaling
    # # for i in range(len(calc_systems)):
    # #     error = calc_systems[i][2] - \
    # #         (fit_fn[1] * calc_systems[i][0] + fit_fn[0])
    # #
    # #     print(error, calc_systems[i])
    #
    # xx = np.arange(x1, x2, delta)
    #
    # # Plotting Scaling Lines
    # ax.plot(xx, fit_fn[1] * xx + fit_fn[0], '--',
    #     lw=1, dashes=(3, 1), c='grey', label='OOH scaling',
    #     )
    #
    # ax.plot(xx, xx + 3.2, '--', lw=1, dashes=(3, 1), c='black')
    #
    # ax.plot(xx, xx, '--', lw=1, dashes=(3, 1), c='black')
    #
    # ax.plot(xx, fit_fn1[1] * xx + fit_fn1[0], '--',
    #     lw=1, dashes=(3, 1), c='red', label='O scaling',
    #     )
    #
    # for i in range(len(calc_systems)):
    #     ax.plot(
    #         calc_systems[i][0],
    #         calc_systems[i][2],
    #         'ro',
    #         ms=3,
    #         marker=calc_systems[i][11],
    #         #alpha=0.2,
    #         color=calc_systems[i][10],
    #         )
    #
    #     ax.plot(
    #         calc_systems[i][0],
    #         calc_systems[i][1],
    #         'ro',
    #         ms=3,
    #         marker=calc_systems[i][11],
    #         #alpha=0.2,
    #         color=calc_systems[i][10],
    #         )
    #
    #     ax.plot(
    #         calc_systems[i][0],
    #         calc_systems[i][0],
    #         calc_systems[i][9],
    #         mec=calc_systems[i][5],
    #         mfc=calc_systems[i][10],
    #         mew=0.8,
    #         zorder=4,
    #         marker=calc_systems[i][11],
    #         label=calc_systems[i][4] + ' : %.2f V' % (calc_systems[i][3]),
    #         )
    #
    #     # ax.text(calc_systems[i][0],
    #     # calc_systems[i][0]+calc_systems[i][7]+0.08,
    #     # calc_systems[i][4]+'(%.2f)' %(calc_systems[i][3]),
    #     # color='black',fontsize=6,horizontalalignment='center',
    #     # rotation=0,zorder=1)
    #
    # ax.legend(
    #     bbox_to_anchor=(1.05, 1.05),
    #     loc=2,
    #     borderaxespad=0.5,
    #     ncol=1,
    #     fancybox=True,
    #     shadow=True,
    #     fontsize='x-small',
    #     handlelength=2,
    #     )
    #
    #
    # fig_path_i = os.path.join(
    #     plots_folder,
    #     "OER_scaling.pdf",
    #     )
    #
    # if not os.path.exists(plots_folder):
    #     os.makedirs(plots_folder)
    #
    # fig.savefig(
    #     fig_path_i,
    #     bbox_inches="tight",
    #     )
    #
    # # fig.savefig('OER_scaling.pdf', bbox_inches='tight')
    #
    # fig.clf()
    #
    #__| **********************************************************************

    # *************************************************************************
    #| - OER_1D_plot **********************************************************
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])

    x1 = axes_ranges["x1"]
    x2 = axes_ranges["x2"]
    y2 = axes_ranges["y2"]
    y1 = axes_ranges["y1"]

    ax.axis([x1, x2, y1, y2])
    delta = 0.01
    x = np.arange(x1, x2, delta)

    ax.set_xlabel(r'$\Delta$G$_{\sf O}-\Delta$G$_{\sf OH}$ (eV)')

    #ax.set_ylabel(r'$\Delta$G$_{\sf O}$ (eV)')
    # ax.set_ylabel(r'U_{\sf OER}$ (V)')
    # ax.set_ylabel(r'$\eta_{\sf OER}$')

    ax.set_ylabel(r'Limiting Potential (V)')
    ax.set_ylim(ax.get_ylim()[::-1])

    #| - Plotting Volcano Lines
    plot(
        x,
        np.maximum(x, 3.2 - x),
        '--',
        color='black',
        lw=0.67,
        dashes=(3, 1),
        zorder=2,
        )
    #__|

    #| - Plotting Data Points
    for i in range(len(calc_systems)):

        marker_color = calc_systems[i][10]
        marker_border_color = calc_systems[i][5]
        marker_border_width = 1.

        if marker_border_color is None:
            marker_border_width = 0.

        x_i = calc_systems[i][1] - calc_systems[i][0]
        y_i = calc_systems[i][3] + 1.23
        size_i = calc_systems[i][9]

        lim_pot_i = calc_systems[i][3] + 1.23
        label_i = calc_systems[i][4] + ' : %.2f V' % (lim_pot_i)
        ax.plot(

            x_i,
            y_i,
            size_i,

            mec=marker_border_color,
            mew=marker_border_width,
            mfc=marker_color,

            zorder=4,
            marker=calc_systems[i][11],
            label=label_i,
            )

        #| - NEW | Adding Labels
        if annotate_data:
            x_tmp = calc_systems[i][1] - calc_systems[i][0]
            y_tmp = calc_systems[i][3] + 1.23

            sys_i_name = calc_systems[i][4]

            xytext_i = annotation_placement_dict.get(sys_i_name, None)

            if xytext_i is None:
                xytext_i = (4, 10)
                # xytext_i = (-10, 10)
            else:
                xytext_i = xytext_i[0]

            facet_i = calc_systems[i][-1]["facet"]

            ax.annotate(
                str(facet_i),

                # xy=(x_tmp, y_tmp), xytext=(-10, 10),
                xy=(x_tmp, y_tmp), xytext=xytext_i,

                textcoords='offset points', ha='right', va='bottom',

                # bbox=dict(
                #     boxstyle='round,pad=0.5',
                #     fc='yellow',
                #     alpha=0.5,
                #     ),

                arrowprops=dict(
                    arrowstyle='->',
                    connectionstyle='arc3,rad=0',
                    ),
                size=font_size - 2,
                )

        #__|

    #__|

    #| - Plotting Experimental Overpotentials <--------------------------------

    # 1.45 for IrO3
    # and 1.6 for IrO2
    # and 1.57 for IrOx (some unknown polymorph from Ir-metal they made)

    if exp_lines:

        #| - IrO3
        pos_IrO3 = 1.45
        color_IrO3 = "red"
        plot(
            [x1, x2],
            2 * [pos_IrO3, ],
            '--',
            color=color_IrO3,
            lw=0.67,
            dashes=(1, 1),
            zorder=2,
            )

        plt.text(
            x1,
            pos_IrO3 - 0.014,
            r'$IrO_3$',
            fontsize=font_size,
            color=color_IrO3,
            horizontalalignment='left',
            verticalalignment='center',
            )

        #__|

        #| - IrO2
        pos_IrO2 = 1.8
        color_IrO3 = "blue"
        plot(
            [x1, x2],
            2 * [pos_IrO2, ],
            '--',
            color=color_IrO3,
            lw=0.67,
            dashes=(1, 1),
            zorder=2,
            )

        plt.text(
            x1,
            pos_IrO2 + 0.024,
            r'$IrO_2$',
            fontsize=font_size,
            color=color_IrO3,
            horizontalalignment='left',
            verticalalignment='center',
            )

        #__|

        #| - IrOx
        pos_IrOx = 1.57
        color_IrOx = "green"
        plot(
            [x1, x2],
            2 * [pos_IrOx, ],
            '--',
            color='green',
            lw=0.67,
            dashes=(1, 1),
            zorder=2,
            )

        plt.text(
            x1,
            pos_IrOx - 0.014,
            r'$IrO_x$',
            fontsize=font_size,
            color=color_IrOx,
            horizontalalignment='left',
            verticalalignment='center',
            # transform=ax.transAxes,
            )
        #__|

    #__|

    ax.legend(
        bbox_to_anchor=(-0.15, 1.425),
        loc=2,
        borderaxespad=0.5,
        ncol=3,
        fancybox=True,
        shadow=False,
        fontsize="x-small",
        handlelength=2,
        )

    # fig.savefig('OER_1D_plot_v13.pdf', bbox_inches='tight')

    fig_path_i = os.path.join(
        plots_folder,
        "OER_1D_Volcano.pdf",
        )

    # fig_path_i_svg = os.path.join(
    #     plots_folder,
    #     "OER_1D_plot_v13.svg",
    #     )

    if not os.path.exists(plots_folder):
        os.makedirs(plots_folder)

    fig.savefig(
        fig_path_i,
        bbox_inches="tight",
        )

    # fig.savefig(
    #     fig_path_i_svg,
    #     bbox_inches="tight",
    #     )

    fig.clf()

    #__| **********************************************************************

    #__|


#| - __old__

    #| - __old__
    #plot(x,1.23,'--',color='black',lw=0.67, dashes=(3,1),zorder=2)
    # xy=np.array([xp for xp in x if 1.55<xp<1.66])
    # ax.fill_between(xy, y2, np.maximum(xy,3.2-xy)-1.23,
    # zorder=1, color='red', alpha=0.3, edgecolor="none")
    # for b in x:
    # if(np.maximum(b, 3.2-b) < 0.44):
    # print(b)

    #plot(x,np.maximum(x,bb-x*(aa)-0.65)-1.23,'--',
    # color='grey',lw=0.67, dashes=(3,1))
    #slope not one
    #plot(x,np.maximum(x,3.18-0.82*x-0.35)-1.23,'--',color='pink',
    # lw=0.67,dashes=(3,1))
    #plot(x,np.maximum(x,2.46-x)-1.23,'-',color='black',lw=0.67)

    #import matplotlib.patches as patches
    #ax.add_patch(
    #    patches.Rectangle(
    #        (calc_systems[1][1]-calc_systems[1][0],
    # calc_systems[1][3]-0.04),   # (x,y)
    #        0.25,          # width
    #        0.05,          # height
    #        fill=True,
    #        edgecolor="none",
    #        facecolor='red',
    #    )
    #)
    #__|

    #| - __old__
    # if(i!=1):
    # ax.text(calc_systems[i][1]-calc_systems[i][0],
    # calc_systems[i][3]-0.02,calc_systems[i][3])
    # color='black',fontsize=6,horizontalalignment='left',rotation=0,zorder=4)
    # else:
    # ax.text(calc_systems[i][1]-calc_systems[i][0],calc_systems[i][3]-0.02,
    # calc_systems[i][4]+'(%.2f)' %(calc_systems[i][3]),
    # color='black',fontsize=6,horizontalalignment='right',rotation=0,zorder=4)
    #__|

    #| - __old__
    #levels = np.arange(0, 2, 0.05)
    #CS = plt.contourf(X,Y,Z, levels, cmap=cm.jet_r, origin='lower')
    #CS = plt.contourf(X,Y,Z, levels, origin='lower')
    #im = plt.imshow(Z, interpolation='bilinear', origin='lower',
    #                cmap=cm.jet, extent=(x1,x2,y1,y2))
    #levels2 = [2.0]
    #CS2 = plt.contour(CS, levels2,
    #                        colors = 'r',
    #                        origin='lower',
    #                        hold='on')
    #CS = plt.contour(Z, levels,
    #                 origin='lower',
    #                 linewidths=0.5,
    #                 extent=(x1,x2,y1,y2))
    ##Thicken the zero contour.
    #zc = CS.collections[6]
    #plt.setp(zc, linewidth=2)
    #__|

    #| - __old__
    #cbar.add_lines(CS2)
    #plt.clabel(CS, levels[1::2],  # label every second level
    #           inline=1,
    #           fmt='%1.1f',
    #           fontsize='x-small')
    #plt.title('Lines with colorbar')
    # We can still add a colorbar for the image, too.
    # This makes the original colorbar look a bit out of place,
    # so let's improve its position.
    #__|

    #| - __old__
    #plot(x,ooh_oh_scaling(x),'--',color='orange',lw=1,
    # dashes=(3,1),label='$\Delta$G$_{\sf OOH}$=0.82G$_{\sf OH}$+3.18 eV')
    #ax.text(x1+0.02,y2-0.3,
    # '$\Delta$G$_{\sf OOH}$=%.2fG$_{\sf OH}$+%.2f eV' %(fit[0],fit[1]),
    # color='orange',fontsize='x-small',zorder=10,horizontalalignment='left')
    #ax.text(x1+0.02,y2-0.3,
    # '$\Delta$G$_{\sf OOH}$=%.2fG$_{\sf OH}$+%.2f eV' %(0.82,3.18),
    # color='orange',fontsize='x-small',zorder=10,horizontalalignment='left')
    #plt.show()
    #__|

#__|
