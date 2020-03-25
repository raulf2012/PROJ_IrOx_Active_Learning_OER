#!/usr/bin/env python

"""Plot electrochemical Pourbaix diagram.

Author: Michal Bajdich
"""


def plot_Pourbaix(

    ):
    """
    """
    # | - plot_Pourbaix

    # | - Import Modules
    from math import pow
    from pylab import *
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    from matplotlib.ticker import *
    from scipy import *
    import subprocess

    #Import plot_settings as ps
    from matplotlib import rc
    #__|

    # | - Setup
    # with golden ration and the whole shebang ...
    # settings size and font for revtex stylesheet
    fig_width_pt = 1.8 * 246.0  # Get this from LaTeX using \showthe\columnwidth
    #fig_width_pt *= 300./72 # convert to 300 dpi
    inches_per_pt = 1.0 / 72.27               # Convert pt to inches
    #inches_per_pt = 1.0/300               # Convert pt to inches
    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height =fig_width*golden_mean       # height in inches
    fig_size = [fig_width,fig_height]

    font_size = 10
    tick_font_size = 10
    xlabel_pad = 8
    ylabel_pad = 18

    matplotlib.rcParams['ps.usedistiller'] = 'xpdf'
    matplotlib.rcParams['font.family'] = 'sans-serif'
    #matplotlib.rcParams['font.family'] = 'serif'
    #matplotlib.rcParams['font.serif'] = 'Arial'
    matplotlib.rcParams['font.sans-serif'] = 'Arial'
    matplotlib.rcParams['font.size'] = 10
    matplotlib.rcParams['axes.labelsize'] = 2*font_size
    matplotlib.rcParams['legend.fontsize'] = font_size
    matplotlib.rcParams['xtick.labelsize'] = tick_font_size
    matplotlib.rcParams['ytick.labelsize'] = tick_font_size
    matplotlib.rcParams['mathtext.default'] = 'regular'


    matplotlib.rcParams['lines.linewidth'] = 1.
    fig = plt.figure(figsize=fig_size,dpi=300)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])

    Umin=0.0
    Umax=2.2



    #Values for X,Y
    pH = arange(0,14,0.10)
    U = arange(Umin,Umax,0.01)
    Umax2=Umax+0.06*14
    U2= arange(Umin,Umax2,0.01)

    ax.axis([0, 14, Umin,Umax])
    ax.set_xlabel(r'pH')
    ax.set_ylabel(r'U/V')

    #plt.xlabel('$pH$',family = 'Helvetica', size='40')
    #plt.ylabel('$U/V$',family = 'Helvetica',size='40')

    #Constants
    kbt = 0.0256
    const = kbt*log(10)

    extraticks=[1.23]
    plt.yticks(list(plt.yticks()[0]) + extraticks)
    #__|

    # | - Methods
    def frhe(x):
        return 1.78-x*const
    def fh2oequi(x):
        return 1.23-x*const

    def addO(x,y):
     return -(h2o-h2)-2*(y+x*const)+dso
    def addOH(x,y):
     return -(h2o-0.5*h2)-(y+x*const)+dsoh
    def addOOH(x,y):
     return -(2*h2o-1.5*h2)-3*(y+x*const)+dsooh

    def addH2O(x,y):
     return -(h2o)-(zpeh2o-tsh2o)

    def addH(x,y):
     return -0.5*h2+1*(y+x*const)+dsh

    #Function to calculate DG
    def dg(i,x,y):
        if(surfs[i][1]+surfs[i][2]+surfs[i][3]>0 and x==0 and y==0):
        	print (surfs[i][0] -surfs[0][0] +surfs[i][1]*addH(x,y) +surfs[i][2]*addO(x,y) +surfs[i][3]*addOH(x,y))/(surfs[i][1]+surfs[i][2]+surfs[i][3])
        return surfs[i][0] -surfs[0][0] +surfs[i][1]*addH(x,y) +surfs[i][2]*addO(x,y) +surfs[i][3]*addOH(x,y)

    #__|

    #--------------------------------------------
    # Total  Energies of the different coverages |
    #--------------------------------------------

    # raw_e, H, O, OH
    surfs = [
    [-178.96067932, 0, 0, 0], #clean
    [-199.44772242, 0, 4, 0], #4O
    [-219.47269277, 0, 0, 4], #4OH*
    ]

    #------------------------
    # Calculation of the DG |
    #------------------------

    #400 eV ,O_s
    #h2=  -6.759300
    #h2o= -14.01977

    #500 eV, O regular
    h2=-6.77014123
    h2o=-14.21744725

    #Max claculated
    #H2	6.770141	-	0.26810	0.09048	0.00136	0.408000	6.720721
    #H2O	-14.217447	-	0.56733	0.00010	0.00186	0.558000	-14.208017

    zpeh2o=0.5741 #exp. NIST 0.5584250
    zpeh2=0.27283 #exp. NIST 0.27283
    tsh2o=0.67 # 1 bar pressure
    tsh2=0.41 #exp.

    #monicas new ZPEs
    zpeo=0.065 #0.113 Max. G IrO2
    zpeoh=0.37 #0.397 Max. G IrO2
    zpeooh=0.44 #stays the same

    dso=zpeo -(zpeh2o-zpeh2 -tsh2o+tsh2)
    dsoh=zpeoh -(zpeh2o -0.5*zpeh2 -tsh2o+ 0.5*tsh2)
    dsooh=zpeooh-(2*zpeh2o -1.5*zpeh2 -2*tsh2o+ 1.5*tsh2)
    dsh = dsoh-dso
    print dsoh, dso, dsh

    #find intersects

    nsurfs=len(surfs)

    i=0  #Ph=0
    lowest_surfaces=[]
    for j in U2:
        #print i,j
        values=[]

        for k in range(nsurfs):
                #print k,dg(k,i,j)
                values.append(dg(k,i,j))
        sorted_values=sorted(range(len(values)), key=lambda k: values[k])
        lowest_surfaces.append(sorted_values[0])
        #print j, values

    #print lowest_surfaces

    crossover=[]
    uniquesurf=[]
    uniquesurf.append(lowest_surfaces[0])
    old_value=lowest_surfaces[0]
    crossover.append(Umin)
    for j in range(len(U2)):
        if(lowest_surfaces[j]!=old_value):
            uniquesurf.append(lowest_surfaces[j])
            crossover.append(U2[j])
            old_value=lowest_surfaces[j]

    crossover.append(Umax2)

    print crossover
    print uniquesurf

    color=['turquoise', 'green', 'red','blue', 'gray', 'gold', 'gray20']
    pH2 = arange(0,14,0.01)


    for i in range(len(uniquesurf)):
        k=uniquesurf[i]
        foo= r"S$_{%i}$(H-%i O-%i OH-%i)" % (k, surfs[k][1], surfs[k][2],surfs[k][3])
        #fbk = {'lw':0.0, 'edgecolor':color[i]}
        fbk = {'lw':0.5, 'edgecolor':'black'}
        fill_between(pH2,crossover[i]-pH2*const, crossover[i+1]-pH2*const, facecolor=color[i],alpha=0.3,**fbk)
        plot([], [], color=color[i],alpha=0.3, linewidth=5,label=foo)

    Vover=0.526
    y=1.23+Vover -pH2*const
    llabel='$\eta$ = '+ repr(Vover)+' V at S$_{12}$'
    plot(pH2,y,'-',color='black', lw=1, dashes=(3,1),label=llabel)

    plot(pH2,1.23-pH2*const,'--',color='blue',lw=1, dashes=(3,1))
    ax.text(0.2,1.25,r'2H$_2$O $\leftrightarrow$ 4H$^+$ +O$_2$+4e$^-$',color='blue',rotation=-11.6,fontsize='x-small')


    #legend()
    legend(ncol=2, fancybox=True, shadow=True, fontsize='x-small',handlelength=3)
    fig.savefig('Pourbaix_simple_generic.pdf', bbox_inches='tight')
    #show()
    #__|
