#!/usr/bin/python

#| - Import Modules
import sys
import os

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "data"))

from proj_data_irox import (
    zpe_h2o,
    cv_h2o,
    ts_h2o,

    zpe_h2,
    cv_h2,
    ts_h2,

    ads_fe_dict,
    )


#| - old imports
# from math import pow
# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.path import Path
# from matplotlib.patches import PathPatch
# from pylab import *
# import subprocess

# from matplotlib.ticker import *
from scipy import *
#__|

#__|

#| - Reassigning my corrections to the correct variable names

zpeh2o = zpe_h2o
cvh2o = cv_h2o
tsh2o = ts_h2o
zpeh2 = zpe_h2
cvh2 = cv_h2
tsh2 = ts_h2

zpeoh = ads_fe_dict["oh"]["zpe"]
cvoh = ads_fe_dict["oh"]["cv"]
tsoh = ads_fe_dict["oh"]["ts"]

zpeo = ads_fe_dict["o"]["zpe"]
cvo = ads_fe_dict["o"]["cv"]
tso = ads_fe_dict["o"]["ts"]

zpeooh = ads_fe_dict["ooh"]["zpe"]
cvooh = ads_fe_dict["ooh"]["cv"]
tsooh = ads_fe_dict["ooh"]["ts"]

#__|

def overpot(deoh, deo, deooh):
    """
    """
    #| - overpot

    #| - Constants
    #Constants
    kbt = 0.0256
    const = kbt*log(10)


    #Contributions to Gibbs Energies for gas molecules (VASP-PBE calculated by Max; T= 300K)
    zpeh2o=0.560    #exp. NIST 0.558
    zpeh2=0.268     #exp. NIST 0.273
    cvh2o=0.103      #From Colin at P = 0.035 bar
    cvh2=0.0905
    tsh2o=0.675     #From Colin at P = 0.035 bar
    tsh2=0.408     #at P = 1bar


    #Contributions to Gibbs Energies for adsorbates (VASP-PBE calculated by Max using Michal data for NiCe; T= 300K)
    #if(PBE_VASP):
    zpeo=0.065 #0.061
    zpeoh=0.344 #0.360
    zpeooh=0.443 #0.468 #0.459 old
    cvo=0.038 #0.0325
    cvoh=0.051 #0.049
    cvooh=0.068 #0.077
    tso=0.080 #0.051    #0.060 From Colin
    tsoh=0.080 #0.085 From Colin
    tsooh=0.116  #0.135    #0.215 From Colin
    #__|

    #| - __old__
    #else:
    #    if (RPBE_GPAW):
    #        zpeo=0.065
    #        zpeoh=0.37
    #        zpeooh=0.44
    #    else:
    #        print 'Dont know ZPEs of adsorbates'
    #        quit()
    #__|

    #Gibbs Energies for the gas molecules
    dgh2o=zpeh2o +cvh2o -tsh2o
    dgh2=zpeh2 +cvh2 -tsh2

    #Gibbs Energy corrections for adsorbates
    dgo=zpeo +cvo -tso -(dgh2o -dgh2)
    dgoh=zpeoh +cvoh -tsoh -(dgh2o -0.5*dgh2)
    dgooh=zpeooh +cvooh -tsooh -(2*dgh2o -1.5*dgh2)
    dgh=dgoh-dgo

    print(dgo, dgoh, dgooh)


    #We calculate the OER steps and the overpotential
    dgde=[dgoh, dgo-dgoh, dgooh-dgo, -dgooh]

    # print("##################################")
    print("    dG corrections to Es in standard OER mechanism")
    print("    dG1-dE1   dG2-dE2   dG3-dE3   dG4-dE4")
    print("    %.3f     %.3f    %.3f     %.3f" %(dgde[0], dgde[1], dgde[2], dgde[3]))
    print("##################################")


    print("---------------------------------------------")
    print('    DEoh      DEo       DEooh     DEooh-DEoh')
    print('    %.3f     %.3f     %.3f     %.3f     ' %(deoh, deo, deooh, deooh-deoh))
    print('    DGoh      DGo       DGooh     DGooh-DGoh')
    print('    %.3f     %.3f     %.3f     %.3f     ' %(deoh+dgoh, deo+dgo, deooh+dgooh, deooh-deoh+(dgooh-dgoh)))

    print('    DG1       DG2       DG3       DG4:')
    dg=[deoh+dgde[0],deo-deoh+dgde[1],deooh-deo+dgde[2],-deooh+4.92+dgde[3]]
    print('    %.3f     %.3f     %.3f     %.3f'%(dg[0], dg[1], dg[2], dg[3]))

    step=0
    max1=max(dg)
    for i in range(len(dg)):
        if(dg[i]==max1):
                step=i+1
    print("---------------------------------------------")
    print("    The Overpotential (step %i): %.3f V " %(step,max1-1.23))
    print("---------------------------------------------")

    return(max1-1.23)
    #__|

#| - __main__

if __name__ == '__main__':

    #| - __old__
    #VASP-PBE Energies for adsorbates
    #PBE_VASP=1
    #RPBE_GPAW=0
    #print "---------------------------------------------"
    #print 'Number of arguments:', len(sys.argv), 'arguments.'
    #__|

    print('Argument List:', str(sys.argv))

    if(len(sys.argv) < 2):
        print('Error: need either raw (-raw -300.05512795 -3    2.845     3.913 ) energies')
        quit()

    if(sys.argv[1] == '-rawsoftsymmetric'):
        #| - rawsoftsymmetric
        clean=float(sys.argv[2])
        eoh=float(sys.argv[3])
        eo=float(sys.argv[4])
        eooh=float(sys.argv[5])

        #400 eV ,O_s
        h2=  -6.759300
        h2o= -14.019771

        deoh=(eoh-clean-2*(h2o-0.5*h2))/2.0
        deo=(eo-clean-2*(h2o-h2))/2.0
        deooh=(eooh-clean-2*(2*h2o-1.5*h2))/2.0
        #__|
    else:

        if(sys.argv[1] == '-rawsymmetric'):
            #| - rawsymmetric
            clean=float(sys.argv[2])
            eoh=float(sys.argv[3])
            eo=float(sys.argv[4])
            eooh=float(sys.argv[5])

            #500 eV, O regular
            h2= -6.77014123
            h2o= -14.21744725


            deoh=(eoh-clean-2*(h2o-0.5*h2))/2.0
            deo=(eo-clean-2*(h2o-h2))/2.0
            deooh=(eooh-clean-2*(2*h2o-1.5*h2))/2.0
            #__|
        else:

            if(sys.argv[1] == '-raw'):
                #| - raw
                clean=float(sys.argv[2])
                eoh=float(sys.argv[3])
                eo=float(sys.argv[4])
                eooh=float(sys.argv[5])

                #500 eV, O regular
                #h2= -6.77014123
                #h2o= -14.21744725
                h2=  -6.77149190
                h2o= -14.23091949
                ###400 eV, O regular (Pallavi)
                #h2= -6.7606
                #h2o=-14.23503

                deoh=eoh-clean-(h2o-0.5*h2)
                deo=eo-clean-(h2o-h2)
                deooh=eooh-clean-(2*h2o-1.5*h2)
                #__|
            else:

                if(sys.argv[1] == '-rawsoft'):
                    #| - rawsoft
                    clean=float(sys.argv[2])
                    eoh=float(sys.argv[3])
                    eo=float(sys.argv[4])
                    eooh=float(sys.argv[5])

                    #400 eV ,O_s
                    h2=  -6.759300
                    h2o= -14.019771

                    deoh=eoh-clean-(h2o-0.5*h2)
                    deo=eo-clean-(h2o-h2)
                    deooh=eooh-clean-(2*h2o-1.5*h2)
                    #__|
                else:

                    if (sys.argv[1] == '-rawsofthse'):
                        #| - rawsofthse
                        clean = float(sys.argv[2])
                        eoh = float(sys.argv[3])
                        eo = float(sys.argv[4])
                        eooh = float(sys.argv[5])

                        #HSE 400 eV ,O_s
                        h2 = -7.82007056
                        h2o = -17.24469897
                        print(h2, h2o)

                        deoh = eoh - clean - (h2o -0.5 * h2)
                        deo = eo - clean - (h2o - h2)
                        deooh = eooh - clean - (2 * h2o - 1.5 * h2)
                        #__|
                    else:

                        if (sys.argv[1]=='-rawsoftrpbe'):
                            #| - rawsoftrpbe
                            clean=float(sys.argv[2])
                            eoh=float(sys.argv[3])
                            eo=float(sys.argv[4])
                            eooh=float(sys.argv[5])

                            #RPBE 400 eV ,O_s
                            h2= -6.97813614
                            h2o= -13.95184458
                            print(h2, h2o)

                            deoh=eoh-clean-(h2o-0.5*h2)
                            deo=eo-clean-(h2o-h2)
                            deooh=eooh-clean-(2*h2o-1.5*h2)
                            #__|
                        else:
                            deoh = float(sys.argv[1])
                            deo = float(sys.argv[2])
                            deooh = float(sys.argv[3])

    print("---------------------------------------------")

    overpot(deoh, deo, deooh)

#__|
