# #!/usr/bin/env python

from scipy import *
from ase.io import *


#------------------------
# Calculation of the DG |
#------------------------

##400 eV ,O_s
#h2=  -6.759300
#h2o= -14.01977

#500 eV, O regular
h2=  -6.77149190
h2o= -14.23091949

##600 eV, O regular
#h2= -6.77412273
#h2o= -14.23797429

#Max claculated
#H2	6.770141	-	0.26810	0.09048	0.00136	0.408000	6.720721
#H2O	-14.217447	-	0.56733	0.00010	0.00186	0.558000	-14.208017

#Contributions to Gibbs Energies for gas molecules (VASP-PBE calculated by Max; T= 300K)
zpeh2o=0.560	#exp. NIST 0.558
zpeh2=0.268 	#exp. NIST 0.273
cvh2o=0.103  	#From Colin at P = 0.035 bar
cvh2=0.0905
tsh2o=0.675 	#From Colin at P = 0.035 bar
tsh2=0.408 	#at P = 1bar

#Contributions to Gibbs Energies for adsorbates (VASP-PBE calculated by Max using Michal data for NiCe; T= 300K)
#if(PBE_VASP):
zpeo=0.065 #0.061
zpeoh=0.344 #0.360
zpeooh=0.443 #0.468 #0.459 old
cvo=0.038 #0.0325
cvoh=0.051 #0.049
cvooh=0.068 #0.077
tso=0.080 #0.051	#0.060 From Colin
tsoh=0.080 #0.085 From Colin
tsooh=0.116  #0.135	#0.215 From Colin

#Gibbs Energies for the gas molecules
dgh2o=zpeh2o +cvh2o -tsh2o
dgh2=zpeh2 +cvh2 -tsh2

#Gibbs Energy corrections for adsorbates
dso=zpeo +cvo -tso -(dgh2o -dgh2)
dsoh=zpeoh +cvoh -tsoh -(dgh2o -0.5*dgh2)
dsooh=zpeooh +cvooh -tsooh -(2*dgh2o -1.5*dgh2)
dsh=dsoh-dso

#dso=zpeo -(zpeh2o-zpeh2 -tsh2o+tsh2)
#dsoh=zpeoh -(zpeh2o -0.5*zpeh2 -tsh2o+ 0.5*tsh2)
#dsooh=zpeooh-(2*zpeh2o -1.5*zpeh2 -2*tsh2o+ 1.5*tsh2)
#dsh = dsoh-dso
print(dsoh, dso, dsh)

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

print(addO(0,0), addOH(0,0), addH(0,0))

# #Function to calculate DG
# def dg(i,x,y):
#     return surfs[i][0] -surfs[0][0] +surfs[i][1]*addH(x,y) +surfs[i][2]*addO(x,y) +surfs[i][3]*addOH(x,y)
