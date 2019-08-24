#------------------------
# Calculation of the DG |
#------------------------

# 400 eV ,O_s
#h2=  -6.759300
#h2o= -14.01977
#h2s= -11.203734 #pallavi's number

# 400 eV ,O regular like CoS2 Pallavi
# h2=  -6.759300
# h2o= -14.242
# h2s= -11.203734 #pallavi's number

# 500 eV, O regular
h2=-6.77014123
h2o=-14.21744725

# Max claculated
# H2     6.770141        -       0.26810 0.09048 0.00136 0.408000        6.720721
# H2O    -14.217447      -       0.56733 0.00010 0.00186 0.558000        -14.208017

#Contributions to Gibbs Energies for gas molecules (VASP-PBE calculated by Max; T= 300K)
zpeh2o=0.560    #exp. NIST 0.558
zpeh2=0.268     #exp. NIST 0.273
zpeh2s=0.3982     #exp. NIST 0.273
cvh2o=0.103     #From Colin at P = 0.035 bar
cvh2=0.0905
cvh2s=0.1032     # from 298.15*205.81/kjmol/1000
tsh2o=0.675     #From Colin at P = 0.035 bar
tsh2s=0.635     # from 298.15*205.753/kjmol/1000 OK
tsh2=0.408      #at P = 1bar


#Contributions to Gibbs Energies for adsorbates (VASP-PBE calculated by Max using Michal data for NiCe; T= 300K)
#if(PBE_VASP):
zpeo=0.065 #0.061
zpeoh=0.344 #0.360
zpeooh=0.443 #0.468 #0.459 old
cvo=0.038 #0.0325
cvoh=0.051 #0.049
cvooh=0.068 #0.077
tso=0.080 #0.051        #0.060 From Colin
tsoh=0.080 #0.085 From Colin
tsooh=0.116  #0.135     #0.215 From Colin

#else:
#    if (RPBE_GPAW):
#        zpeo=0.065
#        zpeoh=0.37
#        zpeooh=0.44
#    else:
#        print 'Dont know ZPEs of adsorbates'
#        quit()


#Gibbs Energies for the gas molecules
dgh2o=zpeh2o +cvh2o -tsh2o
dgh2=zpeh2 +cvh2 -tsh2
dgh2s=zpeh2s +cvh2s -tsh2s


#IrOx_free_energy_calculator.py
#!/usr/bin/python
kjmol=96.485
T=298.15
convert_dS=T/(1000*kjmol)



h2ozpe=0.5584   #zpeh2o=0.560  #exp. NIST 0.558
h2ocp=0.10260 #NIST         0.10      #cvh2o=0.103

h2zpe=0.27283   #zpeh2=0.268  #exp. NIST 0.273
h2cp=0.087785    #=NIST  0.09       #cvh2=0.0905

o2zpe=0.098
o2cp=0.089962

h_ref=(h2+h2zpe+h2cp)/2 #dH at 300 K
h2o_ref=h2o+h2ozpe+h2ocp #dH at 300 K


#o_ref2=h2o_ref-2*h_ref +2.506 -(o2zpe+o2cp)/2  #dE relative to -241.81/kjmol
o_ref2=h2o_ref-2*h_ref +2.506  #dH fitting relative to -241.81/kjmol

#print '1/2 O2 vibrations fo dH ', (o2zpe+o2cp)/2

#o_ref2_de= o_ref2 -0.5*(o2zpe+o2cp) #dE only
#o_ref2=(h2o-h2 +2.506)+(h2ozpe-h2zpe-0.5*o2zpe)+(h2ocp-h2cp-0.5*o2cp) #dH

print 'Error for dH of H2O, 1/2 O2 ref fitted ',h2o_ref-2*h_ref-o_ref2 -(-241.81/kjmol), ' should be close to zero!'

S_o2_gas=205.2 #exp.
S_h2_gas=130.68 #exp.

tsh2o=0.675     #From Colin at P = 0.035 bar
tsh2=S_h2_gas*convert_dS        #at P = 1bar
tso2=S_o2_gas*convert_dS


#check H2O termochemistry
print 'dH molecular refences relative to H2O and H2 at 300K'
print '1/2 O2 ref fitted to dH of H2O ',o_ref2,' 1/2 H2 ref. ', h_ref

#print h2o_ref -2*h_ref -o_ref2
#print tsh2o-tsh2-0.5*tso2
#print h2o_ref -2*h_ref -o_ref2 -(tsh2o-tsh2-0.5*tso2)

#print  (h2o_ref +tsh2o - (2*h_ref+tsh2) -(o_ref2-0.5*(o2zpe+o2cp)+0.5*tso2))
print 'Error for dG of H2O, 1/2 O2 ref fitted: ', (h2o_ref -2*h_ref -o_ref2 -(tsh2o-tsh2-0.5*tso2)) -(-237.140/kjmol), ' should be close to zero!'
#print 'Error for dG of H2O: 1/2 O2 ref PBE   : ', (h2o_ref -2*h_ref -o_ref1 -(tsh2o-tsh2-0.5*tso2)) -(-237.140/kjmol)
