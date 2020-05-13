#!/usr/bin/python

kjmol=96.485
T=298.15
convert_dS=T/(1000*kjmol)


#calculate these in 600 eV
##500 eV, Normal
#h2o=-14.23091949 #PBE, 500 eV
#h2=-6.77149190 #PBE, 500 eV
#o2=-9.88557216 #PBE, 500 eV
#600 eV, Normal
h2o=-14.23835553 #PBE, 600 eV
h2=-6.77402325 #PBE, 600 eV
o2=-9.89149843 #PBE, 600 eV


h2ozpe=0.5584   #zpeh2o=0.560  #exp. NIST 0.558
h2ocp=0.10260 #NIST         0.10      #cvh2o=0.103

h2zpe=0.27283   #zpeh2=0.268  #exp. NIST 0.273
h2cp=0.087785    #=NIST  0.09       #cvh2=0.0905

o2zpe=0.098
o2cp=0.089962

h_ref=(h2+h2zpe+h2cp)/2 #dH at 300 K
h2o_ref=h2o+h2ozpe+h2ocp #dH at 300 K

o_ref1_de=(o2)/2 #dE
o_ref1=(o2+o2zpe+o2cp)/2 #dH
o_ref2_de=h2o_ref-2*h_ref +2.506 -(o2zpe+o2cp)/2  #dE relative to -241.81/kjmol
o_ref2=h2o_ref-2*h_ref +2.506  #dH fitting relative to -241.81/kjmol

#print '1/2 O2 vibrations fo dH ', (o2zpe+o2cp)/2

#o_ref2_de= o_ref2 -0.5*(o2zpe+o2cp) #dE only
#o_ref2=(h2o-h2 +2.506)+(h2ozpe-h2zpe-0.5*o2zpe)+(h2ocp-h2cp-0.5*o2cp) #dH

print('Error for dH of H2O, 1/2 O2 ref fitted ',h2o_ref-2*h_ref-o_ref2 -(-241.81/kjmol), ' should be close to zero!')

S_o2_gas=205.2 #exp.
S_h2_gas=130.68 #exp.

tsh2o=0.675 	#From Colin at P = 0.035 bar
tsh2=S_h2_gas*convert_dS 	#at P = 1bar
tso2=S_o2_gas*convert_dS


#check H2O termochemistry
print('dH molecular refences relative to H2O and H2 at 300K')
print('1/2 O2 ref fitted to dH of H2O ',o_ref2,' PBE 1/2 O2 ref: ', o_ref1, ' 1/2 H2 ref. ', h_ref)

#print h2o_ref -2*h_ref -o_ref2
#print tsh2o-tsh2-0.5*tso2
#print h2o_ref -2*h_ref -o_ref2 -(tsh2o-tsh2-0.5*tso2)

#print  (h2o_ref +tsh2o - (2*h_ref+tsh2) -(o_ref2-0.5*(o2zpe+o2cp)+0.5*tso2))
print('Error for dG of H2O, 1/2 O2 ref fitted: ', (h2o_ref -2*h_ref -o_ref2 -(tsh2o-tsh2-0.5*tso2)) -(-237.140/kjmol), ' should be close to zero!')
print('Error for dG of H2O: 1/2 O2 ref PBE   : ', (h2o_ref -2*h_ref -o_ref1 -(tsh2o-tsh2-0.5*tso2)) -(-237.140/kjmol))

gmu_h2o=h2o_ref-tsh2o
gmu_h2=2*h_ref-tsh2
gmu_o_ref=  o_ref2 -(o2zpe+o2cp)/2 + 0.5*(o2zpe+o2cp)-0.5*tso2

gmu_h2o_oer=gmu_h2o+2*1.229

print('dG muH2O corrected: ', gmu_h2o)
print('dG muH2 corrected: ', gmu_h2)
print('dG muO corrected: ', gmu_o_ref)
print('mu_O under OER conditon in equlibrium with H2O and H2: ',gmu_o_ref, 'equlibrium test ', gmu_h2o_oer-gmu_h2-gmu_o_ref)
print('mu_O under OER conditon in equlibrium with H2O and H2: second version  ', gmu_h2o_oer-gmu_h2)



print('#############now our IrOxHy system########################################')
#now Ir system constants
dh_iro2_exp=-242.672/kjmol
S_ir_metal=35.5 #exp.
S_iro2_solid=58.57 #exp.
iro2=-7.045317*3 #calculated Rutile PBE energy, 600 eV
# iro3=-6.442159*4 #calculated lowest PBE energy , 600 eV
iro3=-6.469847*4 # MINE

irho3=-6.123924*5 #calculated lowest PBE energy , 600 eV

ir_metal=-8.860644725




#molecular references are dH, so they contain vibrations
# here, o_ref2 is a dH of 1/2 O2
ir_metal_fit=iro2-(2*o_ref2)-dh_iro2_exp #fit to exp dH

print('Calculated Ir metal ref : ', ir_metal, ' vs. fitted to IrO2 dH : ', ir_metal_fit)

print('Calculate dHs now')
dh_iro2=iro2-(2*o_ref2+ir_metal_fit)
print('IrO2: dH= ' ,dh_iro2, ' eV ', dh_iro2/3, dh_iro2*kjmol, ' kjmol exp: ', dh_iro2_exp*kjmol, ' kjmol')

dh_iro3=iro3-(3*o_ref2+ir_metal_fit)

print('IrO3: dH= ',dh_iro3, ' eV ',dh_iro3*kjmol, ' kjmol')

dh_irho3=irho3-(3*o_ref2+ir_metal_fit+h_ref)

print('IrHO3: dH= ',dh_irho3, ' eV ',dh_irho3*kjmol, ' kjmol')

print('Calculate dGs now')

TdS_iro2=(S_iro2_solid-S_ir_metal-S_o2_gas)*convert_dS
#print (S_iro2_solid-S_ir_metal-S_o2_gas)*T/(1000)
dg_iro2_exp=-188.386/kjmol

dg_iro2=dh_iro2-TdS_iro2
print('IrO2: dG= ' ,dg_iro2, ' eV ', dg_iro2/3, dg_iro2*kjmol, ' kjmol exp: ', dg_iro2_exp*kjmol, ' kjmol')

factor=3/2 #adjusted to reflect more O in IrO3 vs IrO2
TdS_iro3=(S_iro2_solid*factor-S_ir_metal-S_o2_gas*3/2)*convert_dS
dg_iro3=dh_iro3-TdS_iro3
print('IrO3: dG= ' ,dg_iro3, ' eV ', dg_iro3/4, dg_iro3*kjmol, ' kjmol')

factor=4/2 #adjusted to reflect more O H in IriHO3 vs IrO2, phonons needs to be calculated
TdS_irho3=(S_iro2_solid*factor-S_ir_metal-S_o2_gas*3/2-S_h2_gas/2)*convert_dS
#TdS_irho3=(-S_o2_gas*3/2-S_h2_gas/2)*convert_dS
dg_irho3=dh_irho3-TdS_irho3
print('IrHO3: dG= ' ,dg_irho3, ' eV ', dg_irho3/5, dg_irho3*kjmol, ' kjmol')


print('#############################the end ###################################')
