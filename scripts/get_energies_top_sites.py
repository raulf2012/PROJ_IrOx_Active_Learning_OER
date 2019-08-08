#!/bin/bash

h2o='-14.23091949'

clean=$(grep ' "energy": ' bare.json | gawk '{print substr($2, 1, length($2)-1)}')
oh=$(grep ' "energy": ' oh.json | gawk '{print substr($2, 1, length($2)-1)}')
o=$(grep ' "energy": ' o.json | gawk '{print substr($2, 1, length($2)-1)}')
ooh_up=$(grep ' "energy": ' ooh_up.json | gawk '{print substr($2, 1, length($2)-1)}')
ooh_down=$(grep ' "energy": ' ooh_down.json | gawk '{print substr($2, 1, length($2)-1)}')

# ########################
echo $clean $oh $o $ooh

echo "ooh_up"
echo $ooh_up

echo "ooh_down"
echo $ooh_down

echo ""

echo 'OOH up'
# python ~/bin/get_overpot_MB_MG_vasp_TS_vers2.py -raw  $clean $oh $o $ooh_up
python $PROJ_irox/scripts/02_calc_overpot/get_overpot_MB_MG_vasp_TS_vers2.py -raw  $clean $oh $o $ooh_up

echo 'OOH down'
# python ~/bin/get_overpot_MB_MG_vasp_TS_vers2.py -raw  $clean $oh $o $ooh_down
python $PROJ_irox/scripts/02_calc_overpot/get_overpot_MB_MG_vasp_TS_vers2.py -raw  $clean $oh $o $ooh_down
