#!/usr/bin/env python

"""IrOx Project.

Author: Raul A. Flores

/global/cscratch1/sd/flores12/IrOx_Project/03_OER_Calc/IrO3/111/01_O_covered/04_oh/02_attempt
"""

#| - Import Modules
from os.path import join as join
#__|

#| - Script Parameters
rt = "/global/cscratch1/sd/flores12/IrOx_Project"
surf_e_rt = "01_surface_calcs"
surf_rt = "02_surface_coverage"
ads_rt = "03_OER_Calc"
surf_cov_e_rt = "07_diff_coverages_term"
#__|

#| - dir_list

dir_list = [

    #| - IrO2 *****************************************************************

    #| - Surface Energies
    # 100
    join(rt, surf_e_rt, "IrO2/100/01_layers"),
    join(rt, surf_e_rt, "IrO2/100/03_layers"),
    join(rt, surf_e_rt, "IrO2/100/05_layers"),
    join(rt, surf_e_rt, "IrO2/100/07_layers"),
    join(rt, surf_e_rt, "IrO2/100/09_layers"),
    join(rt, surf_e_rt, "IrO2/100/11_layers"),
    join(rt, surf_e_rt, "IrO2/100/13_layers"),
    join(rt, surf_e_rt, "IrO2/100/15_layers"),

    # 110
    join(rt, surf_e_rt, "IrO2/110/01_layers"),
    join(rt, surf_e_rt, "IrO2/110/03_layers"),
    join(rt, surf_e_rt, "IrO2/110/05_layers"),
    join(rt, surf_e_rt, "IrO2/110/07_layers"),
    join(rt, surf_e_rt, "IrO2/110/09_layers"),
    join(rt, surf_e_rt, "IrO2/110/11_layers"),

    join(rt, surf_e_rt, "IrO2/110/13_layers"),  # Couldn't converge
    #__|

    #| - Surface Coverage Calculations

    # 100

    # join(rt, surf_rt, "IrO2/100/03_bare"),
    # join(rt, surf_rt, "IrO2/100/01_O_covered/1_2_O_covered"),
    # join(rt, surf_rt, "IrO2/100/01_O_covered/2_2_O_covered"),
    # join(rt, surf_rt, "IrO2/100/02_H_covered/1_2_H_covered"),
    # join(rt, surf_rt, "IrO2/100/02_H_covered/2_2_H_covered"),

    join(rt, surf_rt, "IrO2/100_new/01_O_covered"),
    join(rt, surf_rt, "IrO2/100_new/02_H_covered"),
    join(rt, surf_rt, "IrO2/100_new/03_bare"),


    # 110
    join(rt, surf_rt, "IrO2/110/01_O_covered/1_2_O_covered"),
    join(rt, surf_rt, "IrO2/110/01_O_covered/2_2_O_covered"),
    join(rt, surf_rt, "IrO2/110/02_H_covered/1_2_H_covered"),
    join(rt, surf_rt, "IrO2/110/02_H_covered/2_2_H_covered"),
    join(rt, surf_rt, "IrO2/110/03_bare"),

    #__|

    #| - Adsorbate Calculations

    #| - 100

    #| - O-covered
    join(rt, ads_rt, "IrO2/100/01_O_covered/01_bare/01_attempt"),
    join(rt, ads_rt, "IrO2/100/01_O_covered/01_bare/02_attempt"),


    # *OOH Energetics
    join(rt, ads_rt, "IrO2/100/01_O_covered/02_ooh/01_face_up"),
    join(rt, ads_rt, "IrO2/100/01_O_covered/02_ooh/02_face_down"),
    join(rt, ads_rt, "IrO2/100/01_O_covered/02_ooh/02_face_down_2"),

    join(rt, ads_rt,
        "IrO2/100/01_O_covered/02_ooh/03_deprotonated/01_attempt"),
    join(rt, ads_rt,
        "IrO2/100/01_O_covered/02_ooh/03_deprotonated/02_attempt"),


    # *O
    join(rt, ads_rt, "IrO2/100/01_O_covered/03_o/01_attempt"),


    # *OH
    join(rt, ads_rt, "IrO2/100/01_O_covered/04_oh/01_attempt"),
    join(rt, ads_rt, "IrO2/100/01_O_covered/04_oh/02_attempt"),
    #__|

    #| - H-covered
    join(rt, ads_rt, "IrO2/100/02_H_covered/01_bare"),

    join(rt, ads_rt, "IrO2/100/02_H_covered/02_ooh/01_face_up"),
    join(rt, ads_rt, "IrO2/100/02_H_covered/02_ooh/02_face_down"),
    join(rt, ads_rt, "IrO2/100/02_H_covered/02_ooh/03_deprotonated"),

    join(rt, ads_rt, "IrO2/100/02_H_covered/03_o"),
    join(rt, ads_rt, "IrO2/100/02_H_covered/04_oh"),
    #__|

    #__|

    #| - 110

    #| - O-covered
    join(rt, ads_rt, "IrO2/110/01_O_covered/01_bare"),

    # /global/cscratch1/sd/flores12/IrOx_Project/03_OER_Calc/
    # IrO2/110/01_O_covered/02_ooh/03_deprotonated/01_attempt/1-run

    # *OOH Energetics
    join(rt, ads_rt, "IrO2/110/01_O_covered/02_ooh/01_face_up"),
    join(rt, ads_rt, "IrO2/110/01_O_covered/02_ooh/02_face_down"),
    join(rt, ads_rt,
        "IrO2/110/01_O_covered/02_ooh/03_deprotonated/01_attempt"),

    join(rt, ads_rt, "IrO2/110/01_O_covered/03_o"),
    join(rt, ads_rt, "IrO2/110/01_O_covered/04_oh"),

    #__|

    #| - H-covered
    join(rt, ads_rt, "IrO2/110/02_H_covered/01_bare"),

    join(rt, ads_rt, "IrO2/110/02_H_covered/02_ooh/01_face_up"),
    join(rt, ads_rt, "IrO2/110/02_H_covered/02_ooh/02_face_down"),
    join(rt, ads_rt, "IrO2/110/02_H_covered/02_ooh/02_face_down_2"),
    join(rt, ads_rt, "IrO2/110/02_H_covered/02_ooh/03_deprotonated"),
    join(rt, ads_rt, "IrO2/110/02_H_covered/02_ooh/04_sideways"),

    join(rt, ads_rt, "IrO2/110/02_H_covered/03_o"),
    join(rt, ads_rt, "IrO2/110/02_H_covered/04_oh"),
    #__|

    #__|

    #__|

    #| - Surface Energies @ Different Coverages

    # 100
    join(rt, surf_cov_e_rt, "IrO2/100/OH_covered"),
    join(rt, surf_cov_e_rt, "IrO2/100/O_covered"),
    join(rt, surf_cov_e_rt, "IrO2/100/bare_covered"),

    # 110
    join(rt, surf_cov_e_rt, "IrO2/110/OH_covered"),
    join(rt, surf_cov_e_rt, "IrO2/110/O_covered"),
    join(rt, surf_cov_e_rt, "IrO2/110/bare_covered"),
    #__|

    #__| **********************************************************************

    #| - IrO3 (a-AlF3) ********************************************************

    #| - Surface Energies

    # 100
    join(rt, surf_e_rt, "IrO3/100/03_layers"),
    join(rt, surf_e_rt, "IrO3/100/07_layers"),
    join(rt, surf_e_rt, "IrO3/100/11_layers"),
    join(rt, surf_e_rt, "IrO3/100/15_layers"),
    join(rt, surf_e_rt, "IrO3/100/19_layers"),
    join(rt, surf_e_rt, "IrO3/100/23_layers"),

    join(rt, surf_e_rt, "IrO3/100/05_layers"),
    join(rt, surf_e_rt, "IrO3/100/09_layers"),
    join(rt, surf_e_rt, "IrO3/100/13_layers"),
    join(rt, surf_e_rt, "IrO3/100/17_layers"),
    join(rt, surf_e_rt, "IrO3/100/21_layers"),
    join(rt, surf_e_rt, "IrO3/100/25_layers"),
    join(rt, surf_e_rt, "IrO3/100/29_layers"),
    join(rt, surf_e_rt, "IrO3/100/33_layers"),
    join(rt, surf_e_rt, "IrO3/100/37_layers"),
    join(rt, surf_e_rt, "IrO3/100/41_layers"),

    # 110
    join(rt, surf_e_rt, "IrO3/110/01_layers"),
    join(rt, surf_e_rt, "IrO3/110/02_layers"),
    join(rt, surf_e_rt, "IrO3/110/03_layers"),
    join(rt, surf_e_rt, "IrO3/110/04_layers"),
    join(rt, surf_e_rt, "IrO3/110/05_layers"),
    join(rt, surf_e_rt, "IrO3/110/06_layers"),
    join(rt, surf_e_rt, "IrO3/110/07_layers"),
    join(rt, surf_e_rt, "IrO3/110/08_layers"),
    join(rt, surf_e_rt, "IrO3/110/09_layers"),

    join(rt, surf_e_rt, "IrO3/110/10_layers"),
    join(rt, surf_e_rt, "IrO3/110/11_layers"),
    join(rt, surf_e_rt, "IrO3/110/12_layers"),

    # 111
    join(rt, surf_e_rt, "IrO3/111/01_layers"),
    join(rt, surf_e_rt, "IrO3/111/02_layers"),
    join(rt, surf_e_rt, "IrO3/111/03_layers"),
    join(rt, surf_e_rt, "IrO3/111/04_layers"),
    join(rt, surf_e_rt, "IrO3/111/05_layers"),
    join(rt, surf_e_rt, "IrO3/111/06_layers"),
    join(rt, surf_e_rt, "IrO3/111/07_layers"),
    join(rt, surf_e_rt, "IrO3/111/08_layers"),
    join(rt, surf_e_rt, "IrO3/111/09_layers"),
    join(rt, surf_e_rt, "IrO3/111/10_layers"),
    join(rt, surf_e_rt, "IrO3/111/11_layers"),
    join(rt, surf_e_rt, "IrO3/111/12_layers"),
    join(rt, surf_e_rt, "IrO3/111/13_layers"),
    join(rt, surf_e_rt, "IrO3/111/14_layers"),

    # 211
    join(rt, surf_e_rt, "IrO3/211/01_layers"),
    join(rt, surf_e_rt, "IrO3/211/02_layers"),
    join(rt, surf_e_rt, "IrO3/211/03_layers"),
    join(rt, surf_e_rt, "IrO3/211/04_layers"),
    join(rt, surf_e_rt, "IrO3/211/05_layers"),
    join(rt, surf_e_rt, "IrO3/211/06_layers"),
    join(rt, surf_e_rt, "IrO3/211/07_layers"),
    join(rt, surf_e_rt, "IrO3/211/08_layers"),
    join(rt, surf_e_rt, "IrO3/211/09_layers"),
    join(rt, surf_e_rt, "IrO3/211/10_layers"),
    join(rt, surf_e_rt, "IrO3/211/11_layers"),
    join(rt, surf_e_rt, "IrO3/211/12_layers"),
    join(rt, surf_e_rt, "IrO3/211/13_layers"),
    join(rt, surf_e_rt, "IrO3/211/14_layers"),

    #__|

    #| - Surface Coverage Calculations
    # 100
    join(rt, surf_rt, "IrO3/100/01_O_covered"),
    join(rt, surf_rt, "IrO3/100/02_H_covered"),
    join(rt, surf_rt, "IrO3/100/03_bare"),

    # 110
    join(rt, surf_rt, "IrO3/110/01_O_covered"),
    join(rt, surf_rt, "IrO3/110/02_H_covered"),
    join(rt, surf_rt, "IrO3/110/03_bare"),

    # 111
    join(rt, surf_rt, "IrO3/111/01_O_covered/1_6_O_covered"),
    join(rt, surf_rt, "IrO3/111/01_O_covered/2_6_O_covered"),
    join(rt, surf_rt, "IrO3/111/01_O_covered/3_6_O_covered"),
    join(rt, surf_rt, "IrO3/111/01_O_covered/4_6_O_covered"),
    join(rt, surf_rt, "IrO3/111/01_O_covered/5_6_O_covered"),
    join(rt, surf_rt, "IrO3/111/01_O_covered/6_6_O_covered"),

    join(rt, surf_rt, "IrO3/111/02_H_covered/1_6_H_covered"),
    join(rt, surf_rt, "IrO3/111/02_H_covered/2_6_H_covered"),
    join(rt, surf_rt, "IrO3/111/02_H_covered/3_6_H_covered"),
    join(rt, surf_rt, "IrO3/111/02_H_covered/4_6_H_covered"),
    join(rt, surf_rt, "IrO3/111/02_H_covered/5_6_H_covered"),
    join(rt, surf_rt, "IrO3/111/02_H_covered/6_6_H_covered"),

    join(rt, surf_rt, "IrO3/111/03_bare"),

    # 211
    join(rt, surf_rt, "IrO3/211/01_O_covered"),
    join(rt, surf_rt, "IrO3/211/02_H_covered/01_1st_O_covered"),
    join(rt, surf_rt, "IrO3/211/02_H_covered/02_2nd_O_covered"),
    join(rt, surf_rt, "IrO3/211/03_bare"),

    #__|

    #| - Adsorbate Calculations

    #| - 100

    #| - O-covered
    join(rt, ads_rt, "IrO3/100/01_O_covered/01_bare/01_attempt"),
    join(rt, ads_rt, "IrO3/100/01_O_covered/01_bare/02_attempt"),

    # *OOH Energetics
    join(rt, ads_rt, "IrO3/100/01_O_covered/02_ooh/01_face_up/01_attempt"),
    join(rt, ads_rt, "IrO3/100/01_O_covered/02_ooh/02_face_down/01_attempt"),
    join(rt, ads_rt, "IrO3/100/01_O_covered/02_ooh/02_face_down/02_attempt"),
    join(rt, ads_rt, "IrO3/100/01_O_covered/02_ooh/03_deprotonated/01_attempt"),

    # *O Energetics
    join(rt, ads_rt, "IrO3/100/01_O_covered/03_o/01_attempt"),
    join(rt, ads_rt, "IrO3/100/01_O_covered/03_o/02_attempt"),

    # *OH Energetics
    join(rt, ads_rt, "IrO3/100/01_O_covered/04_oh/01_attempt"),
    join(rt, ads_rt, "IrO3/100/01_O_covered/04_oh/02_attempt"),
    #__|

    #| - H-covered
    join(rt, ads_rt, "IrO3/100/02_H_covered/01_bare"),

    join(rt, ads_rt, "IrO3/100/02_H_covered/02_ooh/01_face_up"),
    join(rt, ads_rt, "IrO3/100/02_H_covered/02_ooh/02_face_down"),
    join(rt, ads_rt,
        "IrO3/100/02_H_covered/02_ooh/03_deprotonated/01_attempt"),

    # /global/cscratch1/sd/flores12/IrOx_Project/03_OER_Calc/
    # IrO3/100/02_H_covered/02_ooh/03_deprotonated/01_attempt/1-run

    join(rt, ads_rt, "IrO3/100/02_H_covered/03_o"),
    join(rt, ads_rt, "IrO3/100/02_H_covered/04_oh"),
    #__|

    #__|

    #| - 110

    #| - O-covered
    join(rt, ads_rt, "IrO3/110/01_O_covered/01_bare/01_attempt"),
    join(rt, ads_rt, "IrO3/110/01_O_covered/01_bare/02_attempt"),

    # *OOH Energetics
    join(rt, ads_rt, "IrO3/110/01_O_covered/02_ooh/01_face_up/01_attempt"),
    join(rt, ads_rt, "IrO3/110/01_O_covered/02_ooh/02_face_down/01_attempt"),
    join(rt, ads_rt, "IrO3/110/01_O_covered/02_ooh/02_face_down/02_attempt"),
    join(rt, ads_rt, "IrO3/110/01_O_covered/02_ooh/03_deprotonated/01_attempt"),

    # *O Energetics
    join(rt, ads_rt, "IrO3/110/01_O_covered/03_o/01_attempt"),
    join(rt, ads_rt, "IrO3/110/01_O_covered/03_o/02_attempt"),

    # *OH Energetics
    join(rt, ads_rt, "IrO3/110/01_O_covered/04_oh/01_attempt"),
    join(rt, ads_rt, "IrO3/110/01_O_covered/04_oh/02_attempt"),
    #__|

    #| - H-covered
    join(rt, ads_rt, "IrO3/110/02_H_covered/01_bare"),

    join(rt, ads_rt, "IrO3/110/02_H_covered/02_ooh/01_face_up"),
    join(rt, ads_rt, "IrO3/110/02_H_covered/02_ooh/02_face_down"),
    join(rt, ads_rt, "IrO3/110/02_H_covered/02_ooh/03_deprotonated"),
    join(rt, ads_rt, "IrO3/110/02_H_covered/02_ooh/04_sideways"),

    join(rt, ads_rt, "IrO3/110/02_H_covered/03_o"),
    join(rt, ads_rt, "IrO3/110/02_H_covered/04_oh"),
    #__|

    #__|

    #| - 111

    #| - O-covered
    join(rt, ads_rt, "IrO3/111/01_O_covered/01_bare/01_attempt"),
    join(rt, ads_rt, "IrO3/111/01_O_covered/01_bare/02_attempt"),

    # *OOH Energetics
    join(rt, ads_rt, "IrO3/111/01_O_covered/02_ooh/01_face_up/01_attempt"),
    join(rt, ads_rt, "IrO3/111/01_O_covered/02_ooh/02_face_down/01_attempt"),

    join(rt, ads_rt, "IrO3/111/01_O_covered/02_ooh/03_deprotonated/01_attempt"),

    # *O Energetics
    join(rt, ads_rt, "IrO3/111/01_O_covered/03_o/01_attempt"),
    join(rt, ads_rt, "IrO3/111/01_O_covered/03_o/02_attempt"),

    # *OH Energetics
    join(rt, ads_rt, "IrO3/111/01_O_covered/04_oh/01_attempt"),
    join(rt, ads_rt, "IrO3/111/01_O_covered/04_oh/02_attempt"),
    #__|

    #| - H-covered | Haven't done yet
    # join(rt, ads_rt, "IrO3/111/02_H_covered/01_bare"),
    # join(rt, ads_rt, "IrO3/111/02_H_covered/02_ooh"),
    # join(rt, ads_rt, "IrO3/111/02_H_covered/03_o"),
    # join(rt, ads_rt, "IrO3/111/02_H_covered/04_oh"),
    #__|

    #__|

    #| - 211

    #| - O-covered
    join(rt, ads_rt, "IrO3/211/01_O_covered/01_bare/01_attempt"),
    join(rt, ads_rt, "IrO3/211/01_O_covered/01_bare/02_attempt"),

    join(rt, ads_rt, "IrO3/211/01_O_covered/02_ooh/01_face_up/01_attempt"),
    join(rt, ads_rt, "IrO3/211/01_O_covered/02_ooh/02_face_down/01_attempt"),
    join(rt, ads_rt, "IrO3/211/01_O_covered/02_ooh/03_deprotonated/01_attempt"),

    join(rt, ads_rt, "IrO3/211/01_O_covered/03_o/01_attempt"),
    join(rt, ads_rt, "IrO3/211/01_O_covered/03_o/02_attempt"),

    join(rt, ads_rt, "IrO3/211/01_O_covered/04_oh/01_attempt"),
    join(rt, ads_rt, "IrO3/211/01_O_covered/04_oh/02_attempt"),
    #__|

    #| - H-covered
    join(rt, ads_rt, "IrO3/211/02_H_covered/01_bare"),


    join(rt, ads_rt, "IrO3/211/02_H_covered/02_ooh/01_face_up"),
    join(rt, ads_rt, "IrO3/211/02_H_covered/02_ooh/02_face_down"),


    join(rt, ads_rt, "IrO3/211/02_H_covered/03_o"),
    join(rt, ads_rt, "IrO3/211/02_H_covered/04_oh"),
    #__|

    #__|

    #__|

    #| - Surface Energies @ Different Coverages
    # 100
    join(rt, surf_cov_e_rt, "IrO3/100/OH_covered"),
    join(rt, surf_cov_e_rt, "IrO3/100/O_covered"),
    join(rt, surf_cov_e_rt, "IrO3/100/bare_covered"),

    # 110
    join(rt, surf_cov_e_rt, "IrO3/110/OH_covered"),
    join(rt, surf_cov_e_rt, "IrO3/110/O_covered"),
    join(rt, surf_cov_e_rt, "IrO3/110/bare_covered"),

    # 111
    join(rt, surf_cov_e_rt, "IrO3/111/OH_covered"),
    join(rt, surf_cov_e_rt, "IrO3/111/O_covered"),
    join(rt, surf_cov_e_rt, "IrO3/111/bare_covered"),

    # 211
    join(rt, surf_cov_e_rt, "IrO3/211/OH_covered"),
    join(rt, surf_cov_e_rt, "IrO3/211/O_covered"),
    join(rt, surf_cov_e_rt, "IrO3/211/bare_covered"),

    #__|

    #__| **********************************************************************

    #| - IrO3 (rutile-like) ***************************************************

    #| - Surface Energies

    # 001
    join(rt, surf_e_rt, "IrO3_rutile-like/001/02_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/001/06_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/001/10_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/001/14_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/001/18_layers"),

    # 100
    join(rt, surf_e_rt, "IrO3_rutile-like/100/01_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/100/03_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/100/05_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/100/07_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/100/09_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/100/11_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/100/13_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/100/15_layers"),

    # 110
    join(rt, surf_e_rt, "IrO3_rutile-like/110/01_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/110/03_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/110/05_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/110/07_layers"),
    join(rt, surf_e_rt, "IrO3_rutile-like/110/09_layers"),

    #__|

    #| - Surface Coverage Calculations

    # 001
    join(rt, surf_rt, "IrO3_rutile-like/001/01_O_covered/0_8_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/01_O_covered/1_8_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/01_O_covered/2_8_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/01_O_covered/3_8_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/01_O_covered/4_8_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/01_O_covered/5_8_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/01_O_covered/6_8_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/01_O_covered/7_8_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/01_O_covered/8_8_O_covered"),

    join(rt, surf_rt, "IrO3_rutile-like/001/02_H_covered/0_8_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/02_H_covered/1_8_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/02_H_covered/2_8_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/02_H_covered/3_8_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/02_H_covered/4_8_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/02_H_covered/5_8_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/02_H_covered/6_8_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/02_H_covered/7_8_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/001/02_H_covered/8_8_H_covered"),

    # 100
    join(rt, surf_rt, "IrO3_rutile-like/100/01_O_covered/0_4_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/100/01_O_covered/1_4_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/100/01_O_covered/2_4_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/100/01_O_covered/3_4_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/100/01_O_covered/4_4_O_covered"),

    join(rt, surf_rt, "IrO3_rutile-like/100/02_H_covered/0_4_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/100/02_H_covered/1_4_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/100/02_H_covered/2_4_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/100/02_H_covered/3_4_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/100/02_H_covered/4_4_H_covered"),

    # 110
    join(rt, surf_rt, "IrO3_rutile-like/110/01_O_covered/0_2_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/110/01_O_covered/1_2_O_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/110/01_O_covered/2_2_O_covered"),

    join(rt, surf_rt, "IrO3_rutile-like/110/02_H_covered/0_2_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/110/02_H_covered/1_2_H_covered"),
    join(rt, surf_rt, "IrO3_rutile-like/110/02_H_covered/2_2_H_covered"),

    #__|

    #| - Adsorbate Calculations

    # TODO
    #| - 001

    #__|

    #| - 100

    #| - 01_O-4_OH-0
    join(rt, ads_rt, "IrO3_rutile-like/100/01_O-4_OH-0/01_bare/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/100/01_O-4_OH-0/01_bare/02_attempt"),

    # *OOH Energetics
    join(rt, ads_rt, "IrO3_rutile-like/100/01_O-4_OH-0/02_ooh/01_face_up/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/100/01_O-4_OH-0/02_ooh/02_face_down/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/100/01_O-4_OH-0/02_ooh/03_deprotonated/01_attempt"),

    join(rt, ads_rt, "IrO3_rutile-like/100/01_O-4_OH-0/03_o/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/100/01_O-4_OH-0/03_o/02_attempt"),

    join(rt, ads_rt, "IrO3_rutile-like/100/01_O-4_OH-0/04_oh/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/100/01_O-4_OH-0/04_oh/02_attempt"),
    #__|

    #| - 02_O-2_OH-0
    join(rt, ads_rt, "IrO3_rutile-like/100/02_O-2_OH-0/01_bare/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/100/02_O-2_OH-0/01_bare/02_attempt"),


    # *OOH Energetics
    join(rt, ads_rt, "IrO3_rutile-like/100/02_O-2_OH-0/02_ooh", "01_face_up/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/100/02_O-2_OH-0/02_ooh", "02_face_down/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/100/02_O-2_OH-0/02_ooh", "02_face_down/02_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/100/02_O-2_OH-0/02_ooh", "03_deprotonated/01_attempt"),


    join(rt, ads_rt, "IrO3_rutile-like/100/02_O-2_OH-0/03_o/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/100/02_O-2_OH-0/03_o/02_attempt"),

    join(rt, ads_rt, "IrO3_rutile-like/100/02_O-2_OH-0/04_oh/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/100/02_O-2_OH-0/04_oh/02_attempt"),
    #__|

    #| - 03_O-2_OH-2
    join(rt, ads_rt, "IrO3_rutile-like/100/03_O-2_OH-2/01_bare"),
    join(rt, ads_rt, "IrO3_rutile-like/100/03_O-2_OH-2/02_ooh"),
    join(rt, ads_rt, "IrO3_rutile-like/100/03_O-2_OH-2/03_o"),
    join(rt, ads_rt, "IrO3_rutile-like/100/03_O-2_OH-2/04_oh"),
    #__|

    #__|

    #| - 110

    #| - O_covered
    join(rt, ads_rt, "IrO3_rutile-like/110/01_O_covered/01_bare/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/110/01_O_covered/01_bare/02_attempt"),


    # *OOH Energetics
    join(rt, ads_rt, "IrO3_rutile-like/110/01_O_covered/02_ooh", "01_face_up/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/110/01_O_covered/02_ooh", "02_face_down/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/110/01_O_covered/02_ooh", "03_deprotonated/01_attempt"),


    join(rt, ads_rt, "IrO3_rutile-like/110/01_O_covered/03_o/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/110/01_O_covered/03_o/02_attempt"),

    join(rt, ads_rt, "IrO3_rutile-like/110/01_O_covered/04_oh/01_attempt"),
    join(rt, ads_rt, "IrO3_rutile-like/110/01_O_covered/04_oh/02_attempt"),
    #__|

    #| - H_covered
    join(rt, ads_rt, "IrO3_rutile-like/110/02_H_covered/01_bare"),
    join(rt, ads_rt, "IrO3_rutile-like/110/02_H_covered/02_ooh"),
    join(rt, ads_rt, "IrO3_rutile-like/110/02_H_covered/03_o"),
    join(rt, ads_rt, "IrO3_rutile-like/110/02_H_covered/04_oh"),
    #__|

    #__|

    #__|

    #| - Surface Energies @ Different Coverages

    # 001
    join(rt, surf_cov_e_rt, "IrO3_rutile-like/001/OH_covered"),
    join(rt, surf_cov_e_rt, "IrO3_rutile-like/001/O_covered"),
    join(rt, surf_cov_e_rt, "IrO3_rutile-like/001/bare_covered"),

    # 100
    join(rt, surf_cov_e_rt, "IrO3_rutile-like/100/OH_covered"),
    join(rt, surf_cov_e_rt, "IrO3_rutile-like/100/O_covered"),
    join(rt, surf_cov_e_rt, "IrO3_rutile-like/100/bare_covered"),

    # 110
    join(rt, surf_cov_e_rt, "IrO3_rutile-like/110/OH_covered"),
    join(rt, surf_cov_e_rt, "IrO3_rutile-like/110/O_covered"),
    join(rt, surf_cov_e_rt, "IrO3_rutile-like/110/bare_covered"),

    #__|

    #__|

    #| - IrO3 (battery) *******************************************************

    #| - Adsorbate Calculations

    #| - 001

    #__|

    #| - 010

    #| - 01_surface_type_a

    #| - O-covered
    join(rt, ads_rt, "IrO3_battery/010/01_surface_type_a/01_O_covered/01_bare/01_attempt"),
    join(rt, ads_rt, "IrO3_battery/010/01_surface_type_a/01_O_covered/01_bare/02_attempt"),


    # *OOH Energetics
    join(rt, ads_rt, "IrO3_battery/010/01_surface_type_a/01_O_covered/02_ooh/01_face_up/01_attempt"),
    join(rt, ads_rt, "IrO3_battery/010/01_surface_type_a/01_O_covered/02_ooh/02_face_down/01_attempt"),
    join(rt, ads_rt, "IrO3_battery/010/01_surface_type_a/01_O_covered/02_ooh/03_deprotonated/01_attempt"),


    join(rt, ads_rt, "IrO3_battery/010/01_surface_type_a/01_O_covered/03_o/01_attempt"),
    join(rt, ads_rt, "IrO3_battery/010/01_surface_type_a/01_O_covered/03_o/02_attempt"),

    join(rt, ads_rt, "IrO3_battery/010/01_surface_type_a/01_O_covered/04_oh/01_attempt"),
    join(rt, ads_rt, "IrO3_battery/010/01_surface_type_a/01_O_covered/04_oh/02_attempt"),
    #__|

    #__|

    #| - 02_surface_type_b

    #| - O-covered
    join(rt, ads_rt, "IrO3_battery/010/01_O_covered/01_bare/01_attempt"),
    join(rt, ads_rt, "IrO3_battery/010/01_O_covered/01_bare/02_attempt"),


    # *OOH Energetics
    join(rt, ads_rt, "IrO3_battery/010/01_O_covered/02_ooh/01_face_up/01_attempt"),
    join(rt, ads_rt, "IrO3_battery/010/01_O_covered/02_ooh/02_face_down/01_attempt"),
    join(rt, ads_rt, "IrO3_battery/010/01_O_covered/02_ooh/03_deprotonated/01_attempt"),


    join(rt, ads_rt, "IrO3_battery/010/01_O_covered/03_o/01_attempt"),
    join(rt, ads_rt, "IrO3_battery/010/01_O_covered/03_o/02_attempt"),


    join(rt, ads_rt, "IrO3_battery/010/01_O_covered/04_oh/01_attempt"),
    join(rt, ads_rt, "IrO3_battery/010/01_O_covered/04_oh/02_attempt"),
    #__|

    #__|

    #__|

    #__|

    #| - Surface Energies @ Different Coverages

    #| - 001
    join(rt, surf_cov_e_rt, "IrO3_battery/001/01_surface_type_a/12_layers/O_covered/01_attempt"),
    # join(rt, surf_cov_e_rt, ""),
    #__|

    #| - 010
    # join(rt, surf_cov_e_rt, "IrO3_battery/010/01_surface_type_a/04_layers/O_covered"),

    join(rt, surf_cov_e_rt, "IrO3_battery/010/01_surface_type_a/02_layers/O_covered/01_attempt"),
    join(rt, surf_cov_e_rt, "IrO3_battery/010/01_surface_type_a/02_layers/half_O_covered/01_attempt"),
    join(rt, surf_cov_e_rt, "IrO3_battery/010/01_surface_type_a/02_layers/OH_covered/01_attempt"),
    join(rt, surf_cov_e_rt, "IrO3_battery/010/01_surface_type_a/02_layers/bare_covered/01_attempt"),

    # join(rt, surf_cov_e_rt, ""),
    #__|

    #__|

    #__| **********************************************************************

    ]

#__|
