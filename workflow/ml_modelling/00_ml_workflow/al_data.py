#!/usr/bin/env python

"""TEMP.

Author(s): Raul A. Flores
"""


#| - Import Modules
import os
import sys
#__|



#| - AL Settings and Variables

main_AB2_run = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/01_abx_al_runs_new/out_data/AB2/gp_ucb_True",
    "AB2_AL_02_ropeseso.pickle")

    #| - __old__
    #  "workflow/ml_modelling/00_ml_workflow",
    #  "191102_new_workflow/01_abx_al_runs_new/out_data/AB2/gp_ucb_True",
    #  "TEST_AL_6_radetama.pickle",


    #  "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/01_abx_al_runs_new/out_data/AB2/gp_ucb_True",
    #  "TEST_AL_6_radetama.pickle",
    #__|

main_AB3_run = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB3/gp_ucb_True",
    "TEST_AL_2_fugunefo.pickle")

    #| - __old__
    # "01_attempt/AL_geheneva.pickle",
    # "01_attempt/AL_pifehohu.pickle",
    # "01_attempt/AL_suturomo.pickle",

    # NEW RUNS
    # "TEST_AL_2_fugunefo.pickle",
    # "TEST_AL_2_seruladi.pickle",
    # "TEST_AL_masahiti.pickle",  # This is a good one to pick ************

    # #################################################
    #  'workflow/ml_modelling/00_ml_workflow',
    #  '191102_new_workflow/00_abx_al_runs',
    #  'out_data/AB3/gp_ucb_True',
    #  'AL_kupugeti.pickle',
    #__|

#__|


# | - Generations to plot
# generations to plot

#| - AB3
# if stoich_i == "AB3":

#     gen_0 = 0
#     gen_1 = 2
#     gen_2 = 9  # 5/10 top structures found
#     gen_3 = 14 # 10/10 top structures found
#     gen_4 = "last"

#     gen_0 = 0
#     gen_1 = 3
#     gen_2 = 6  # 7/10 top structures found
#     gen_3 = 20 # 10/10 top structures found
#     gen_4 = -10

gen_0 = 0   # 00/10
gen_1 = 3   # 00/10
gen_2 = 5   # 04/10
gen_3 = 14  # 09/10
gen_4 = 20 # 10/10

gens_to_plot_ab3 = [gen_0, gen_1, gen_2, gen_3, gen_4]
#__|

#| - AB2
# AB2  # 92 total gens
# 01 found | gen 3
# 03 found | gen 5
# 04 found | gen 9
# 06 found | gen 12
# 07 found | gen 15
# 08 found | gen 16
# 09 found | gen 27
# 10 found | gen 75


# if stoich_i == "AB2":
gen_0 = 0
gen_1 = 5
gen_2 = 12
gen_3 = 16
gen_4 = 27

gens_to_plot_ab2 = [gen_0, gen_1, gen_2, gen_3, gen_4]

main_gen = gen_2


gens_to_plot_dict = {
    "AB2": gens_to_plot_ab2,
    "AB3": gens_to_plot_ab3,
    }

#__|

#__|


#| - Files to read ############################################################
pre_path_all = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow",
    )

#| - AB3

#| - Random Aquisition

pre_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB3/random_True",
    )

files_list_ab3_random = [
    # pre_path + "/AL_duhakiro.pickle",
    # pre_path + "/AL_firegohi.pickle",
    # pre_path + "/AL_givohegu.pickle",
    # pre_path + "/AL_tiweluku.pickle",
    # pre_path + "/AL_vevenuwa.pickle",

    pre_path + "/TEST_AL_2_fibataha.pickle",
    pre_path + "/TEST_AL_2_kitagego.pickle",
    pre_path + "/TEST_AL_3_bafigika.pickle",
    pre_path + "/TEST_AL_3_bafomepo.pickle",
    pre_path + "/TEST_AL_3_besurogi.pickle",
    pre_path + "/TEST_AL_3_dotupibu.pickle",
    pre_path + "/TEST_AL_3_fuliguso.pickle",
    pre_path + "/TEST_AL_3_gosirinu.pickle",
    pre_path + "/TEST_AL_3_gowutoga.pickle",
    pre_path + "/TEST_AL_3_hiroguwe.pickle",
    pre_path + "/TEST_AL_3_kavasaki.pickle",
    pre_path + "/TEST_AL_3_mutefoti.pickle",
    pre_path + "/TEST_AL_3_nivalula.pickle",
    pre_path + "/TEST_AL_3_pasimuha.pickle",
    pre_path + "/TEST_AL_3_sukiwalo.pickle",
    pre_path + "/TEST_AL_3_temonofo.pickle",
    pre_path + "/TEST_AL_3_vesudewa.pickle",
    pre_path + "/TEST_AL_3_vuwugupi.pickle",
    pre_path + "/TEST_AL_3_walipebi.pickle",
    pre_path + "/TEST_AL_3_wetipotu.pickle",
    pre_path + "/TEST_AL_3_wusabupa.pickle",
    pre_path + "/TEST_AL_3_wutonovi.pickle",

    pre_path + "/TEST_AL_4_benegeka.pickle",
    pre_path + "/TEST_AL_4_dehebiko.pickle",
    pre_path + "/TEST_AL_4_dinisefa.pickle",
    pre_path + "/TEST_AL_4_fefefigi.pickle",
    pre_path + "/TEST_AL_4_fefesama.pickle",
    pre_path + "/TEST_AL_4_fivokito.pickle",
    pre_path + "/TEST_AL_4_fuwasufi.pickle",
    pre_path + "/TEST_AL_4_gekuporu.pickle",
    pre_path + "/TEST_AL_4_gepapeba.pickle",
    pre_path + "/TEST_AL_4_gerisiwe.pickle",
    pre_path + "/TEST_AL_4_goderiwo.pickle",
    pre_path + "/TEST_AL_4_hofavasi.pickle",
    pre_path + "/TEST_AL_4_kavadosu.pickle",
    pre_path + "/TEST_AL_4_kudadega.pickle",
    pre_path + "/TEST_AL_4_masufika.pickle",
    pre_path + "/TEST_AL_4_menireve.pickle",
    pre_path + "/TEST_AL_4_metekovo.pickle",
    pre_path + "/TEST_AL_4_mibunova.pickle",
    pre_path + "/TEST_AL_4_milesumi.pickle",
    pre_path + "/TEST_AL_4_nepubene.pickle",
    pre_path + "/TEST_AL_4_nuvopeki.pickle",
    pre_path + "/TEST_AL_4_petosuso.pickle",
    pre_path + "/TEST_AL_4_pokikugi.pickle",
    pre_path + "/TEST_AL_4_ragifipa.pickle",
    pre_path + "/TEST_AL_4_rodadibi.pickle",
    pre_path + "/TEST_AL_4_rovukuma.pickle",
    pre_path + "/TEST_AL_4_sanaruri.pickle",
    pre_path + "/TEST_AL_4_sapaveme.pickle",
    pre_path + "/TEST_AL_4_sawuhewe.pickle",
    pre_path + "/TEST_AL_4_sifisulu.pickle",
    pre_path + "/TEST_AL_4_sifobuwe.pickle",
    pre_path + "/TEST_AL_4_suhegope.pickle",
    pre_path + "/TEST_AL_4_teruwufo.pickle",
    pre_path + "/TEST_AL_4_tikeluvi.pickle",
    pre_path + "/TEST_AL_4_vabesipo.pickle",
    pre_path + "/TEST_AL_4_visisese.pickle",
    pre_path + "/TEST_AL_4_wemawumu.pickle",
    pre_path + "/TEST_AL_4_wibemafu.pickle",
    pre_path + "/TEST_AL_4_wobumone.pickle",
    pre_path + "/TEST_AL_4_wohaguha.pickle",

    pre_path + "/TEST_AL_5_bimugese.pickle",
    pre_path + "/TEST_AL_5_hufalage.pickle",
    pre_path + "/TEST_AL_5_humekofa.pickle",
    pre_path + "/TEST_AL_5_kebufetu.pickle",
    pre_path + "/TEST_AL_5_lewehiho.pickle",
    pre_path + "/TEST_AL_5_morosebi.pickle",
    pre_path + "/TEST_AL_5_nehemahe.pickle",
    pre_path + "/TEST_AL_5_remimafu.pickle",
    pre_path + "/TEST_AL_5_siwodonu.pickle",
    pre_path + "/TEST_AL_5_tufekudi.pickle",
    pre_path + "/TEST_AL_5_vapasuda.pickle",
    pre_path + "/TEST_AL_5_watimeni.pickle",
    pre_path + "/TEST_AL_5_wogibalu.pickle",

    pre_path + "/TEST_AL_6_budowilu.pickle",
    pre_path + "/TEST_AL_6_dusideti.pickle",
    pre_path + "/TEST_AL_6_fipawefu.pickle",
    pre_path + "/TEST_AL_6_fisisuhi.pickle",
    pre_path + "/TEST_AL_6_giwomedi.pickle",
    pre_path + "/TEST_AL_6_gosuligu.pickle",
    pre_path + "/TEST_AL_6_guvisabu.pickle",
    pre_path + "/TEST_AL_6_kemelofu.pickle",
    pre_path + "/TEST_AL_6_kihukava.pickle",
    pre_path + "/TEST_AL_6_lonabovu.pickle",
    pre_path + "/TEST_AL_6_pasavona.pickle",
    pre_path + "/TEST_AL_6_ridagimi.pickle",
    pre_path + "/TEST_AL_6_rotadomo.pickle",
    pre_path + "/TEST_AL_6_tabotefi.pickle",
    pre_path + "/TEST_AL_6_tagonowa.pickle",
    pre_path + "/TEST_AL_6_tahavomu.pickle",
    pre_path + "/TEST_AL_6_tumifehu.pickle",
    pre_path + "/TEST_AL_6_vetabari.pickle",
    pre_path + "/TEST_AL_6_wufubusi.pickle",
    pre_path + "/TEST_AL_6_wurapali.pickle",

    ]

#__|

#| - GP-UCB Aquisition
pre_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB3/gp_ucb_True",
    )


files_list_ab3_gp_ucb = [
    pre_path + "/01_attempt/AL_geheneva.pickle",
    pre_path + "/01_attempt/AL_nisoponi.pickle",
    pre_path + "/01_attempt/AL_pifehohu.pickle",
    pre_path + "/01_attempt/AL_suturomo.pickle",
    pre_path + "/01_attempt/AL_vobifoko.pickle",
    pre_path + "/TEST_AL_masahiti.pickle",

    # pre_path + "/TEST_AL_2_devehowo.pickle",  # Not a good run for some reason
    pre_path + "/TEST_AL_2_fugunefo.pickle",
    pre_path + "/TEST_AL_2_hilerika.pickle",
    pre_path + "/TEST_AL_2_pomogobu.pickle",
    pre_path + "/TEST_AL_2_seruladi.pickle",


    pre_path + "/TEST_AL_3_bikufupi.pickle",
    pre_path + "/TEST_AL_3_dakubiku.pickle",
    pre_path + "/TEST_AL_3_duloputo.pickle",
    pre_path + "/TEST_AL_3_fahovara.pickle",
    pre_path + "/TEST_AL_3_fulomoto.pickle",
    pre_path + "/TEST_AL_3_laburike.pickle",
    pre_path + "/TEST_AL_3_libapidi.pickle",
    pre_path + "/TEST_AL_3_nikimido.pickle",
    pre_path + "/TEST_AL_3_raluduhu.pickle",
    pre_path + "/TEST_AL_3_supemono.pickle",

    pre_path + "/TEST_AL_4beradeka.pickle",
    pre_path + "/TEST_AL_4buruduwe.pickle",
    pre_path + "/TEST_AL_4degekoku.pickle",
    pre_path + "/TEST_AL_4deromeru.pickle",
    pre_path + "/TEST_AL_4forafago.pickle",
    pre_path + "/TEST_AL_4hefehepa.pickle",
    pre_path + "/TEST_AL_4kaveboma.pickle",
    pre_path + "/TEST_AL_4kihalage.pickle",
    pre_path + "/TEST_AL_4lidirope.pickle",
    pre_path + "/TEST_AL_4mebetige.pickle",
    pre_path + "/TEST_AL_4megimodi.pickle",
    pre_path + "/TEST_AL_4mohomato.pickle",
    pre_path + "/TEST_AL_4moponuso.pickle",
    pre_path + "/TEST_AL_4mukigapa.pickle",
    pre_path + "/TEST_AL_4nekumuno.pickle",
    pre_path + "/TEST_AL_4nipagula.pickle",
    pre_path + "/TEST_AL_4nupetofi.pickle",
    pre_path + "/TEST_AL_4rifibume.pickle",
    pre_path + "/TEST_AL_4somageho.pickle",
    pre_path + "/TEST_AL_4wolewoba.pickle",

    pre_path + "/TEST_AL_5bofufada.pickle",
    pre_path + "/TEST_AL_5delepaku.pickle",
    pre_path + "/TEST_AL_5derohebi.pickle",
    pre_path + "/TEST_AL_5dotesiga.pickle",
    pre_path + "/TEST_AL_5dumokeru.pickle",
    pre_path + "/TEST_AL_5fapehudo.pickle",
    pre_path + "/TEST_AL_5fenumanu.pickle",
    pre_path + "/TEST_AL_5gomememu.pickle",
    pre_path + "/TEST_AL_5gulumobu.pickle",
    pre_path + "/TEST_AL_5huwihime.pickle",
    pre_path + "/TEST_AL_5kanototo.pickle",
    pre_path + "/TEST_AL_5kibapeto.pickle",
    pre_path + "/TEST_AL_5kifomimi.pickle",
    pre_path + "/TEST_AL_5rukeraku.pickle",
    pre_path + "/TEST_AL_5sipeteni.pickle",
    pre_path + "/TEST_AL_5sokepefo.pickle",
    pre_path + "/TEST_AL_5vikowagu.pickle",
    pre_path + "/TEST_AL_5volibavo.pickle",
    pre_path + "/TEST_AL_5wigidipu.pickle",
    pre_path + "/TEST_AL_5wolowewu.pickle",

    pre_path + "/TEST_AL_6_dodemuho.pickle",
    pre_path + "/TEST_AL_6_fisopova.pickle",
    pre_path + "/TEST_AL_6_gelabere.pickle",
    pre_path + "/TEST_AL_6_gelenuni.pickle",
    pre_path + "/TEST_AL_6_haligagu.pickle",
    pre_path + "/TEST_AL_6_higusare.pickle",
    pre_path + "/TEST_AL_6_kagesinu.pickle",
    pre_path + "/TEST_AL_6_liwuderu.pickle",
    pre_path + "/TEST_AL_6_lopukipu.pickle",
    pre_path + "/TEST_AL_6_nubinada.pickle",
    pre_path + "/TEST_AL_6_pitumimi.pickle",
    pre_path + "/TEST_AL_6_pudapepi.pickle",
    pre_path + "/TEST_AL_6_rividisi.pickle",
    pre_path + "/TEST_AL_6_rovomama.pickle",
    pre_path + "/TEST_AL_6_sovawafo.pickle",
    pre_path + "/TEST_AL_6_takoliwo.pickle",
    pre_path + "/TEST_AL_6_tehileme.pickle",
    pre_path + "/TEST_AL_6_tibuduka.pickle",
    pre_path + "/TEST_AL_6_warusiha.pickle",
    pre_path + "/TEST_AL_6_wevuwofu.pickle",

    pre_path + "/TEST_AL_7_begopuwa.pickle",
    pre_path + "/TEST_AL_7_dakonene.pickle",
    pre_path + "/TEST_AL_7_depapuwa.pickle",
    pre_path + "/TEST_AL_7_dimemihu.pickle",
    pre_path + "/TEST_AL_7_dosegasu.pickle",
    pre_path + "/TEST_AL_7_dubuseko.pickle",
    pre_path + "/TEST_AL_7_ganefusi.pickle",
    pre_path + "/TEST_AL_7_gitewora.pickle",
    pre_path + "/TEST_AL_7_gofihedi.pickle",
    pre_path + "/TEST_AL_7_howurohi.pickle",
    pre_path + "/TEST_AL_7_kunuhefa.pickle",
    pre_path + "/TEST_AL_7_muhubila.pickle",
    pre_path + "/TEST_AL_7_nudevima.pickle",
    pre_path + "/TEST_AL_7_pavonodi.pickle",
    pre_path + "/TEST_AL_7_rurawovo.pickle",
    pre_path + "/TEST_AL_7_sisisodo.pickle",
    pre_path + "/TEST_AL_7_solefata.pickle",
    pre_path + "/TEST_AL_7_tolamika.pickle",
    pre_path + "/TEST_AL_7_vifupebe.pickle",
    pre_path + "/TEST_AL_7_viwovave.pickle",

    ]
#__|

#__|

#| - AB2

#| - Random Aquisition
pre_path = os.path.join(
    pre_path_all,
    "00_abx_al_runs/out_data/AB2/random_True")

pre_path_new = os.path.join(
    pre_path_all,
    "01_abx_al_runs_new/out_data/AB2/random_True")

files_list_ab2_random = [

    pre_path_new + "/AB2_AL_00_didilavo.pickle",
    pre_path_new + "/AB2_AL_00_dipuwera.pickle",
    pre_path_new + "/AB2_AL_00_fakalawu.pickle",
    pre_path_new + "/AB2_AL_00_falokeko.pickle",
    pre_path_new + "/AB2_AL_00_ginukuvi.pickle",
    pre_path_new + "/AB2_AL_00_gufomobu.pickle",
    pre_path_new + "/AB2_AL_00_havonuso.pickle",
    pre_path_new + "/AB2_AL_00_hulifofi.pickle",
    pre_path_new + "/AB2_AL_00_netomewu.pickle",
    pre_path_new + "/AB2_AL_00_novusabu.pickle",
    pre_path_new + "/AB2_AL_00_pebefepo.pickle",
    pre_path_new + "/AB2_AL_00_redokefe.pickle",
    pre_path_new + "/AB2_AL_00_resadiso.pickle",
    pre_path_new + "/AB2_AL_00_tatefulo.pickle",
    pre_path_new + "/AB2_AL_00_tokitiri.pickle",
    pre_path_new + "/AB2_AL_00_tosiwiho.pickle",
    pre_path_new + "/AB2_AL_00_vakavugo.pickle",
    pre_path_new + "/AB2_AL_00_vebegeve.pickle",
    pre_path_new + "/AB2_AL_00_vodosupi.pickle",
    pre_path_new + "/AB2_AL_00_wulesifo.pickle",

    pre_path_new + "/AB2_AL_01_daburalo.pickle",
    pre_path_new + "/AB2_AL_01_dataritu.pickle",
    pre_path_new + "/AB2_AL_01_debatihe.pickle",
    pre_path_new + "/AB2_AL_01_dodukugu.pickle",
    pre_path_new + "/AB2_AL_01_dorimopi.pickle",
    pre_path_new + "/AB2_AL_01_dosesabi.pickle",
    pre_path_new + "/AB2_AL_01_fogedopo.pickle",
    pre_path_new + "/AB2_AL_01_gebohohe.pickle",
    pre_path_new + "/AB2_AL_01_lubulosi.pickle",
    pre_path_new + "/AB2_AL_01_masomapi.pickle",
    pre_path_new + "/AB2_AL_01_nagoduro.pickle",
    pre_path_new + "/AB2_AL_01_nohobune.pickle",
    pre_path_new + "/AB2_AL_01_nositofo.pickle",
    pre_path_new + "/AB2_AL_01_nunasihi.pickle",
    pre_path_new + "/AB2_AL_01_nupidolo.pickle",
    pre_path_new + "/AB2_AL_01_petugagu.pickle",
    pre_path_new + "/AB2_AL_01_rofudopi.pickle",
    pre_path_new + "/AB2_AL_01_sobigodo.pickle",
    pre_path_new + "/AB2_AL_01_tasaperu.pickle",
    pre_path_new + "/AB2_AL_01_vukutala.pickle",
    pre_path_new + "/AB2_AL_02_bakovovo.pickle",
    pre_path_new + "/AB2_AL_02_dasahedi.pickle",
    pre_path_new + "/AB2_AL_02_davoseba.pickle",
    pre_path_new + "/AB2_AL_02_deragubo.pickle",
    pre_path_new + "/AB2_AL_02_dilekoma.pickle",
    pre_path_new + "/AB2_AL_02_fovobewo.pickle",
    pre_path_new + "/AB2_AL_02_gagotuve.pickle",
    pre_path_new + "/AB2_AL_02_gewikule.pickle",
    pre_path_new + "/AB2_AL_02_giduvotu.pickle",
    pre_path_new + "/AB2_AL_02_gopudavu.pickle",
    pre_path_new + "/AB2_AL_02_hagaseha.pickle",
    pre_path_new + "/AB2_AL_02_hanihusa.pickle",
    pre_path_new + "/AB2_AL_02_hekehoti.pickle",
    pre_path_new + "/AB2_AL_02_hitutihi.pickle",
    pre_path_new + "/AB2_AL_02_holipafe.pickle",
    pre_path_new + "/AB2_AL_02_honovino.pickle",
    pre_path_new + "/AB2_AL_02_kavugedu.pickle",
    pre_path_new + "/AB2_AL_02_kenupuwu.pickle",
    pre_path_new + "/AB2_AL_02_kiwetori.pickle",
    pre_path_new + "/AB2_AL_02_komekima.pickle",
    pre_path_new + "/AB2_AL_02_laribosi.pickle",
    pre_path_new + "/AB2_AL_02_lawufita.pickle",
    pre_path_new + "/AB2_AL_02_lekuguku.pickle",
    pre_path_new + "/AB2_AL_02_luwodihi.pickle",
    pre_path_new + "/AB2_AL_02_movalagu.pickle",
    pre_path_new + "/AB2_AL_02_nafisuta.pickle",
    pre_path_new + "/AB2_AL_02_nafumoda.pickle",
    pre_path_new + "/AB2_AL_02_nifirifi.pickle",
    pre_path_new + "/AB2_AL_02_nigovuku.pickle",
    pre_path_new + "/AB2_AL_02_nusebipi.pickle",
    pre_path_new + "/AB2_AL_02_nutohuvi.pickle",
    pre_path_new + "/AB2_AL_02_pamivili.pickle",
    pre_path_new + "/AB2_AL_02_pedoratu.pickle",
    pre_path_new + "/AB2_AL_02_rodedami.pickle",
    pre_path_new + "/AB2_AL_02_sofaruhi.pickle",
    pre_path_new + "/AB2_AL_02_tofafanu.pickle",
    pre_path_new + "/AB2_AL_02_verepepo.pickle",
    pre_path_new + "/AB2_AL_02_vupodiba.pickle",
    pre_path_new + "/AB2_AL_02_wiberola.pickle",
    pre_path_new + "/AB2_AL_02_wuheruge.pickle",
    pre_path_new + "/AB2_AL_03_bokuvoki.pickle",
    pre_path_new + "/AB2_AL_03_didutoga.pickle",
    pre_path_new + "/AB2_AL_03_duvisuku.pickle",
    pre_path_new + "/AB2_AL_03_ganigitu.pickle",
    pre_path_new + "/AB2_AL_03_kakimuni.pickle",
    pre_path_new + "/AB2_AL_03_lamebume.pickle",
    pre_path_new + "/AB2_AL_03_lelaniso.pickle",
    pre_path_new + "/AB2_AL_03_lomiwini.pickle",
    pre_path_new + "/AB2_AL_03_miherepi.pickle",
    pre_path_new + "/AB2_AL_03_nipisali.pickle",
    pre_path_new + "/AB2_AL_03_norafida.pickle",
    pre_path_new + "/AB2_AL_03_pivihuvo.pickle",
    pre_path_new + "/AB2_AL_03_rireseka.pickle",
    pre_path_new + "/AB2_AL_03_riresori.pickle",
    pre_path_new + "/AB2_AL_03_setorike.pickle",
    pre_path_new + "/AB2_AL_03_sipomude.pickle",
    pre_path_new + "/AB2_AL_03_sudipaku.pickle",
    pre_path_new + "/AB2_AL_03_vevutigi.pickle",
    pre_path_new + "/AB2_AL_03_vonisusi.pickle",
    pre_path_new + "/AB2_AL_03_vivorolu.pickle",

    ]

#__|

#| - GP-UCB Aquisition
pre_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB2/gp_ucb_True",
    )
pre_path_new = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/01_abx_al_runs_new/out_data/AB2/gp_ucb_True"
    )


files_list_ab2_gp_ucb = [
    pre_path_new + "/AB2_AL_00_binavepo.pickle",
    pre_path_new + "/AB2_AL_00_duvewire.pickle",
    pre_path_new + "/AB2_AL_00_gihamoba.pickle",
    pre_path_new + "/AB2_AL_00_gusahose.pickle",
    pre_path_new + "/AB2_AL_00_hidawone.pickle",
    pre_path_new + "/AB2_AL_00_ketegopi.pickle",
    pre_path_new + "/AB2_AL_00_kupekomu.pickle",
    pre_path_new + "/AB2_AL_00_lafunida.pickle",
    pre_path_new + "/AB2_AL_00_lalefilu.pickle",
    pre_path_new + "/AB2_AL_00_leriheso.pickle",
    pre_path_new + "/AB2_AL_00_lorawovu.pickle",
    pre_path_new + "/AB2_AL_00_nehuderu.pickle",
    pre_path_new + "/AB2_AL_00_nerakipe.pickle",
    pre_path_new + "/AB2_AL_00_nusepuba.pickle",
    pre_path_new + "/AB2_AL_00_rakamidi.pickle",
    pre_path_new + "/AB2_AL_00_tifosede.pickle",
    pre_path_new + "/AB2_AL_00_vanetesu.pickle",
    pre_path_new + "/AB2_AL_00_vevovofo.pickle",
    pre_path_new + "/AB2_AL_00_wanufome.pickle",
    pre_path_new + "/AB2_AL_00_wihelawi.pickle",

    pre_path_new + "/AB2_AL_01_futakita.pickle",
    pre_path_new + "/AB2_AL_01_gomohuhi.pickle",
    pre_path_new + "/AB2_AL_01_hibevewe.pickle",
    pre_path_new + "/AB2_AL_01_lavuwili.pickle",
    pre_path_new + "/AB2_AL_01_lelefubi.pickle",
    pre_path_new + "/AB2_AL_01_lutodudi.pickle",
    pre_path_new + "/AB2_AL_01_malasoma.pickle",
    pre_path_new + "/AB2_AL_01_manefeno.pickle",
    pre_path_new + "/AB2_AL_01_nimowagi.pickle",
    pre_path_new + "/AB2_AL_01_pobavavu.pickle",
    pre_path_new + "/AB2_AL_01_podegopi.pickle",
    pre_path_new + "/AB2_AL_01_ruboreno.pickle",
    pre_path_new + "/AB2_AL_01_rugatumu.pickle",
    pre_path_new + "/AB2_AL_01_sanovawa.pickle",
    pre_path_new + "/AB2_AL_01_sibubilu.pickle",
    pre_path_new + "/AB2_AL_01_tasorutu.pickle",
    pre_path_new + "/AB2_AL_01_tetoware.pickle",
    pre_path_new + "/AB2_AL_01_tofolami.pickle",
    pre_path_new + "/AB2_AL_01_vunokiti.pickle",
    pre_path_new + "/AB2_AL_01_wepevese.pickle",
    pre_path_new + "/AB2_AL_02_baweteho.pickle",
    pre_path_new + "/AB2_AL_02_dalafoda.pickle",
    pre_path_new + "/AB2_AL_02_fudunihi.pickle",
    pre_path_new + "/AB2_AL_02_gekipafo.pickle",
    pre_path_new + "/AB2_AL_02_hahapohu.pickle",
    pre_path_new + "/AB2_AL_02_honitavu.pickle",
    pre_path_new + "/AB2_AL_02_hukugite.pickle",
    pre_path_new + "/AB2_AL_02_kokakeri.pickle",
    pre_path_new + "/AB2_AL_02_kolinele.pickle",
    pre_path_new + "/AB2_AL_02_latelosi.pickle",
    pre_path_new + "/AB2_AL_02_loguraka.pickle",
    pre_path_new + "/AB2_AL_02_mesapovo.pickle",
    pre_path_new + "/AB2_AL_02_misasoro.pickle",
    pre_path_new + "/AB2_AL_02_newuveha.pickle",
    pre_path_new + "/AB2_AL_02_pirovibu.pickle",
    pre_path_new + "/AB2_AL_02_rawosiga.pickle",
    pre_path_new + "/AB2_AL_02_riwosewo.pickle",
    pre_path_new + "/AB2_AL_02_ropeseso.pickle",
    pre_path_new + "/AB2_AL_02_silopuvu.pickle",
    pre_path_new + "/AB2_AL_02_wuladugo.pickle",
    pre_path_new + "/AB2_AL_03_bahuruba.pickle",
    pre_path_new + "/AB2_AL_03_bowakatu.pickle",
    pre_path_new + "/AB2_AL_03_bubirevu.pickle",
    pre_path_new + "/AB2_AL_03_dufavage.pickle",
    pre_path_new + "/AB2_AL_03_dupapebi.pickle",
    pre_path_new + "/AB2_AL_03_fodoguno.pickle",
    pre_path_new + "/AB2_AL_03_gafifori.pickle",
    pre_path_new + "/AB2_AL_03_gavoheta.pickle",
    pre_path_new + "/AB2_AL_03_gigomuso.pickle",
    pre_path_new + "/AB2_AL_03_gowivugo.pickle",
    pre_path_new + "/AB2_AL_03_kabefira.pickle",
    pre_path_new + "/AB2_AL_03_kebosime.pickle",
    pre_path_new + "/AB2_AL_03_keparudo.pickle",
    pre_path_new + "/AB2_AL_03_kudagawu.pickle",
    pre_path_new + "/AB2_AL_03_lemumimu.pickle",
    pre_path_new + "/AB2_AL_03_lidelaso.pickle",
    pre_path_new + "/AB2_AL_03_lotuhalu.pickle",
    pre_path_new + "/AB2_AL_03_mafapile.pickle",
    pre_path_new + "/AB2_AL_03_megolawe.pickle",
    pre_path_new + "/AB2_AL_03_minonalo.pickle",
    pre_path_new + "/AB2_AL_03_musihedi.pickle",
    pre_path_new + "/AB2_AL_03_nuwerugo.pickle",
    pre_path_new + "/AB2_AL_03_palomogo.pickle",
    pre_path_new + "/AB2_AL_03_papademi.pickle",
    pre_path_new + "/AB2_AL_03_pevinudu.pickle",
    pre_path_new + "/AB2_AL_03_pikavalu.pickle",
    pre_path_new + "/AB2_AL_03_pisahife.pickle",
    pre_path_new + "/AB2_AL_03_pofutefa.pickle",
    pre_path_new + "/AB2_AL_03_sadiwoni.pickle",
    pre_path_new + "/AB2_AL_03_saduhida.pickle",
    pre_path_new + "/AB2_AL_03_sapunumo.pickle",
    pre_path_new + "/AB2_AL_03_todalaka.pickle",
    pre_path_new + "/AB2_AL_03_vevufulo.pickle",
    pre_path_new + "/AB2_AL_03_vonemire.pickle",
    pre_path_new + "/AB2_AL_03_vovekaba.pickle",
    pre_path_new + "/AB2_AL_03_vulowohe.pickle",
    pre_path_new + "/AB2_AL_03_wemifopi.pickle",
    pre_path_new + "/AB2_AL_03_wenulidu.pickle",
    pre_path_new + "/AB2_AL_03_wibukepa.pickle",
    pre_path_new + "/AB2_AL_03_winegifa.pickle",

    # pre_path_new + "/TEST_AL_6_nihagipo.pickle",
    # pre_path_new + "/TEST_AL_6_radetama.pickle",
    # pre_path_new + "/TEST_AL_6_rimuwumi.pickle",


    ]
#__|

#__|



al_data_files_dict = {

    "AB3": {
        "files_list_gp_ucb": files_list_ab3_gp_ucb,
        "files_list_random": files_list_ab3_random,
        },

    "AB2": {
        "files_list_random": files_list_ab2_random,
        "files_list_gp_ucb": files_list_ab2_gp_ucb,
        },

    }

#__|
