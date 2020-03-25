"""
TEMP

"""

#| - Import Modules
import os
import sys
#__|

#  stoich_i = "AB2"
stoich_i = "AB3"

#| - Top IDs to Track
top_ids_to_track_ab2 = [
    '64cg6j9any',
    'n36axdbw65',
    'clc2b1mavs',
    'ck638t75z3',
    'mkbj6e6e9p',
    'b49kx4c19q',
    '85z4msnl6o',
    'bpc2nk6qz1',
    '926dnunrxf',
    'mwmg9p7s6o',
    ]
#__|

# | - top_ids_to_track_ab3
top_ids_to_track_ab3 = [

    '8p8evt9pcg',
    'zimixdvdxd',
    'xw9y6rbkxr',
    'b5cgvsb16w',
    'mj7wbfb5nt',

    '949rnem5z2',
    '6pvt6r95ve',
    '8l919k6s7p',
    '7ic1vt7pz4',
    'zwvqnhbk7f',

    ]
#__|


#| - Files to read ############################################################

if stoich_i == "AB3":

    #| - Random Aquisition
    pre_path = os.path.join(
        os.environ["PROJ_irox"],
        "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB3/random_True",
        )

    files_list_random = [
        pre_path + "/AL_duhakiro.pickle",
        pre_path + "/AL_firegohi.pickle",
        pre_path + "/AL_givohegu.pickle",
        pre_path + "/AL_tiweluku.pickle",
        pre_path + "/AL_vevenuwa.pickle",
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

    files_list_gp_ucb = [
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

elif stoich_i == "AB2":

    #| - Random Aquisition
    pre_path = os.path.join(
        os.environ["PROJ_irox"],
        "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB2/random_True",
        )

    files_list_random = [
        pre_path + "/AL_fodagafu.pickle",
        pre_path + "/AL_ginetape.pickle",
        pre_path + "/AL_litupita.pickle",
        pre_path + "/AL_povibanu.pickle",
        pre_path + "/AL_towohegu.pickle",
        ]

    #__|

    #| - GP-UCB Aquisition
    pre_path = os.path.join(
        os.environ["PROJ_irox"],
        "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB2/gp_ucb_True",
        )

    files_list_gp_ucb = [
        pre_path + "/AL_bapihato.pickle",
        pre_path + "/AL_bavapive.pickle",
        pre_path + "/AL_piritapo.pickle",
        pre_path + "/AL_ralusubi.pickle",
        pre_path + "/AL_tifotogi.pickle",
        pre_path + "/TEST_AL_leteruma.pickle",
        pre_path + "/TEST_AL_mabulute.pickle",
        pre_path + "/TEST_AL_muvasubu.pickle",
        pre_path + "/TEST_AL_wakuhifa.pickle",
        pre_path + "/TEST_AL_wivuheda.pickle",

        ]
    #__|

#__|





#| - __old__
    # '8wxibl7lm4',
    # 'ngn4xec1mo',
    # 'xlziv2zr9g',
    # 'mp6lno9jzr',
    # 'zgxg9o7kny',

    # 'v2blxebixh',
    # '6384vt7pml',
    # '9i6ixublcr',
    # 'xlbfb49wml',
    # 'cqbrnhbacg',
    #
    # 'bgcpc2vabf',
    # 'nrml6dms9l',
    # '9pb4c1927h',
    # 'xozr8f7p7g',
    # 'vgzkzqx5cf',
    #
    # 'mfme9lbtxq',
    # 'mtv3ca7585',
    # 'mkmsvkcyc5',
    # '7h7yns937p',
    # '8gnt9rnd92',
    #
    # 'xfvhck9gcs',
    # 'mrbine8k72',
    # '6qvlcl6iv2',
    # 'vtnfnrbh6s',
    # 'njntmu9w93',
    #
    # 'xy6kzjninu',
    # 'zimuby8uzj',
    # '6tmjv4myvg',
    # 'cy94mecq6g',
    # '6f9fzy7ac5',
    #
    # 'mjctxrx3zf',
    # 'me8d9sx47e',
    # '826d7a7y6t',
    # 'b1mh6qmj62',
    # 'zk9q9yn3b2',
    #
    # 'ntx38tmfxi',
    # 'n5858d8o9j',
    # '82bucfb28a',
    # 'vwxfn3blxi',
    # '71ndxsch8y',

#__|
