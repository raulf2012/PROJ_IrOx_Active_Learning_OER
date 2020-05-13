
structure_id_map = {

    # #####################################################
    #| - IrO2 Structures
    '64cg6j9any': 'i (rutile)',  # rutile
    'cg8p7fxq65': 'anatase', # anatase
    'm2bs8w82x5': 'brookite', # Brookite
    'n36axdbw65': 'ii', # 2nd stable columbite like?
    '85z4msnl6o': 'iii (pyrite)', # Pyrite
    #'myc4ng73xh': 'v', # Fm3m
    'zizr7rvpxs': 'vi', # Porous
    'b49kx4c19q': 'v (columbite)', # Columbite
    'nscdbpmdct': 'iv',  # P63 (layered)
    #'m2bs8w82x5': 'vi',
    #__|

    # #####################################################
    #| - IrO3
    'mp6lno9jzr': 'i', # 482_2d
    'v2blxebixh': 'ii', # sg=2
    '9i6ixublcr': 'iii', # porous
    'nrml6dms9l': 'iv',   # 472_mplowest _63
    #'xozr8f7p7g': 'iv',  # Mp 2nd sg=38
    '6tmjv4myvg': 'v',  # 1D sg=1
    #__|

    # #####################################################
    #| - Top 8 IrO3 structures
    '8p8evt9pcg': '(1)',
    'zimixdvdxd': '(2)',
    # '6fcdbh9fz2': '(3)',
    'xw9y6rbkxr': '(4)',
    'b5cgvsb16w': '(5)',
    'mj7wbfb5nt': '(6)',
    '949rnem5z2': '(7)',
    '6pvt6r95ve': '(8)',

    #| - __old__
    # #'9lmkmh8s8r': '', # 489_alpha
    # #'zimixdvdxd': '', #492_alpha_like
    # '8p8evt9pcg': '(1)', #'alpha',
    # 'zimixdvdxd': '(2)', #'P6_322',
    # 'b5cgvsb16w': '(3)', #'rutile-like',
    # 'mj7wbfb5nt': '(4)', #'sg=52, battery?',
    # '949rnem5z2': '(5)',   #'sg=53',
    #__|

    #__|

    }


annot_offset_dict = {
    #| - IrO2
    '64cg6j9any': dict( ax=15,   ay=-1,    xanchor="left",  ),  # rutile
    'm2bs8w82x5': dict( ax=15,   ay=-10,   xanchor="left",  ),  # Brookite
    'n36axdbw65': dict( ax=-15,  ay=2,     xanchor="right", ),
    'b49kx4c19q': dict( ax=15,   ay=1,     xanchor="left",  ),  # Columbite
    '85z4msnl6o': dict( ax=20,   ay=-10,   xanchor="left",  ),  # pyrite
    'nscdbpmdct': dict( ax=15,   ay=0,     xanchor="left",  ),
    'cg8p7fxq65': dict( ax=-17,  ay=0,     xanchor="right", ),  # Anatase

    'zizr7rvpxs': dict(ax=15, ay=5),

    # '': dict(ax=15, ay=2, xanchor="left", ),

    #__|

    # | - IrO3
    'mp6lno9jzr': dict(ax=10, ay=10, xanchor="left", ),  # i
    'v2blxebixh': dict(ax=15, ay=0, xanchor="center", ),  # ii
    '9i6ixublcr': dict(ax=0, ay=-15, xanchor="center", ),  # iii
    'nrml6dms9l': dict(ax=-15, ay=0, xanchor="center", ),  # iv
    '6tmjv4myvg': dict(ax=15, ay=10, xanchor="center", ),  # v
    #__|

    #| - Top 8 IrO3 structures
    # '6fcdbh9fz2': dict(ax=10, ay=10),

    '8p8evt9pcg': dict(ax=15,  ay=3,    xanchor="left",  ),  # (1)
    'zimixdvdxd': dict(ax=-7, ay=3,    xanchor="right",  ),  # (2)
    'xw9y6rbkxr': dict(ax=-15, ay=-15,  xanchor="center",  ),  # (4)
    'b5cgvsb16w': dict(ax=0,  ay=-18,    xanchor="center",  ),  # (5)
    'mj7wbfb5nt': dict(ax=15,  ay=10,   xanchor="center",  ),  # (6)
    '949rnem5z2': dict(ax=-5,  ay=-15,  xanchor="center",  ),  # (7)
    '6pvt6r95ve': dict(ax=15,  ay=5,    xanchor="left",  ),  # (8)
    #__|

    }



#| - __old__
# together = True

# See energy_references_dH.py
# Ir_ref = -9.316166236367316
# O_ref = -4.657947279999998

# Ir_ref = -9.32910211636731
# O_ref = -4.64915959

# marker = 'o'
# markersize = 6

#  coord_env_style = {
#      'O:6': "rgb(177,83,101)",
#      'mixed': "rgb(86,193,94)",
#      'T:4': "rgb(192,87,197)",
#      'SS:4': "rgb(133,185,55)",
#      'S:4': "rgb(116,99,206)",
#      'L:2': "rgb(200,187,57)",
#      'S:1': "rgb(100,132,200)",
#      'C:8': "rgb(219,148,53)",
#      'HB:8': "rgb(73,185,209)",
#      'S:5': "rgb(215,83,45)",
#      'T:6': "rgb(87,193,151)",
#      'T:5': "rgb(215,74,145)",
#      'C:12': "rgb(68,136,49)",
#      'A:2':" rgb(213,65,86)",
#      'BS_1:10': "rgb(54,129,91)",
#      'TT_1:9': "rgb(214,144,207)",
#      'TY:3': "rgb(167,156,54)",
#      'TL:3': "rgb(151,81,141)",
#      'SBT:8': "rgb(146,180,102)",
#      'I:12': "rgb(163,84,41)",
#      'SY:4': "rgb(93,112,54)",
#      'BO_2:8': "rgb(225,139,112)",
#      'ST:7': "rgb(106,115,30)",
#      'AC:12': "pink",
#      }
#__|
