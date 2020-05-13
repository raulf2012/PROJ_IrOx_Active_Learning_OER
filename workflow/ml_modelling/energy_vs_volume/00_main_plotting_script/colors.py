import numpy as np

import colormap as cmaps

import colorlover as cl
from IPython.display import HTML



def hex_to_rgb(hex):
    h = hex
    rgb_str = 'rgb(' + ",".join([str(i) for i in tuple(int(h[i:i+2], 16) for i in (0, 2, 4))]) + ")"
    return(rgb_str)

# hex_to_rgb("d93535d1")


def get_color_scale(df=None, dx=None):
    """
    """

    # #############################################################################
    # Using built in colormaps ####################################################
    c = cmaps.Colormap()
    CM_i = c.cmap("viridis")

    color_pall_tmp = CM_i(np.linspace(0, 1., num=200), alpha=0.6)

    color_pall_tmp2 = []
    for i in color_pall_tmp:

        color_i = "rgba(" + \
        str(255 * i[0]) + "," + \
        str(255 * i[1]) + "," + \
        str(255 * i[2]) + "," + \
        str(i[3]) + \
        ")"

        color_pall_tmp2.append(color_i)
    color_pall_0 = color_pall_tmp2


    # #############################################################################
    # #############################################################################
    # #############################################################################

    # Defining custom colormap
    # color_pall = cl.scales['9']['seq']['BuPu'][::-1][2:-1]
    color_pall_1 = [

        # Grey scale gradient
        hex_to_rgb("d5d5d5d1"),  # Lighter
        hex_to_rgb("929292d1"),  # Darker

        # hex_to_rgb("afafafd1"),
        # hex_to_rgb("1f1f1fd1"),

        # #####################################################################
        # hex_to_rgb("b6d935d1"),
        # hex_to_rgb("e07a58d1"),

        # hex_to_rgb("c9c9c9"),
        # hex_to_rgb("6e6e6e"),

        # hex_to_rgb("afafafd1"),
        # hex_to_rgb("1f1f1fd1"),

        # "rgb(175,175,175)",
        # "rgb(31,31,31)",
        ]

    # #


    color_pall = color_pall_1

    # #############################################################################
    # #############################################################################
    # #############################################################################



    color_pall_interp = cl.interp(color_pall, 80) # Map color scale to 500 bins

    color_pall_interp_new = []
    for color_i in color_pall_interp:
        color_i = color_i.replace("%", "")
        color_i = color_i.replace(" ", "")
        color_pall_interp_new.append(color_i)
    color_pall_interp = color_pall_interp_new
    color_pall_interp = cl.to_rgb(color_pall_interp)

    # #############################################################################
    num_bins = len(color_pall_interp)
    d_step = (1 / (num_bins - 1))

    custom_color_pall_final = []
    for i_cnt, color_i in enumerate(color_pall_interp):
        custom_color_pall_final.append(
            [i_cnt * d_step, color_i]
            )

    #print(custom_color_pall_final)

    colorscale_i = custom_color_pall_final




    # #############################################################################
    # #############################################################################
    # #############################################################################


    # COMBAK AB2/3 have the same min/max ave. coord. so this works, but otherwise we would have to do it separately for each stoich (I think)
    stoich = "AB2"
    df_i = df[df.stoich == stoich]

    z_min = df_i.mean_coor.min()
    z_max = df_i.mean_coor.max()

    b = z_min
    m = (z_max - z_min) / 1.

    color_tmp = "grey"


    special_coord_dict = {

        # Blue and green
        4: dict(color=hex_to_rgb("d93535d1")),
        6: dict(color=hex_to_rgb("35c0d9d1")),

        # Blue and green
        # 4: dict(color=hex_to_rgb("3cc9e3ff")),
        # 6: dict(color=hex_to_rgb("a4e168ff")),

        # 4: dict(color=hex_to_rgb("afafafd1")),
        # 6: dict(color=hex_to_rgb("1f1f1fd1")),
        }

    for coord_i, val in special_coord_dict.items():
        spec_color_i = val["color"]

        lower_bound = (coord_i - dx - b) / m
        upper_bound = (coord_i + dx - b) / m

        # #############################################################################
        for i_cnt, color_i in enumerate(colorscale_i):
            scale_pos_i = color_i[0]
            if scale_pos_i > lower_bound:
                lower_bound_ind = i_cnt
                break
        # #############################################################################
        for i_cnt, color_i in enumerate(colorscale_i):
            scale_pos_i = color_i[0]
            if scale_pos_i > upper_bound:
                upper_bound_ind = i_cnt
                break

        del colorscale_i[lower_bound_ind:upper_bound_ind]

        colorscale_i.insert(lower_bound_ind, [(coord_i + dx - b) / m, spec_color_i])
        colorscale_i.insert(lower_bound_ind, [(coord_i - b) / m, spec_color_i])
        colorscale_i.insert(lower_bound_ind, [(coord_i - dx - b) / m, spec_color_i])

    return(colorscale_i)
