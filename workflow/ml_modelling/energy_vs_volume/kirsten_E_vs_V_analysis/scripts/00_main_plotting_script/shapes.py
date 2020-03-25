#| - Import Modules
import plotly.graph_objs as go
#__|


def get_plot_shapes(
    df=None,
    inset_range_0_x=None,
    inset_range_1_x=None,
    ):
    """
    """
    #| - get_plot_shapes

    #| - Inset zoom in rectangles
    ab2_min_e = df[df.stoich == "AB2"].dH.min()
    ab3_min_e = df[df.stoich == "AB3"].dH.min()

    # inset_range_0_x = [9.5, 17.]
    inset_range_0_y = [ab2_min_e - 0.02, ab2_min_e + 0.3]

    # inset_range_1_x = [9.5, 17.]
    inset_range_1_y = [ab3_min_e - 0.02, ab3_min_e + 0.12]

    shape_0 = go.layout.Shape(
        type="rect",
        xref="x",
        yref="y",

        x0=inset_range_0_x[0],
        y0=inset_range_0_y[0],
        x1=inset_range_0_x[1],
        y1=inset_range_0_y[1],

        line=dict(
            color="black",
            width=0.5,
            ),
        # fillcolor="LightSkyBlue",
        )

    shape_1 = go.layout.Shape(
        type="rect",
        xref="x2",
        yref="y2",

        x0=inset_range_1_x[0],
        y0=inset_range_1_y[0],
        x1=inset_range_1_x[1],
        y1=inset_range_1_y[1],

        line=dict(
            color="black",
            width=0.5,
            ),
        # fillcolor="LightSkyBlue",
        )


    #| - Lines connecting insets to main plots

    shared_a = dict(
        x1=40.,
        y1=0.52,
        )
    shared_b = dict(
        x1=23.7,
        y1=2.,
        )


    shape_inset_line_ab2_0 = go.layout.Shape(
        type="line",
        xref="x",
        yref="y",
        x0=inset_range_0_x[0],
        y0=inset_range_0_y[1],
        # x1=23.,
        # y1=1.92,
        line=dict(
            color="black",
            width=0.5,
            ),
        **shared_b,
        )
    shape_inset_line_ab2_1 = go.layout.Shape(
        type="line",
        xref="x",
        yref="y",
        x0=inset_range_0_x[1],
        y0=inset_range_0_y[0],
        line=dict(
            color="black",
            width=0.5,
            ),
        **shared_a,
        )


    shape_inset_line_ab3_0 = go.layout.Shape(
        type="line",
        xref="x2",
        yref="y2",
        x0=inset_range_1_x[0],
        y0=inset_range_1_y[1],
        # x1=23.,
        # y1=1.92,
        line=dict(
            color="black",
            width=0.5,
            ),
        **shared_b,
        )
    shape_inset_line_ab3_1 = go.layout.Shape(
        type="line",
        xref="x2",
        yref="y2",
        x0=inset_range_1_x[1],
        y0=inset_range_1_y[0],
        line=dict(
            color="black",
            width=0.5,
            ),
        **shared_a,
        )

    #__|

    #__|

    #| - Metastability horizontal lines
    metastability_limit = 0.2
    # #############################################################################
    shape_inset_metastability_ab2 = go.layout.Shape(
        type="line",
        x0=0,
        y0=ab2_min_e + metastability_limit,
        x1=20,
        y1=ab2_min_e + metastability_limit,

        xref="x3",
        yref="y3",
        line=dict(
            color="grey",
            width=1,
            )
        )

    shape_inset_metastability_ab3 = go.layout.Shape(
        type="line",
        x0=0,
        y0=ab3_min_e + metastability_limit,
        x1=20,
        y1=ab3_min_e + metastability_limit,

        xref="x4",
        yref="y4",
        line=dict(
            color="grey",
            width=1,
            )
        )
    #__|

    # #############################################################################

    shapes_list = [
        shape_0,
        shape_1,
        #shape_inset_metastability_ab2,
        #shape_inset_metastability_ab3,
        shape_inset_line_ab2_0,
        shape_inset_line_ab2_1,
        shape_inset_line_ab3_0,
        shape_inset_line_ab3_1,
        ]

    out_dict = {
        "inset_range_0_y": inset_range_0_y,
        "inset_range_1_y": inset_range_1_y,
        "shapes_list": shapes_list,
        }

    return(out_dict)
    #__|
