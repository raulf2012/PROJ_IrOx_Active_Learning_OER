#!/usr/bin/env python

"""Preprocess data from various sources.

Author: Raul A. Flores
"""

#| - Import Modules
import os
import sys

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "data",
        ),
    )

import pickle
import copy
import pandas as pd
pd.options.mode.chained_assignment = None

import numpy as np

from ase_modules.ase_methods import max_force
from dft_job_automat.job_analysis import DFT_Jobs_Analysis
from oxr_reaction.oxr_methods import df_calc_adsorption_e
from energetics.dft_energy import Element_Refs

from ase_modules.ase_methods import create_species_element_dict
from vasp.vasp_methods import parse_incar

from proj_data_irox import (
    h2_ref,
    h2o_ref,
    # h2o_corr,
    # h2_corr,
    )

from proj_data_irox import (
    corrections_dict,
    IrO2_bulk_e_dft,
    IrO3_bulk_e_dft,
    IrO3_rutile_like_bulk_e_dft,
    IrO3_battery_bulk_e_dft,
    )
#__|

#| - Elemental References

# For now I'm applying the entire free energy correction to the electronic
# adsorption energies, so keep the energies of gas references as only
# electronic for now
h2o_ref = h2o_ref
h2_ref = h2_ref

# h2o_ref = h2o_ref + h2o_corr
# h2_ref = h2_ref + h2_corr

Elem_Refs = Element_Refs(
    H2O_dict={
        "gibbs_e": h2o_ref,
        "electronic_e": h2o_ref,
        },

    H2_dict={
        "gibbs_e": h2_ref,
        "electronic_e": h2_ref,
        },
    )

oxy_ref, hyd_ref = Elem_Refs.calc_ref_energies()

oxy_ref = oxy_ref.gibbs_e
hyd_ref = hyd_ref.gibbs_e

# print(20 * "TEMP TEMP TEMP")
# print(oxy_ref); print(hyd_ref)
#__|


def load_df(
    from_file=False,
    root_dir=".",
    data_dir=".",
    file_name="df_master.pickle",
    process_df=True,
    filter_early_revisions=True,
    unique_params=["facet", "coverage_type", "bulk_system", "surface_type"],
    name_list=["facet", "coverage_type", "bulk_system"],
    ):
    """Load dataframe and perform some preprocessing.

    unique_params=["facet", "coverage_type", "bulk_system"],

    Usage:

    df_pourbaix, df_ads, df_surf = load_df(
        from_file=False,
        root_dir=data_dir,
        data_dir=data_dir + "/190103_new_job_df",
        file_name="df_master.pickle",
        process_df=True,
        )

    Args:
        from_file:
    """
    #| - load_df
    if from_file:

        #| - From Saved Pickle File
        print("Attempting to load df from pickle")
        with open(data_dir + "/" + file_name, "rb") as fle:
            df_master = pickle.load(fle)
        return(df_master)
        #__|

    else:

        #| - Process Data Frame
        Jobs = DFT_Jobs_Analysis(
            update_job_state=False,
            job_type_class=None,
            load_dataframe=True,
            root_dir=root_dir,
            working_dir=root_dir,
            dataframe_dir=data_dir,
            )

        if filter_early_revisions:
            df_master = Jobs.filter_early_revisions(Jobs.data_frame)
        else:
            df_master = Jobs.data_frame


        #| - TEMP
        # TEMP_PRINT
        # print("_S)DFJ()SDIFU)_SDF")
        # print(len(df_master))
        #
        # df_master = df_master[
        #     (df_master["dopant"] == "Ni") & \
        #     (df_master["bulk_system"] == "IrO2") & \
        #     (df_master["facet"] == "110") & \
        #     (df_master["site"] == "ir_site")
        #     ]
        #
        # print(len(df_master))
        #__|

        if process_df:

            #| - Short Path
            # Shorter path
            root_dir = "/global/cscratch1/sd/flores12/IrOx_Project"

            def calc_short_path(row, root_dir):
                """Remove root path from string and return short path."""
                short_path = row["path"].replace(root_dir, "")[1:]
                return(short_path)

            df_master["path_short"] = df_master.apply(
                calc_short_path,
                args=(root_dir,),
                axis=1,
                )
            #__|

            #| - Entry Name Column
            # System Specific Name (for legends)
            df_master["name_i"] = df_master["facet"] + \
                " | " + df_master["coverage_type"] + " | " + \
                df_master["bulk_system"]

            df_master["name_i_2"] = df_master["facet"] + \
                "_" + df_master["coverage_type"] + "_" + \
                df_master["bulk_system"]


            # name_list = ["coverage_type", "bulk_system"]
            df_master["name_i_3"] = ""
            for i in name_list:
                df_master["name_i_3"] += df_master[i]
                df_master["name_i_3"] += ", "
            #__|

            #| - NEW | 181226 | Updating "surface_type" column
            # np.nan -> 'NaN'
            if "surface_type" in list(df_master):
                df_master["surface_type"] = df_master["surface_type"].replace(
                    np.nan,
                    "NaN",
                    regex=True,
                    )
            #__|

            #| - Extract Chemical Formula
            def get_chemical_formula(row, index=-1):
                """Extract chemical formula from atoms object.

                Args:
                    row:
                    index:
                """
                chem_form = row["atoms_object"][index].get_chemical_formula()

                return(chem_form)

            # df_master["chem_formula"] = df_master.apply(
            #     get_chemical_formula,
            #     index=-1,
            #     axis=1,
            #     )
            #__|

            #| - Max Force
            def get_max_force(row):
                """
                """
                #| - get_max_force
                if row["atoms_object"] is not None:
                    try:
                        max_force_tmp = max_force(row["atoms_object"][-1])
                        return(max_force_tmp[0])
                    except:
                        return(None)

                else:
                    return(None)
                #__|

            def get_sum_force(row):
                """
                """
                #| - get_sum_force
                if row["atoms_object"] is not None:
                    try:
                        max_force_tmp = max_force(row["atoms_object"][-1])
                        return(max_force_tmp[1])

                    except:
                        return(None)

                else:
                    return(None)
                #__|

            df_master["max_force"] = df_master.apply(
                get_max_force,
                axis=1,
                )

            df_master["sum_force"] = df_master.apply(
                get_sum_force,
                axis=1,
                )
            #__|

            #| - Number of Atoms
            def N_atoms(row):
                """
                """
                #| - N_atoms
                if row["atoms_object"] is not None:

                    try:
                        atoms_i = row.atoms_object[-1]
                        elem_dict = create_species_element_dict(atoms_i)
                        N_atoms = int(sum(list(elem_dict.values())))

                        return(N_atoms)

                    except:

                        try:
                            atoms_i = row.init_atoms
                            elem_dict = create_species_element_dict(atoms_i)
                            N_atoms = int(sum(list(elem_dict.values())))

                            return(N_atoms)

                        except:
                            return(None)

                else:
                    return(None)
                #__|

            df_master["N_atoms"] = df_master.apply(
                N_atoms,
                axis=1,
                )
            #__|

            #| - Element Number Dict
            def get_elem_num_dict(row):
                """
                """
                #| - get_elem_num_dict
                if row["atoms_object"] is not None:

                    try:
                        atoms_i = row.atoms_object[-1]
                        elem_dict = create_species_element_dict(atoms_i)
                        # N_atoms = int(sum(list(elem_dict.values())))

                        return(elem_dict)

                    except:

                        try:
                            atoms_i = row.init_atoms
                            elem_dict = create_species_element_dict(atoms_i)
                            # N_atoms = int(sum(list(elem_dict.values())))

                            return(elem_dict)

                        except:
                            return(None)

                else:
                    return(None)
                #__|

            df_master["elem_num_dict"] = df_master.apply(
                get_elem_num_dict,
                axis=1,
                )
            #__|

            #| - NEW | INCAR Processing
            def parse_incar_tmp(row_i):
                #| - parse_incar_tmp

                # print(row_i.incar)
                # print(row_i.path)
                # print("__(*&__----)")

                if type(row_i.incar) is not list:
                    if pd.isna(row_i.incar):
                        return({})

                incar_dict = parse_incar(row_i.incar)

                return(incar_dict)
                #__|

            df_master["incar_parsed"] = df_master.apply(
                parse_incar_tmp,
                axis=1,
                )

            def ldipol(row_i):
                """
                """
                #| - ldipol
                return(row_i.incar_parsed.get("LDIPOL", np.nan))
                #__|

            df_master["dipole_correction"] = df_master.apply(
                ldipol,
                axis=1,
                )

            def ldau(row_i):
                """
                """
                #| - ldau
                return(row_i.incar_parsed.get("LDAU", np.nan))
                #__|

            df_master["u_correction"] = df_master.apply(
                ldau,
                axis=1,
                )

            #__|

            #| - NEW | 181107 | Getting magmoms From Atoms Calc Object

            def get_final_magmoms(row):
                """
                """
                #| - get_magmoms
                atoms_i = row["atoms_object"]
                if atoms_i is None:
                    return(None)


                if len(atoms_i) == 0:
                    return(None)
                elif len(atoms_i) > 0:
                    final_image = atoms_i[-1]
                    magmoms_i = final_image.get_magnetic_moments()

                    return(magmoms_i)
                #__|

            df_master["magmoms"] = df_master.apply(
                get_final_magmoms,
                axis=1,
                )


            def total_magmom(row):
                """
                """
                #| - total_magmom

                magmoms_i = row["magmoms"]
                if magmoms_i is not None:
                    total_magmom_out = np.sum(magmoms_i)
                else:
                    total_magmom_out = None

                return(total_magmom_out)
                #__|

            df_master["total_magmom"] = df_master.apply(
                total_magmom,
                axis=1,
                )

            def abs_magmom(row):
                """
                """
                #| - abs_magmom
                magmoms_i = row["magmoms"]
                if magmoms_i is not None:
                    abs_magmom_out = np.sum(np.abs(magmoms_i))
                else:
                    abs_magmom_out = None

                return(abs_magmom_out)
                #__|

            df_master["abs_magmom"] = df_master.apply(
                abs_magmom,
                axis=1,
                )

            #__|

            df_m = df_master

            # Drop jobs which have manually been marked as "not successful"
            if "success" in df_m.columns:
                df_m = df_m.drop(df_m[df_m["success"] == False].index)

                # This doesn't work for some reason
                # df_m = df_m.drop(df_m[df_m["success"] is False].index)


            df_surf = df_m[df_m["job_type"] == "surface_energy"]
            df_surf_cov = df_m[df_m["job_type"] == "surface_coverage_energy"]

            frames = [df_surf, df_surf_cov]
            df_surf = pd.concat(frames, sort=False)

            # TODO
            # Combine df_surf and df_surf_cov into single DF

            df_pourbaix = df_m[df_m["job_type"] == "surface_coverage"]
            df_ads = df_m[df_m["job_type"] == "ORR_adsorption"]


            #| - DF Adsorption ************************************************

            if not df_ads.empty:

                #| - TEMP TEMP TEMP
                # print(5 * "Fixing Elec E Values | ")
                # print(5 * "Fixing Elec E Values | ")
                # print(5 * "Fixing Elec E Values | ")
                # print(5 * "Fixing Elec E Values | ")
                #
                # df_tmp = df_ads[
                #     (df_ads["bulk_system"] == "IrO2") & \
                #     (df_ads["dopant"] == "Cr") & \
                #     (df_ads["facet"] == "100") & \
                #     (df_ads["site"] == "ir_site")
                #     ]
                #
                #
                # # bare
                # df_tmp_bare = df_tmp[(df_tmp["adsorbate"] == "bare")]
                # ind_bare = df_tmp_bare.index[0]
                #
                # # -425.344083 is the original number
                # # -424.69857736 is with low dipole moment
                # df_ads.at[ind_bare, "elec_energy"] = -424.69857736
                # # df_ads.at[ind_bare, "elec_energy"] = -425.344083
                #
                #
                # # # ooh
                # # df_tmp_ooh = df_tmp[(df_tmp["adsorbate"] == "ooh")]
                # # ind_ooh = df_tmp_bare.index[0]
                # # df_ads.at[ind_ooh, "elec_energy"] = -440.314515  # -440.314515
                # #
                # # # o
                # # df_tmp_o = df_tmp[(df_tmp["adsorbate"] == "o")]
                # # ind_o = df_tmp_bare.index[0]
                # # df_ads.at[ind_o, "elec_energy"] = -440.314515  # -440.314515
                # #
                # # # oh
                # # df_tmp_oh = df_tmp[(df_tmp["adsorbate"] == "oh")]
                # # ind_oh = df_tmp_bare.index[0]
                # # df_ads.at[ind_oh, "elec_energy"] = -435.784688  # -435.784688
                #__|

                #| - Calculate Adsorption Energy
                def calc_ads_e(group):
                    """Calculate species adsorption energy.

                    Args:
                        group
                    """
                    #| - calc_ads_e
                    df_calc_adsorption_e(
                        group,
                        oxy_ref,
                        hyd_ref,
                        group[
                            group["adsorbate"] == "bare"
                            ]["elec_energy"].iloc[0],
                        corrections_mode="corr_dict",
                        #     corrections_column="gibbs_correction",
                        corrections_dict=corrections_dict,
                        )

                    return(group)
                    #__|

                def calc_ads_e_2(group):
                    """Calculate species adsorption energy.

                    Args:
                        group
                    """
                    #| - calc_ads_e
                    row_i = group[group["coverage_type"] == "bare"]
                    bare_e = row_i["elec_energy"].iloc[0]

                    df_calc_adsorption_e(
                        group,
                        oxy_ref,
                        hyd_ref,
                        bare_e,
                        corrections_mode="corr_dict",
                        #     corrections_column="gibbs_correction",
                        # corrections_dict=corrections_dict,
                        )

                    return(group)
                    #__|


                for param_i in unique_params:
                    if param_i not in list(df_ads):
                        unique_params.remove(param_i)

                grouped = df_ads.groupby(unique_params)
                df_ads = grouped.apply(calc_ads_e)


                # NOTE This will not work with my 'df_calc_adsorption_e' method
                # because it relies on there being an 'adsorbate' attribute to
                # count the number of O and H atoms necessary

                # grouped = df_pourbaix.groupby(["facet", "bulk_system"])
                # df_pourbaix = grouped.apply(calc_ads_e_2)
                #__|

                #| - Ordering Columns

                col_order_list = [
                    # Main system variables
                    "bulk_system",
                    "facet",
                    "adsorbate",
                    "coverage_type",
                    "ooh_direction",

                    # Energetics
                    "ads_e",
                    "elec_energy",

                    # Magnetic Moments
                    "magmoms",
                    "total_magmom",
                    "abs_magmom",

                    "path_short",

                    "name_i",

                    "max_force",
                    "sum_force",
                    "elem_num_dict",

                    "incar",
                    "incar_parsed",

                    # Atom properties
                    "init_atoms",
                    "atoms_object",
                    "N_atoms",

                    "dipole_correction",
                    "u_correction",

                    # Low priority
                    "job_type",
                    "max_revision",
                    "revision_number",
                    "success",
                    "coverage",


                    "path",
                    "name_i_2",
                    "name_i_3",


                    # Not needed columns
                    "Job",
                    "layers",
                    ]

                # TEMP | Made df_reorder into a method
                from misc_modules.pandas_methods import reorder_df_columns
                df_ads = reorder_df_columns(col_order_list, df_ads)

                df_ads = df_ads.drop(
                    [
                        "magmoms",
                        "layers",
                        "Job",
                        "coverage",
                        "success",
                        "revision_number",
                        "max_revision",
                        "job_type",
                        "u_correction",
                        "incar",
                        ],
                    axis=1)

                #__|

                #| - Filtering Out Extra Calculations
                """
                Each OER set (bare, O, OH, OOH) should only have 1 calculation
                corresponding to the 4 intermediate structures

                If there are more than 1 calculation available, then the "best"
                one should be selected.

                This is often based on the criteria of energy (most stable)
                """
                from proj_data_irox import groupby_props
                groupby_props = copy.deepcopy(groupby_props)


                df_ads.loc[df_ads["coverage_type"] == "O-4_OH-0", "coverage_type"] = "o_covered"
                df_ads.loc[df_ads["coverage_type"] == "O-2_OH-0", "coverage_type"] = "o_covered_2"
                df_ads.loc[df_ads["coverage_type"] == "O-2_OH-2", "coverage_type"] = "h_covered"


                groupby_props.append("adsorbate")
                grouped = df_ads.groupby(groupby_props)

                ignore_indices = np.array([])
                for i_ind, (name, group) in enumerate(grouped):
                    props_i = dict(zip(groupby_props, list(name)))
                    df_i = group

                    if len(df_i) > 1:
                        print(""); print("_____")
                        print("more than 1 structure here")
                        if props_i["adsorbate"] == "ooh":
                            if "up" in df_i["ooh_direction"].tolist():
                                ignore_indices_i = list(df_i[df_i["ooh_direction"] != "up"].index.values)
                                ignore_indices = np.append(ignore_indices, ignore_indices_i)

                            elif "sideways" in df_i["ooh_direction"].tolist():
                                ignore_indices_i = list(df_i[df_i["ooh_direction"] != "sideways"].index.values)
                                ignore_indices = np.append(ignore_indices, ignore_indices_i)
                            else:
                                tmp = 42

                        elif props_i["adsorbate"] == "bare":
                            df_copy_i = df_i.copy(deep=True)
                            min_e_ind = df_copy_i["elec_energy"].idxmin()

                            ignore_indices_i = df_copy_i.drop(min_e_ind).index.values
                            ignore_indices = np.append(ignore_indices, ignore_indices_i)

                        elif props_i["adsorbate"] == "o":
                            df_copy_i = df_i.copy(deep=True)
                            min_e_ind = df_copy_i["elec_energy"].idxmin()

                            ignore_indices_i = df_copy_i.drop(min_e_ind).index.values
                            ignore_indices = np.append(ignore_indices, ignore_indices_i)

                        elif props_i["adsorbate"] == "oh":
                            df_copy_i = df_i.copy(deep=True)
                            min_e_ind = df_copy_i["elec_energy"].idxmin()

                            ignore_indices_i = df_copy_i.drop(min_e_ind).index.values
                            ignore_indices = np.append(ignore_indices, ignore_indices_i)

                        else:
                            tmp = 42

                df_ads = df_ads.drop(labels=ignore_indices)
                #__|

            #__| **************************************************************


            #| - DF Surface Energies ******************************************

            if not df_surf.empty:

                #| - Slab Area
                def slab_area(row):
                    """
                    """
                    #| - slab_area
                    if row["atoms_object"] is not None:
                        try:
                            atoms_i = row.atoms_object[-1]
                            cell_i = atoms_i.cell

                            cross_prod_i = np.cross(cell_i[0], cell_i[1])
                            area_i = np.linalg.norm(cross_prod_i)
                            return(area_i)

                        except:

                            try:
                                atoms_i = row.init_atoms
                                cell_i = atoms_i.cell

                                cross_prod_i = np.cross(cell_i[0], cell_i[1])
                                area_i = np.linalg.norm(cross_prod_i)
                                return(area_i)

                            except:
                                return(None)

                    else:
                        return(None)
                    #__|

                df_surf["slab_area"] = df_surf.apply(
                    slab_area,
                    axis=1,
                    )
                #__|

                #| - Bulk Energy per Atom (DFT)
                def bulk_elec_e(row):
                    if row["bulk_system"] == "IrO3":
                        return(IrO3_bulk_e_dft)
                    elif row["bulk_system"] == "IrO2":
                        return(IrO2_bulk_e_dft)
                    elif row["bulk_system"] == "IrO3_rutile-like":
                        return(IrO3_rutile_like_bulk_e_dft)
                    elif row["bulk_system"] == "IrO3_battery":
                        return(IrO3_battery_bulk_e_dft)
                    else:
                        raise ValueError(
                            "Didn't assign reference bulk energy",
                            )

                df_surf["bulk_e_per_atom_DFT"] = df_surf.apply(
                    bulk_elec_e,
                    axis=1,
                    )
                #__|

                #| - Surface Energy (Averaged Bulk Reference State)
                # Calculated surface energy using only the bulk phase as a
                # reference frame. This scheme does not treat
                # non-stoicheometric oxygens, instead just having a
                # (Total Energy) / (Total Atoms) term

                def surface_e_ave_bulk_ref(row):
                    surf_e_i = row.elec_energy \
                        - row.N_atoms * row.bulk_e_per_atom_DFT
                    surf_e_i = surf_e_i / row.slab_area
                    return(surf_e_i)

                df_surf["surface_e_ave_bulk_ref"] = df_surf.apply(
                    surface_e_ave_bulk_ref,
                    axis=1,
                    )
                #__|

                #| - Nonstoicheometric Oxygen Count
                def nonstoich_Os(row_i):
                    bulk_i = row_i.bulk_system
                    atoms_i = row_i.init_atoms

                    if bulk_i == "IrO2":
                        O_Ir_ratio = 2
                    elif bulk_i == "IrO3":
                        O_Ir_ratio = 3
                    elif bulk_i == "IrO3_rutile-like":
                        O_Ir_ratio = 3
                    elif bulk_i == "IrO3_battery":
                        O_Ir_ratio = 3
                    else:
                        raise ValueError(
                            "Not expected bulk encountered!",
                            )

                        # print("Not expected bulk encountered!")


                    elems_dict = create_species_element_dict(
                        atoms_i,
                        include_all_elems=False,
                        elems_to_always_include=None,
                        )

                    N_O_stoich = elems_dict["Ir"] * O_Ir_ratio
                    N_O_nonstoich = elems_dict["O"] - N_O_stoich

                    return(N_O_nonstoich)

                df_surf["nonstoich_Os"] = df_surf.apply(nonstoich_Os, axis=1)
                #__|

            #__| **************************************************************


            #| - Saving Dataframe to Pickle
            with open(data_dir + "/" + file_name, "wb") as fle:
                pickle.dump((df_pourbaix, df_ads, df_surf), fle)
            #__|


            #| - Constructing Output
            out_set = ()

            if not df_pourbaix.empty:
                out_set = out_set + (df_pourbaix,)
            else:
                out_set = out_set + (None,)

            if not df_ads.empty:
                out_set = out_set + (df_ads,)
            else:
                out_set = out_set + (None,)

            if not df_surf.empty:
                out_set = out_set + (df_surf,)
            else:
                out_set = out_set + (None,)
            #__|


            return(out_set)

        else:
            return(df_master)
    #__|

    #__|
