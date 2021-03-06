"""
"""

# | - Import Modules
import os
import sys

import pickle

import pandas as pd
#__|


def get_ml_dataframes(
    names=[
        "bulk_dft_data_path",
        "unique_ids_path",
        "prototypes_data_path",
        "static_irox_structures_path",
        "static_irox_structures_kirsten_path",
        "oqmd_irox_data_path",
        "df_features_pre_opt_path",
        "df_features_pre_opt_kirsten_path",
        "df_features_post_opt_path",
        # "df_features_path",
        # "df_features_cleaned_path",
        # "df_features_cleaned_pca_path",
        "oer_bulk_structures_path",
        "df_ccf_path",
        "df_dij_path",
        "ids_to_discard__too_many_atoms_path",
        "df_prototype_dft",
        "df_prototype_static",
        "df_dft_final_final_path",
        "ids_duplicates_path",
        ],

    ):
    """
    """
    # | - get_ml_dataframes

    # | - Import paths from data file
    sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
    from proj_data_irox import (
        bulk_dft_data_path,
        unique_ids_path,
        prototypes_data_path,
        static_irox_structures_path,
        static_irox_structures_kirsten_path,
        oqmd_irox_data_path,
        fp_base_path,
        df_features_pre_opt_path,
        df_features_pre_opt_kirsten_path,
        df_features_post_opt_path,
        # df_features_path,
        # df_features_cleaned_path,
        # df_features_cleaned_pca_path,
        oer_bulk_structures_path,
        df_ccf_path,
        df_dij_path,
        ids_to_discard__too_many_atoms_path,

        df_prototype_dft_path,
        df_prototype_static_path,
        df_dft_final_final_path,
        ids_duplicates_path,
        )

    data_paths = {
        "bulk_dft_data_path": bulk_dft_data_path,
        "unique_ids_path": unique_ids_path,
        "prototypes_data_path": prototypes_data_path,
        "static_irox_structures_path": static_irox_structures_path,
        "static_irox_structures_kirsten_path": static_irox_structures_kirsten_path,
        "oqmd_irox_data_path": oqmd_irox_data_path,
        "df_features_pre_opt_path": df_features_pre_opt_path,
        "df_features_pre_opt_kirsten_path": df_features_pre_opt_kirsten_path,
        "df_features_post_opt_path": df_features_post_opt_path,
        # "df_features_path": df_features_path,
        # "df_features_cleaned_path": df_features_cleaned_path,
        # "df_features_cleaned_pca_path": df_features_cleaned_pca_path,
        "oer_bulk_structures_path": oer_bulk_structures_path,
        "df_ccf_path": df_ccf_path,
        "df_dij_path": df_dij_path,
        "ids_to_discard__too_many_atoms_path": ids_to_discard__too_many_atoms_path,

        #  "df_prototype_dft": df_prototype_dft,
        #  "df_prototype_static": df_prototype_static,
        "df_prototype_dft_path": df_prototype_dft_path,
        "df_prototype_static_path": df_prototype_static_path,
        "df_dft_final_final_path": df_dft_final_final_path,
        "ids_duplicates_path": ids_duplicates_path,
        }
    #__|

    temp = data_paths
    data_dict = dict()
    for key, path in temp.items():

        #  print("key:", key)
        #  print("path:", path)

        if key not in names:
            continue

        if key[-5:] == "_path":
            key_new = key[:-5]
        else:
            key_new = key

        is_pickle = False
        is_csv = False
        try:
            ext = path.split("/")[-1].split(".")[-1]
            if ext == "pickle":
                is_pickle = True
            if ext == "csv":
                is_csv = True
        except:
            pass

        data_i = None
        if is_pickle:
            try:
                # print("TEMP sidjfijsd9i", path)
                with open(path, "rb") as fle:
                    data_i = pickle.load(fle)
            except:
                data_i = None
        elif is_csv:
            data_i = pd.read_csv(path)
        else:
            pass

        if data_i is None:
            print(key)
            print(ext)

        data_dict[key_new] = data_i

    temp = data_dict
    return(temp)
    #__|


def get_data_for_al(
    stoich="AB2",
    verbose=True,
    drop_too_many_atoms=True,
    ):
    """
    """
    # | - get_data_for_al
    sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
    from proj_data_irox import (ids_to_discard__too_many_atoms_path)

    # | - Get all necessary dfs
    df_dict = get_ml_dataframes(
        names=[
            "bulk_dft_data_path",
            "unique_ids_path",
            # "prototypes_data_path",
            "static_irox_structures_path",
            # "static_irox_structures_kirsten_path",
            # "oqmd_irox_data_path",
            "df_features_pre_opt_path",
            "df_features_pre_opt_kirsten_path",
            "df_features_post_opt_path",
            # "oer_bulk_structures_path",
            # "df_ccf_path",
            "df_dij_path",
            # "ids_to_discard__too_many_atoms_path",
            ],
        )

    df_ids = df_dict.get("unique_ids", None)
    df_bulk_dft = df_dict.get("bulk_dft_data", None)
    df_features_pre = df_dict.get("df_features_pre_opt", None)
    # df_features_pre = df_dict.get("df_features_pre_opt_kirsten", None)
    df_features_post = df_dict.get("df_features_post_opt", None)

    df_dij = df_dict.get("df_dij", None)

    print("ISDFIODISFIDS*F*SDF*SDYUGFSODIUFG")
    print("6fcdbh9fz2 in df_bulk_dft", "6fcdbh9fz2" in df_bulk_dft.index)


    df_static_irox = df_dict.get("static_irox_structures", None)
    #__|

    # | - Filter ids to user specifications
    df_ids = df_ids[
        (df_ids["stoich"] == stoich) & \
        (df_ids["source"] != "oqmd") & \
        (df_ids["source"] != "raul_oer") & \
        [True for i in range(len(df_ids))]]
    ids = df_ids["unique_ids"]
    #__|

    # #####################################################

    # | - DFT dataframe
    df_i = df_bulk_dft

    # print("isidfjisdjifjsidjf8yu2894h90832uy4908tyu98023wht0982quj098gtfujw3e")
    # print(df_i.index.shape)
    # print(df_i.index.unique().shape)

    # Common ids between user ids and df
    common_ids = list(set(df_i.index) & set(ids))

    ids_not_in__df_i = [i for i in ids if i not in common_ids]

    df_i = df_i.loc[common_ids]

    if verbose:
        print("len(ids):", len(ids))
        print("len(common_ids)", len(common_ids))
        print("len(ids_not_in__bulk_dft_data):", len(ids_not_in__df_i))
        print("\n", "df_i.shape: ", df_i.shape, sep="")

    df_i = df_i[df_i.source == "raul"]

    df_bulk_dft = df_i

    print("ISDFIODISFIDS*F*SDF*SDYUGFSODIUFG")
    print("6fcdbh9fz2 in df_bulk_dft", "6fcdbh9fz2" in df_bulk_dft.index)

    # print("TEMP TEMP TEMP 89ihsjdgf", "6dzhcimdxs" in df_bulk_dft.index)
    #__|

    # | - Featurs pre-DFT
    df_i = df_features_pre

    # Common ids between user ids and df
    common_ids = list(set(df_i.index) & set(ids))

    ids_not_in__df_i = [i for i in ids if i not in common_ids]

    df_i = df_i.loc[common_ids]

    if verbose:
        print("len(ids):", len(ids))
        print("len(common_ids)", len(common_ids))
        print("len(ids_not_in__bulk_dft_data):", len(ids_not_in__df_i))
        print("\n", "df_i.shape: ", df_i.shape, sep="")

    df_features_pre = df_i

    df_features_pre = df_features_pre["voronoi"]
    #__|

    # | - Features post-DFT
    df_i = df_features_post

    # Common ids between user ids and df
    common_ids = list(set(df_i.index) & set(ids))

    ids_not_in__df_i = [i for i in ids if i not in common_ids]

    df_i = df_i.loc[common_ids]

    if verbose:
        print("len(ids):", len(ids))
        print("len(common_ids)", len(common_ids))
        print("len(ids_not_in__bulk_dft_data):", len(ids_not_in__df_i))
        print("\n", "df_i.shape: ", df_i.shape, sep="")

    df_features_post = df_i

    # Only use post-DFT features from my data set
    df_features_post = \
        df_features_post[df_features_post["data"]["source"] == "raul"]

    df_features_post = df_features_post["voronoi"]
    #__|


    # | - Dropping certain rows
    all_ids = list(set(
        df_bulk_dft.index.tolist() + \
        df_features_pre.index.tolist() + \
        df_features_post.index.tolist() ))


    ids_to_drop = []

    # #########################################################################
    # path_i = os.path.join(
    #     os.environ["PROJ_irox"],
    #     "workflow/ml_modelling/processing_bulk_dft/static_prototypes_structures/out_data",
    #     "ids_to_discard__proto_dupl.pickle")
    # with open(path_i, "rb") as fle:
    #     ids_to_discard__proto_dupl = pickle.load(fle)
    #     # ids_to_drop.extend(ids_to_discard__proto_dupl)
    # #########################################################################

    if drop_too_many_atoms:
        # #####################################################################
        with open(ids_to_discard__too_many_atoms_path, "rb") as fle:
            ids_to_drop__too_many_atoms = pickle.load(fle)
            ids_to_drop.extend(ids_to_drop__too_many_atoms)

            # ids_to_drop = ids_to_drop__too_many_atoms
            # ids_to_drop = [i for i in ids_to_drop if i in all_ids]
        # #####################################################################


    ids_to_drop = [i for i in ids_to_drop if i in all_ids]

    # print("in ids to drop", "6fcdbh9fz2" in ids_to_drop)

    df_features_pre = df_features_pre.drop(
        labels=ids_to_drop, axis=0)

    df_i = df_features_post
    df_features_post = df_i.loc[
        df_i.index.intersection(
            df_features_pre.index.unique()
            ).unique()
        ]


    # print("TEMP TEMP TEMP 89ihsjdgf", "6dzhcimdxs" in df_bulk_dft.index)

    df_bulk_dft = df_bulk_dft.loc[
        df_bulk_dft.index.intersection(
            df_features_pre.index
            ).unique()
        ]

    print("ISDFIODISFIDS*F*SDF*SDYUGFSODIUFG")
    print("6fcdbh9fz2 in df_bulk_dft", "6fcdbh9fz2" in df_bulk_dft.index)

    # print("TEMP TEMP TEMP 89ihsjdgf", "6dzhcimdxs" in df_bulk_dft.index)

    df_static_irox = df_static_irox.loc[
        df_static_irox.index.intersection(
            df_features_pre.index
            ).unique()
        ]

    ids_static = df_dij.index.intersection(df_static_irox["static_id"])
    ids_completed_post_dft = \
        df_dij.index.intersection(df_features_pre.index)


    # print("TEMP TEMP TEMP", "6fcdbh9fz2" in df_dij.index)

    ids_dij = ids_static.tolist() + ids_completed_post_dft.tolist()
    df_dij = df_dij.loc[ids_dij, ids_dij]

    # print("TEMP TEMP TEMP", "6fcdbh9fz2" in df_dij.index)

    #__|

    out_dict = dict()

    out_dict["df_features_post"] = df_features_post
    out_dict["df_features_pre"] = df_features_pre
    out_dict["df_bulk_dft"] = df_bulk_dft

    # TEMP
    out_dict["df_ids"] = df_ids
    out_dict["df_dij"] = df_dij
    out_dict["df_static_irox"] = df_static_irox


    return(out_dict)
    #__|


def create_mixed_df(
    all_ids,
    computed_ids,
    df_post,
    df_pre,

    verbose=True,
    ):
    """Creates dataframe by pulling rows from the df_pre and df_post dfs.

    Resuling df will have the same rows as the all_ids list

    ids in computed_ids will be attempted to be pulled from df_post, otherwise
    the row will be pulled from df_pre
    """
    # | - create_mixed_df
    # #########################################################################

    # TODO Check that computed_ids are in all_ids
    # TODO Check that computed_ids are in df_post
    # TODO Check that all_ids are all in df_pre

    rows_list = []
    for id_i in all_ids:

        if id_i in computed_ids:
            if id_i in df_post.index:
                row_i = df_post.loc[id_i]
                if verbose:
                    print("Got a row from the post df")
                rows_list.append(row_i)
            else:
                if verbose:
                    print("NO GOOD!!!!!! ijjstfy8wjsifd")
        else:
            if id_i in df_pre.index:
                row_i = df_pre.loc[id_i]
                rows_list.append(row_i)
            else:
                # if verbose:
                print("id not in df_pre, should be this is the fall back df")


    df_out = pd.DataFrame(rows_list)

    return(df_out)
    #__|
