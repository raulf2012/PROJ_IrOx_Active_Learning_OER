
# | - Import Modules
import os
import sys

import copy

import time
import pickle

import numpy as np
import pandas as pd

import gpflow

from sklearn.decomposition import PCA

# #############################################################################
from catlearn.regression.gaussian_process import GaussianProcess
from catlearn.preprocess.clean_data import (
    clean_infinite,
    clean_variance,
    clean_skewness)
from catlearn.preprocess.scaling import standardize

# #############################################################################
from active_learning.al_bulkopt import ALBulkOpt
from active_learning.active_learning import (
    RegressionModel,
    FingerPrints,
    CandidateSpace,
    )
from active_learning.al_analysis import ALAnalysis, ALAnimation

# #############################################################################
from IPython.display import display
# -

# from inputs import (
#     stoich_i,
#     verbose,
#     gp_settings,
#     name_i,
#     )
#__|

def run_al_i(
    stoich_i=None,
    verbose=None,
    gp_settings=None,
    name_i=None,
    save_dir_extra=None,
    acquisition_method=None,
    duplicate_analysis=None,
    seed=None,
    ):
    """
    """
    # | - run_al_i
    # # Read Data

    # + {"jupyter": {"source_hidden": true}}
    sys.path.insert(0, os.path.join(
        os.environ["PROJ_irox"],
        "workflow/ml_modelling"))
    from ml_methods import get_data_for_al

    out_dict = get_data_for_al(
        stoich=stoich_i, verbose=False,
        drop_too_many_atoms=True,
    #     drop_too_many_atoms=False,
        )

    df_bulk_dft = out_dict["df_bulk_dft"]
    df_bulk_dft = df_bulk_dft[df_bulk_dft["source"] == "raul"]

    # df_bulk_dft = df_bulk_dft[["atoms", "energy_pa"]]
    df_bulk_dft = df_bulk_dft[["atoms", "dH"]]
    df_bulk_dft.columns.values[1] = "y_real"

    df_features_pre = out_dict["df_features_pre"]
    df_features_post = out_dict["df_features_post"]

    df_ids = out_dict["df_ids"]


    df_static_irox = out_dict["df_static_irox"]
    df_dij = out_dict["df_dij"]
    # -

    # # Filter to candidates w/ DFT energy

    # + {"jupyter": {"source_hidden": true}}
    ids_w_dft = df_bulk_dft.index

    # TEMP | Reduce size of candidate space
    # np.random.seed(8)
    # ids_w_dft = np.sort(np.random.choice(np.sort(ids_w_dft), size=200))
    ids_w_dft = list(set(ids_w_dft))
    # print("ids_w_dft:", ids_w_dft)

    df_bulk_dft = df_bulk_dft.loc[ids_w_dft]

    df_features_pre = df_features_pre.loc[ids_w_dft]
    df_features_post = df_features_post.loc[ids_w_dft]

    # + {"active": ""}
    #
    #
    #
    #
    # -

    # # CCF Class

    # + {"jupyter": {"source_hidden": true}}
    sys.path.insert(0, os.path.join(
        os.environ["PROJ_irox"],
        "python_classes"))
    from ccf_similarity.ccf import CCF

    d_thresh = 0.02
    CCF = CCF(
        df_dij=df_dij,
        d_thresh=d_thresh)

    # + {"jupyter": {"source_hidden": true}}
    RM = RegressionModel(
        opt_hyperparameters=True,
        gp_settings_dict=gp_settings,
        verbose=verbose,
        )

    FP = FingerPrints(
        df_features_pre,
        df_features_post=df_features_post,
        pca_mode="num_comp",  # 'num_comp' or 'perc'
        pca_comp=10,
        pca_perc=None,
        verbose=verbose,
        )

    CS = CandidateSpace(
        Y_data=df_bulk_dft,
        Y_key="y_real",
        FingerPrints=FP,
        )

    # +
    AL = ALBulkOpt(
        CandidateSpace=CS,
        RegressionModel=RM,
        DuplicateFinder=CCF,  # Optional
        duplicate_analysis=duplicate_analysis,
        # num_seed_calcs=11,
        num_seed_calcs=5,
        acquisition_bin=5,
        # stop_mode="num_generations",
        stop_mode=None,
        stop_num_generations=3,
        name=name_i,
        save_dir_extra=save_dir_extra,
        verbose=verbose,
        # acquisition_method="gp_ucb",
        acquisition_method=acquisition_method,
        seed=seed,
        )

    run_al = True
    if run_al:
        AL.run_AL()
        AL.duplicate_system_history_analysis()
        AL.__save_state__()

    #__|
