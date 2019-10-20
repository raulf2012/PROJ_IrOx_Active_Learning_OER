#!/usr/bin/env python

"""Module to TEMP TEMP.

Author: Raul A. Flores
"""

#| - IMPORT MODULES
import os
import copy

import pickle

import numpy as np
import pandas as pd

# SciKitLearn
from sklearn.decomposition import PCA

# Catlearn
from catlearn.regression.gaussian_process import GaussianProcess
from catlearn.preprocess.clean_data import (
    clean_infinite,
    clean_variance,
    clean_skewness)
from catlearn.preprocess.scaling import standardize
#__|



class ALAnalysis:
    """
    """

    #| - ALAnalysis ***********************************************************
    _TEMP = "TEMP"


    def __init__(self,
        ALBulkOpt=None,
        ):
        """
        """
        #| - __init__

        #| - Setting Argument Instance Attributes
        self.ALBulkOpt = ALBulkOpt
        #__|

        #| - Initializing Internal Instance Attributes

        #__|

        #__|

    #__| **********************************************************************



#  █████  ██      ██████  ██    ██ ██      ██   ██  ██████  ██████  ████████
# ██   ██ ██      ██   ██ ██    ██ ██      ██  ██  ██    ██ ██   ██    ██
# ███████ ██      ██████  ██    ██ ██      █████   ██    ██ ██████     ██
# ██   ██ ██      ██   ██ ██    ██ ██      ██  ██  ██    ██ ██         ██
# ██   ██ ███████ ██████   ██████  ███████ ██   ██  ██████  ██         ██

class ALBulkOpt:
    """
    """

    #| - ALBulkOpt ************************************************************
    _TEMP = "TEMP"


    def __init__(self,
        CandidateSpace=None,
        RegressionModel=None,
        num_seed_calcs=11,
        acquisition_bin=10,
        mode="ATF",  # 'ATF' (after the fact)

        stop_mode=None,  # None, 'num_generations'
        stop_num_generations=None,

        verbose=True,
        name="AL_temp"
        ):
        """Initialize ALBulkOpt class instance.

        Args:
            CandidateSpace: <CandidateSpace>
            RegressionModel: <RegressionModel>
            num_seed_calcs: <int>
                Number of initial DFT calcs to be seeded for the AL loop
            acquisition_bin: <int>
                Number of systems to acquire per generation of AL
            mode: <str>
                Mode of operation for AL campaign
                Currently only 'ATF' (after the fact) is allowed
            stop_mode: None or <str>
                None | Doesn't stop AL until data is exhausted
                num_generations | Stops after certain number of AL generations
            stop_num_generations: <int>
                Number of generations to run before stopping
                Only activated if stop_mode == 'num_generations'
            verbose: <bool>
                Controlls verbosity of modules
            name: <str>
                Name of active learning campaign, will be used for file name(s)
        """
        #| - __init__

        #| - Setting Argument Instance Attributes
        self.CandidateSpace = CandidateSpace
        self.RegressionModel = RegressionModel
        self.num_seed_calcs = num_seed_calcs
        self.acquisition_bin = acquisition_bin
        self.mode = mode
        self.stop_mode = stop_mode
        self.stop_num_generations = stop_num_generations
        self.verbose = verbose
        self.name = name
        #__|

        #| - Initializing Internal Instance Attributes
        self.seed_ids = None
        self.al_gen = 0
        self.al_converged = False
        self.completed_ids = []
        self.al_gen_dict = dict()
        #__|

        # Initialize seed ids
        self.seed_ids = self.get_seed_ids()
        self.completed_ids.extend(self.seed_ids)


        # Checking inputs
        stop_mode = self.stop_mode
        stop_num_generations = self.stop_num_generations

        if stop_mode == "num_generations":
            mess_i = "stop_mode='num_generations', \
                Must pass int to 'stop_num_generations'"
            assert type(stop_num_generations) == type(1), mess_i

        #__|


    def run_AL(self):
        """
        """
        #| - run_AL

        #| - class attributes #################################################
        al_gen = self.al_gen
        verbose = self.verbose
        # seed_ids = self.seed_ids
        acquisition_bin = self.acquisition_bin
        completed_ids = self.completed_ids
        CandidateSpace = self.CandidateSpace
        RegressionModel = self.RegressionModel
        al_gen_dict = self.al_gen_dict

        stop_mode = self.stop_mode
        stop_num_generations = self.stop_num_generations
        # __| #################################################################

        while not self.al_converged:
            print(str(self.al_gen).zfill(3), " | init  | ", 64 * "*")

            #| - Testing whether to break AL loop
            # Calculate remaining ids #########################################
            remaining_ids = \
                CandidateSpace.Y_data.index.difference(completed_ids)

            if len(remaining_ids) == 0:
                print("No remaining ids, breaking loop")
                self.al_converged = True


            # Stop criteria based on stop_num_generations
            if verbose:
                print("self.al_gen:", self.al_gen)
                print("stop_num_generations:", stop_num_generations)

            if self.al_gen >= stop_num_generations - 1:
                self.al_converged = True

                if stop_num_generations == 0:
                    continue
            # __|

            #| - ALGeneration #################################################
            ALGen_i = ALGeneration(
                completed_ids=completed_ids,
                acquisition_bin=acquisition_bin,
                CandidateSpace=CandidateSpace,
                RegressionModel=RegressionModel,
                verbose=verbose)

            al_gen_dict[self.al_gen] = ALGen_i
            self.completed_ids.extend(ALGen_i.new_acquisition)
            # __|

            # Check that the completed_ids are within the available data
            for id_i in completed_ids:
                if id_i not in CandidateSpace.Y_data.index:
                    raise Exception("ISJIFDSJFISDJIfj")

            print(str(self.al_gen).zfill(3), " | final | ", 64 * "*"); print()
            self.al_gen += 1
        # __|

    def get_seed_ids(self):
        """Retrieve the ids for the initial seed calculations."""
        #| - get_seed_ids
        mode = self.mode
        num_seed_calcs = self.num_seed_calcs
        CandidateSpace = self.CandidateSpace

        ids_candidate_space = CandidateSpace.FingerPrints.df_pre.index.tolist()

        # Randomize ids in candidate space
        np.random.shuffle(ids_candidate_space)
        ids_all_randomized = ids_candidate_space
        ids_for_seed = ids_all_randomized

        if mode == "ATF":
            ids_post = CandidateSpace.FingerPrints.df_post.index.tolist()
            ids_in_post_randomized = \
                [i for i in ids_all_randomized if i in ids_post]
            ids_for_seed = ids_in_post_randomized

        seed_ids = ids_for_seed[0:num_seed_calcs]

        return(seed_ids)
        #__|

    def __save_state__(self):
        """
        """
        #| - __save_state__
        AL = self
        name = self.name

        directory = "out_data"
        if not os.path.exists(directory): os.makedirs(directory)
        with open(os.path.join(directory, name + ".pickle"), "wb") as fle:
            pickle.dump(AL, fle)
        # __|

    #__| **********************************************************************



#  █████  ██       ██████  ███████ ███    ██
# ██   ██ ██      ██       ██      ████   ██
# ███████ ██      ██   ███ █████   ██ ██  ██
# ██   ██ ██      ██    ██ ██      ██  ██ ██
# ██   ██ ███████  ██████  ███████ ██   ████

class ALGeneration:
    """
    """

    #| - ALGeneration ************************************************************
    _TEMP = "TEMP"


    def __init__(self,
        completed_ids=None,
        acquisition_bin=None,
        CandidateSpace=None,
        RegressionModel=None,
        verbose=True,
        ):
        """
        """
        #| - __init__

        #| - Setting Argument Instance Attributes
        self.completed_ids = completed_ids
        self.acquisition_bin = acquisition_bin
        self.CandidateSpace = CandidateSpace
        self.RegressionModel = RegressionModel
        self.verbose = verbose
        #__|

        # Copy completed_ids to this instance so it doesn't change
        completed_ids = copy.copy(completed_ids)
        self.completed_ids = completed_ids

        #| - Initializing Internal Instance Attributes
        self.df_train = None
        self.df_test = None

        self.new_acquisition = None
        #__|



        # #####################################################################
        CS = self.CandidateSpace
        FP = CS.FingerPrints
        completed_ids = self.completed_ids


        df_mixed = CS.create_mixed_candidate_space(completed_ids)
        self.df_mixed = df_mixed


        Y_data_series = CS.Y_data_series
        Y_data_series_completed = Y_data_series.loc[completed_ids]

        d = {
            "Y": pd.DataFrame(Y_data_series_completed),
            "X": df_mixed.loc[completed_ids]}
        df_train = pd.concat(d.values(), keys=d.keys(), axis=1, sort=False)

        df_test = df_mixed
        df_test = pd.concat([df_mixed], keys=["X"], axis=1, sort=False)

        FP.clean_data(df_train["X"], df_test["X"])

        FP.pca_analysis()

        df_train = FP.df_train
        df_test = FP.df_test

        self.df_train = df_train
        self.df_test = df_test

        # #####################################################################
        # #####################################################################
        self.run_regression_model()
        self.new_acquisition = self.acquisition(acquisition_method="gp_ucb")

        # #####################################################################
        # Save instances of CandidateSpace and FingerPrints AL_gen instance
        # CS_i = copy.deepcopy(CS)
        # self.CandidateSpace = CS_i
        #
        # FP_i = copy.deepcopy(FP)
        # self.FingerPrints = FP_i
        #__|

    def run_regression_model(self):
        """
        """
        #| - run_regression_model
        # #####################################################################
        RegressionModel = self.RegressionModel
        CandidateSpace = self.CandidateSpace
        df_train = self.df_train
        df_test = self.df_test
        completed_ids = self.completed_ids
        # #####################################################################

        train_x = df_train
        train_targets_all = CandidateSpace.Y_data_series
        train_targets = train_targets_all.loc[completed_ids]

        RM = RegressionModel

        RM.set_df_train(df_train)
        RM.set_train_targets(train_targets)
        RM.set_df_test(df_test)

        RM.run_regression()
        self.RegressionModel = RM


        model = pd.concat([
            CandidateSpace.Y_data_series,
            RM.model,
            ], axis=1, sort=False)

        self.model = model
        # __|

    def acquisition(self,
        acquisition_method="gp_ucb",  # 'gp_ucb' or 'random'
        ):
        """
        """
        #| - acquisition
        # #####################################################################
        acquisition_bin = self.acquisition_bin
        model = self.model
        # #####################################################################
        if acquisition_method == "gp_ucb":
            acquisition_ids_ordered = self.acquisition_gp_ucb(model, kappa=1.)
        elif acquisition_method == "random":
            acquisition_ids_ordered = self.acquisition_random(model)


        # Model ordered based on acquisition function
        model_tmp = model.loc[acquisition_ids_ordered]

        # Remove rows for which DFT data is not available
        model_data_avail = model_tmp[~model_tmp["y_real"].isna()]

        # Remove rows which have already been acquired
        # Only acquire what hasn't already been acquired
        df = model_data_avail
        model_data_avail = df[df["acquired"] == False]

        new_acquis_ids = model_data_avail.index[0:acquisition_bin].tolist()

        return(new_acquis_ids)
        # __|


    def acquisition_gp_ucb(self,
        model_i,
        kappa=1.,
        ):
        """GP-UCB aquisition criteria.

        u = Target - kappa * uncertainty
        """
        #| - acquisition_gp_ucb
        model_acq = copy.deepcopy(model_i)

        pred = model_acq["y"]
        err = model_acq["err"]

        model_acq["acquisition"] = pred - kappa * err

        model_acq = model_acq.sort_values("acquisition", ascending=True)

        self.TEMP = model_acq

        acquisition_index_ordered = model_acq.index

        return(acquisition_index_ordered)


        # | - __old__
        # df_i = df_i[df_i["computed"] == False]

        # if df_bulk_dft_all is not None:
        #     # Doing the aquistition
        #     def method(row_i):
        #         actually_computed = False
        #         if pd.isnull(row_i[y_train_key]): actually_computed = False
        #         else: actually_computed = True
        #         return(actually_computed)
        #
        #     df_i = model_i
        #     df_i["actually_computed"] = df_i.apply(
        #         method,
        #         axis=1)
        #     model_i = df_i


        # df_tmp = df_i[0:aqs_bin_size]
        # df_not_comp_but_needed = df_tmp[df_tmp["actually_computed"] == False]

        # df_tmp = df_i[df_i["actually_computed"] == True][0:aqs_bin_size]
        # new_ids_to_compute = df_tmp.index.tolist()

        # out_dict = {
        #     "new_ids_to_compute": new_ids_to_compute,
        #     "ids_needed_but_not_avail": df_not_comp_but_needed["id"].tolist(),
        #     }

        # return(out_dict)
        # __|

        #__|


    def acquisition_random(
        model_i,
        aqs_bin_size=5,
        df_bulk_dft_all=None,
        y_train_key="energy_pa",
        ):
        """
        """
        #| - job_aquisition
        if df_bulk_dft_all is not None:
            # Doing the aquistition
            def method(row_i):
                actually_computed = False
                if pd.isnull(row_i[y_train_key]): actually_computed = False
                else: actually_computed = True
                return(actually_computed)

            df_i = model_i
            df_i["actually_computed"] = df_i.apply(
                method,
                axis=1)
            model_i = df_i


        # pred = model_i["prediction_unstandardized"]
        # uncert = model_i["uncertainty_unstandardized"]
        # model_i["tmp"] = pred - uncert

        df_i = model_i
        df_i = df_i[df_i["computed"] == False]

        # df_i = df_i[df_i["actually_computed"] == True]
        df_i = df_i.sample(frac=1.)

        # new_ids_to_compute = df_i.index.tolist()
        # new_ids_to_compute = df_i.iloc[0:aqs_bin_size]

        df_tmp = df_i.iloc[0:aqs_bin_size]
        df_not_comp_but_needed = df_tmp[df_tmp["actually_computed"] == False]

        df_2 = df_i[df_i["actually_computed"] == True][0:aqs_bin_size]
        new_ids_to_compute = df_2.index.tolist()

        # df_tmp = df_i[0:aqs_bin_size]
        # df_not_comp_but_needed = df_tmp[df_tmp["actually_computed"] == False]
        # df_tmp = df_i[df_i["actually_computed"] == True][0:aqs_bin_size]
        # new_ids_to_compute = df_tmp.index.tolist()

        out_dict = {
            "new_ids_to_compute": new_ids_to_compute,
            "ids_needed_but_not_avail": df_not_comp_but_needed["id"].tolist(),
            }

        return(out_dict)
        #__|


    #__| **********************************************************************



# ██████  ███████  ██████  ██████  ███████ ███████ ███████
# ██   ██ ██      ██       ██   ██ ██      ██      ██
# ██████  █████   ██   ███ ██████  █████   ███████ ███████
# ██   ██ ██      ██    ██ ██   ██ ██           ██      ██
# ██   ██ ███████  ██████  ██   ██ ███████ ███████ ███████

class RegressionModel:
    """
    """

    #| - RegressionModel ******************************************************
    _TEMP = "TEMP"


    def __init__(self,
        df_train=None,
        train_targets=None,
        df_test=None,

        opt_hyperparameters=True,
        gp_settings_dict=None,
        uncertainty_type="regular",
        verbose=True,
        ):
        """

        Args:
            uncertainty_type: 'regular' or 'with_reg'
                Whether to use the regular uncertainty from the GP or the
                "regulization" corrected one
        """
        #| - __init__

        #| - Setting Argument Instance Attributes
        self.df_train = df_train
        self.train_targets = train_targets
        self.df_test = df_test
        self.opt_hyperparameters = opt_hyperparameters
        self.gp_settings_dict = gp_settings_dict
        self.uncertainty_type = uncertainty_type
        self.verbose = verbose
        #__|

        #| - Initializing Internal Instance Attributes
        self.model = None
        #__|

        #__|

    def set_df_train(self, df_train):
        """
        """
        #| - set_df_train
        self.df_train = df_train
        # __|

    def set_train_targets(self, train_targets):
        """
        """
        #| - set_train_targets
        self.train_targets = train_targets
        # __|

    def set_df_test(self, df_test):
        """
        """
        #| - set_df_test
        self.df_test = df_test
        # __|


    def run_regression(self):
        """
        """
        #| - run_regression
        # #####################################################################
        df_train = self.df_train
        train_targets = self.train_targets
        df_test = self.df_test
        opt_hyperparameters = self.opt_hyperparameters
        gp_settings_dict = self.gp_settings_dict
        uncertainty_type = self.uncertainty_type
        verbose = self.verbose
        # #####################################################################

        train_x = df_train
        train_y = train_targets
        train_y_standard = (train_y - train_y.mean()) / train_y.std()

        gp_model = self.gp_model_catlearn

        gp_model_out, m = gp_model(
            train_x,
            train_y_standard,
            df_predict=df_test,
            opt_hyperparameters=opt_hyperparameters,
            gp_settings_dict=gp_settings_dict)

        # TEMP
        self.gp_model_out = gp_model_out

        if uncertainty_type == "regular":
            gp_model_out = {
                "y": gp_model_out["prediction"],
                "err": gp_model_out["uncertainty"]}
        elif uncertainty_type == "with_reg":
            gp_model_out = {
                "y": gp_model_out["prediction"],
                "err": gp_model_out["uncertainty"]}

        model_0 = pd.DataFrame(
            gp_model_out,
            index=df_test.index)


        #| - Add column to model df that indicates the acquired points
        df_acquired = pd.DataFrame(index=df_train.index.unique())
        df_acquired["acquired"] = True

        model_i = pd.concat(
            [model_0, df_acquired],
            axis=1,
            sort=False)

        model_i = model_i.fillna(value={'acquired': False})
        # __|

        #  ####################################################################
        # Unstandardizing the output ##########################################
        y_std = train_y.std().values[0]
        y_mean = train_y.mean().values[0]

        model_i["y"] = (model_i["y"] * y_std) + y_mean
        model_i["err"] = (model_i["err"] * y_std)

        self.model = model_i
        # __|

    def gp_model_catlearn(self,
        train_features,
        train_target,
        df_predict=None,
        gp_settings_dict={},
        opt_hyperparameters=False,
        ):
        """test_features
        """
        #| - gp_model_catlearn
        test_features = df_predict

        noise_default = 0.01  # Regularisation parameter.
        sigma_l_default = 0.8  # Length scale parameter.
        sigma_f_default = 0.2337970892240513  # Scaling parameter.
        alpha_default = 2.04987167  # Alpha parameter.

        noise = gp_settings_dict.get("noise", noise_default)
        sigma_l = gp_settings_dict.get("sigma_l", sigma_l_default)
        sigma_f = gp_settings_dict.get("sigma_f", sigma_f_default)
        alpha = gp_settings_dict.get("alpha", alpha_default)

        #| - Jose Optimized GP
        # Define initial prediction parameters.
        #
        # noise = 0.0042  # Regularisation parameter.
        # sigma_l = 6.3917  # Length scale parameter.
        # sigma_f = 0.5120  # Scaling parameter.
        # alpha = 0.3907  # Alpha parameter.
        #
        # kdict = [
        #     {
        #         'type': 'quadratic',
        #         'dimension': 'single',
        #         # 'dimension': 'features',
        #         'slope': sigma_l,
        #         'scaling': sigma_f,
        #         'degree': alpha,
        #         }
        #     ]
        #
        # GP = GaussianProcess(
        #     kernel_list=kdict, regularization=noise, train_fp=train_features,
        #     train_target=train_target, optimize_hyperparameters=True,
        #     scale_data=False,
        #     )
        #__|

        #| - Sandbox parameters

        #| - HIDE
        # noise = 0.0042  # Regularisation parameter.
        # sigma_l = 6.3917  # Length scale parameter.
        # sigma_f = 0.5120  # Scaling parameter.
        # alpha = 0.3907  # Alpha parameter.

        # noise = 0.00042  # Regularisation parameter.
        # sigma_l = 3.3917  # Length scale parameter.
        # sigma_f = 1.5120  # Scaling parameter.
        # alpha = 0.1907  # Alpha parameter.

        # noise = 0.01  # Regularisation parameter.
        # sigma_l = 0.8  # Length scale parameter.
        # sigma_f = 0.2337970892240513  # Scaling parameter.
        # alpha = 2.04987167  # Alpha parameter.
        #__|

        kdict = [

            #| - Rational Quadratic Kernel
            # {
            #     'type': 'quadratic',
            #     'dimension': 'single',
            #     # 'dimension': 'features',
            #     'slope': sigma_l,
            #     'scaling': sigma_f,
            #     'degree': alpha,
            #     },
            #__|

            #| - Guassian Kernel (RBF)
            {
                'type': 'gaussian',
                'dimension': 'single',
                # 'dimension': 'features',
                'width': sigma_l,
                'scaling': sigma_f,
                'bounds': ((0.0001, 10.),),
                'scaling_bounds': ((0.0001, 10.),),
                # 'scaling_bounds': (0.0001, 100.),
                },

            {
                'type': 'gaussian',
                'dimension': 'single',
                # 'dimension': 'features',
                'width': sigma_l / 10,
                'scaling': sigma_f / 10,
                'bounds': ((0.0001, 10.),),
                'scaling_bounds': ((0.0001, 10.),),
                },
            #__|

            #| - Constant Kernel
            # {
            #     "type": "constant",
            #     # "operation": 0.2,
            #     # "features": ,
            #     "dimension": "single",
            #     # "dimension": "features",
            #     "const": 0.1,
            #     # "bound": ,
            #     },
            #__|

            ]

        GP = GaussianProcess(
            kernel_list=kdict, regularization=noise, train_fp=train_features,
            train_target=train_target, optimize_hyperparameters=False,
            scale_data=False,
            )


        if opt_hyperparameters:
            GP.optimize_hyperparameters(
                global_opt=False,

                algomin='L-BFGS-B',  # The standard one ***************************

                #| - algomin
                # algomin='Nelder-Mead',  # Seems to work well **********************
                # algomin='Newton-CG',  # Doesn't work
                # algomin='BFGS',  # Didn't work
                # algomin='CG',  # Didn't work
                # algomin='dogleg',  # Didn't work
                # algomin='Powell',  # Didn't work
                # algomin='TNC',  # Does work ***************************************
                # algomin='COBYLA',  # Didn't work
                # algomin='SLSQP  # Does work ***************************************
                # algomin='trust-constr',
                # algomin='trust-ncg',  # Didn't work
                # algomin='trust-krylov',  # Didn't work
                # algomin='trust-exact',  # Didn't work
                # algomin='',
                #__|

                eval_jac=False,
                loss_function='lml',
                # loss_function='rmse',
                )
        #__|

        pred = GP.predict(test_fp=test_features, uncertainty=True)

        pred["prediction"] = pred["prediction"].flatten()

        return(pred, GP)

        #__|

    #__| **********************************************************************



# ███████ ███████  █████  ████████ ██    ██ ██████  ███████ ███████
# ██      ██      ██   ██    ██    ██    ██ ██   ██ ██      ██
# █████   █████   ███████    ██    ██    ██ ██████  █████   ███████
# ██      ██      ██   ██    ██    ██    ██ ██   ██ ██           ██
# ██      ███████ ██   ██    ██     ██████  ██   ██ ███████ ███████

class FingerPrints:
    """
    """

    #| - FingerPrints *********************************************************
    _TEMP = "TEMP"

    def __init__(self,
        df_features_pre,
        df_features_post=None,
        pca_mode="num_comp",  # 'num_comp' or 'perc'
        pca_comp=None,
        pca_perc=None,
        verbose=True,
        ):
        """
        """
        #| - __init__

        #| - Setting Argument Instance Attributes
        self.df_pre = df_features_pre
        self.df_post = df_features_post

        self.pca_mode = pca_mode
        self.pca_comp = pca_comp
        self.pca_perc = pca_perc

        self.verbose = verbose
        #__|

        #| - Initializing Internal Instance Attributes
        self.df_test = None
        self.df_train = None
        #__|

        #__|


    def clean_data(self,
        df_train, df_test,

        clean_variance_flag=True,
        clean_skewness_flag=True,
        clean_infinite_flag=True,
        standardize_data_flag=True,
        ):
        """
        """
        #| - clean_data
        mess_i = "Train and test data need to have the same columns"
        assert list(df_test) == list(df_train), mess_i
        labels = list(df_test)

        # Test and train indices
        train_index = df_train.index
        test_index = df_test.index

        # Clean Variance ######################################################
        if clean_variance_flag:
            out_dict_i = self.clean_variance(df_train, df_test, labels)

            df_train = out_dict_i["df_train"]
            df_test = out_dict_i["df_test"]
            labels = out_dict_i["labels"]

        # Clean Skewness ######################################################
        if clean_skewness_flag:
            out_dict_i = self.clean_skewness(df_train, df_test, labels)

            df_train = out_dict_i["df_train"]
            df_test = out_dict_i["df_test"]
            labels = out_dict_i["labels"]

        # Clean Infinite ######################################################
        if clean_infinite_flag:
            out_dict_i = self.clean_infinite(df_train, df_test, labels)

            df_train = out_dict_i["df_train"]
            df_test = out_dict_i["df_test"]
            labels = out_dict_i["labels"]

        # Standardize #########################################################
        if standardize_data_flag:
            out_dict_i = self.clean_standardize(df_train, df_test)

            df_train = out_dict_i["df_train"]
            df_test = out_dict_i["df_test"]

        # Construct DataFrames ################################################
        df_train = pd.DataFrame(df_train, columns=labels, index=train_index)
        df_test = pd.DataFrame(df_test, columns=labels, index=test_index)

        self.df_test = df_test
        self.df_train = df_train
        #__|


    def clean_variance(self, df_train, df_test, labels):
        """
        """
        #| - clean_variance
        verbose = self.verbose
        # #####################################################################

        if verbose:
            print("Cleaning variance:")
            print("train_data.shape:", df_train.shape)

        cleaned_data = clean_variance(
            df_train,
            test=df_test,
            labels=labels,
            )

        df_train = cleaned_data["train"]
        df_test = cleaned_data["test"]
        labels = cleaned_data["labels"]

        if verbose:
            print("df_train.shape:", df_train.shape); print("")

        out_dict = {
            "df_train": df_train,
            "df_test": df_test,
            "labels": labels,
            }

        return(out_dict)
        # __|

    def clean_skewness(self, df_train, df_test, labels):
        """
        """
        #| - clean_skewness
        verbose = self.verbose
        # #####################################################################

        if verbose:
            print("Cleaning skewness:")
            print("train_data.shape:", df_train.shape)

        cleaned_data = clean_skewness(
            df_train,
            test=df_test,
            labels=labels,
            skewness=2.0)

        df_train = cleaned_data["train"]
        df_test = cleaned_data["test"]
        labels = cleaned_data["labels"]

        if verbose:
            print("train_data.shape:", df_train.shape); print("")

        out_dict = {
            "df_train": df_train,
            "df_test": df_test,
            "labels": labels,
            }

        return(out_dict)
        #__|

    def clean_infinite(self, df_train, df_test, labels):
        """
        """
        #| - clean_infinite
        verbose = self.verbose
        # #####################################################################

        if verbose:
            print("Cleaning infinite:")
            print("train_data.shape:", df_train.shape)
        cleaned_data = clean_infinite(
            df_train,
            test=df_test,
            targets=None,
            labels=labels,
            mask=None,
            max_impute_fraction=0,
            strategy='mean')

        df_train = cleaned_data["train"]
        df_test = cleaned_data["test"]
        labels = cleaned_data["labels"]
        if verbose:
            print("train_data.shape:", df_train.shape); print("")

        out_dict = {
            "df_train": df_train,
            "df_test": df_test,
            "labels": labels,
            }

        return(out_dict)
        # __|

    def clean_standardize(self, df_train, df_test):
        """
        """
        #| - clean_standardize
        cleaned_data = standardize(
            df_train,
            test_matrix=df_test,
            mean=None,
            std=None,
            # COMBAK
            local=True)

        df_train = cleaned_data["train"]
        df_test = cleaned_data["test"]

        out_dict = {
            "df_train": df_train,
            "df_test": df_test,
            }

        return(out_dict)
        #__|



    def pca_analysis(self):
        """
        """
        #| - pca_analysis
        df_pre = self.df_pre
        df_post = self.df_post
        verbose = self.verbose
        df_test = self.df_test
        df_train = self.df_train

        pca_mode = self.pca_mode
        pca_comp = self.pca_comp
        pca_perc = self.pca_perc

        # #####################################################################

        df = df_train
        if pca_mode == "num_comp":
            pca = PCA(
                n_components=pca_comp,
                svd_solver="full",
                whiten=False)
        elif pca_mode == "perc":
            pca = PCA(
                n_components=pca_perc,
                svd_solver="full",
                whiten=True)
        else:
            print("ISDJFIESIFJ NO GOODD")

        pca.fit(df)


        #| - Transforming the training data set
        pca_features_cleaned = pca.transform(df)

        num_pca_comp = pca_features_cleaned.shape[-1]

        if verbose:
            print("num_pca_comp: ", num_pca_comp)
            print(df.shape)

        df_pca = pd.DataFrame(
            pca_features_cleaned,
            columns=['PCA%i' % i for i in range(num_pca_comp)],
            index=df.index)

        if verbose:
            print(df_pca.shape)

        df_train = df_pca
        # __|

        #| - Transforming the test data set
        df = df_test
        pca_features_cleaned = pca.transform(df)

        num_pca_comp = pca_features_cleaned.shape[-1]
        # print("num_pca_comp: ", num_pca_comp)

        df_pca = pd.DataFrame(
            pca_features_cleaned,
            columns=['PCA%i' % i for i in range(num_pca_comp)],
            index=df.index)

        df_test = df_pca
        # __|

        # out_dict["pca"] = pca


        self.df_test = df_test
        self.df_train = df_train
        # __|

    #__| **********************************************************************



#  ██████  █████  ███    ██ ██████  ███████ ██████   █████   ██████ ███████
# ██      ██   ██ ████   ██ ██   ██ ██      ██   ██ ██   ██ ██      ██
# ██      ███████ ██ ██  ██ ██   ██ ███████ ██████  ███████ ██      █████
# ██      ██   ██ ██  ██ ██ ██   ██      ██ ██      ██   ██ ██      ██
#  ██████ ██   ██ ██   ████ ██████  ███████ ██      ██   ██  ██████ ███████

class CandidateSpace:
    """
    """

    #| - CandidateSpace *******************************************************
    _TEMP = "TEMP"

    def __init__(self,
        Y_data=None,
        Y_key="energy_pa",
        FingerPrints=None,
        restrict_to_validated_space=False,
        ):
        """

        Args:
            Y_data: <pandas dataframe>
                Dataframe with Y data, indices must match that of FingerPrints
              data.
            Y_key: <string>
                Column name of Y_data df that correponds to the target value
        """
        #| - __init__

        #| - Setting Argument Instance Attributes
        self.Y_data = Y_data
        self.Y_key = Y_key
        self.FingerPrints = FingerPrints

        self.restrict_to_validated_space = restrict_to_validated_space
        #__|

        #| - Initializing Internal Instance Attributes

        #__|

        self.check_input_data()

        # #####################################################################
        Y_data = self.Y_data
        Y_key = self.Y_key

        self.Y_data_series = Y_data[[Y_key]]
        #__|

    def set_fingerprints(self, FingerPrints):
        """Set instance of FingerPrints class."""
        #| - set_fingerprints
        restrict_to_validated_space = self.restrict_to_validated_space
        # #####################################################################

        self.FingerPrints = FingerPrints

        if restrict_to_validated_space:
            tmp = 42

        # __|

    def create_mixed_candidate_space(self,
        computed_ids,
        ):
        """Create mixed candidate space from pre and post dataframes.

        Creates dataframe by pulling rows from the df_pre and df_post dfs.

        Resuling df will have the same rows as the all_ids list

        ids in computed_ids will be attempted to be pulled from df_post, otherwise
        the row will be pulled from df_pre
        """
        #| - create_mixed_candidate_space
        FingerPrints = self.FingerPrints
        df_pre = FingerPrints.df_pre
        df_post = FingerPrints.df_post
        # #####################################################################

        rows_list = []
        for id_i in df_pre.index:
            if id_i in computed_ids and id_i in df_post.index:
                row_i = df_post.loc[id_i]
            else:
                row_i = df_pre.loc[id_i]

            rows_list.append(row_i)

        df_mixed = pd.DataFrame(rows_list)

        return(df_mixed)

        #| - __old__
        # FingerPrints = self.FingerPrints
        # # TODO Check that computed_ids are in all_ids
        # # TODO Check that computed_ids are in df_post
        # # TODO Check that all_ids are all in df_pre
        # rows_list = []
        # for id_i in all_ids:
        #     if id_i in computed_ids:
        #         if id_i in df_post.index:
        #             row_i = df_post.loc[id_i]
        #             if verbose:
        #                 print("Got a row from the post df")
        #             rows_list.append(row_i)
        #         else:
        #             if verbose:
        #                 print("NO GOOD!!!!!! ijjstfy8wjsifd")
        #     else:
        #         if id_i in df_pre.index:
        #             row_i = df_pre.loc[id_i]
        #             rows_list.append(row_i)
        #         else:
        #             # if verbose:
        #             print("id not in df_pre, should be this is the fall back df")
        # df_out = pd.DataFrame(rows_list)
        # return(df_out)
        # __|

        #__|

    def check_input_data(self):
        """
        """
        #| - check_input_data
        Y_data = self.Y_data

        if len(Y_data.index) != len(Y_data.index.unique()):
            print("Y_data has duplicate rows!!!")

        # __|

    #__| **********************************************************************
