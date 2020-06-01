#!/usr/bin/env python

"""Module to TEMP TEMP.
TEMP NEW
Author: Raul A. Flores
"""

# | - IMPORT MODULES
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

from IPython.display import display
#__|

# This is here so that previous import statements work properly
# from active_learning.al_bulkopt import ALBulkOpt



# ██████  ███████  ██████  ██████  ███████ ███████ ███████
# ██   ██ ██      ██       ██   ██ ██      ██      ██
# ██████  █████   ██   ███ ██████  █████   ███████ ███████
# ██   ██ ██      ██    ██ ██   ██ ██           ██      ██
# ██   ██ ███████  ██████  ██   ██ ███████ ███████ ███████

class RegressionModel:
    """
    """

    # | - RegressionModel ******************************************************
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
        # | - __init__

        # | - Setting Argument Instance Attributes
        self.df_train = df_train
        self.train_targets = train_targets
        self.df_test = df_test
        self.opt_hyperparameters = opt_hyperparameters
        self.gp_settings_dict = gp_settings_dict
        self.uncertainty_type = uncertainty_type
        self.verbose = verbose
        #__|

        # | - Initializing Internal Instance Attributes
        self.model = None
        #__|

        #__|

    def set_df_train(self, df_train):
        """
        """
        # | - set_df_train
        self.df_train = df_train
        #__|

    def set_train_targets(self, train_targets):
        """
        """
        # | - set_train_targets
        self.train_targets = train_targets
        #__|

    def set_df_test(self, df_test):
        """
        """
        # | - set_df_test
        self.df_test = df_test
        #__|


    def run_regression(self):
        """
        """
        # | - run_regression
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
        # TEMP
        # print("train_y.describe():", train_y.describe())
        train_y_standard = (train_y - train_y.mean()) / train_y.std()

        # TEST PRINT TEMP
        # print("Ijfefh69w6y7")
        # print("train_y_standard.describe():", train_y_standard.describe())

        gp_model = self.gp_model_catlearn

        gp_model_out, m = gp_model(
            train_x,
            train_y_standard,
            df_predict=df_test,
            opt_hyperparameters=opt_hyperparameters,
            gp_settings_dict=gp_settings_dict)


        # TEMP
        self.gp_model_out = gp_model_out
        self.gp_model = m

        if uncertainty_type == "regular":
            gp_model_out = {
                "y": gp_model_out["prediction"],
                "err": gp_model_out["uncertainty"]}
        elif uncertainty_type == "with_reg":
            gp_model_out = {
                "y": gp_model_out["prediction"],
                "err": gp_model_out["uncertainty_with_reg"]}

        model_0 = pd.DataFrame(
            gp_model_out,
            index=df_test.index)


        # | - Add column to model df that indicates the acquired points
        df_acquired = pd.DataFrame(index=df_train.index.unique())
        df_acquired["acquired"] = True

        model_i = pd.concat(
            [model_0, df_acquired],
            axis=1,
            sort=False)

        model_i = model_i.fillna(value={'acquired': False})
        #__|

        #  ####################################################################
        # Unstandardizing the output ##########################################

        y_std = train_y.std()

        #  print("")
        #  print("siofisjdf890gsddf898s0ddgfs")
        #  print("y_std:", y_std)
        #  print("isinstance(y_std, np.float64):", isinstance(y_std, np.float64))
        #  print("type(y_std) != float:", type(y_std) != float)
        #  print("siofisjdf890gsddf898s0ddgfs")
        #  print("")

        if type(y_std) != float and not isinstance(y_std, np.float64):
            print("This if is True")
            y_std = y_std.values[0]

        y_mean = train_y.mean()
        # if type(y_mean) != float and type(y_mean) != isinstance(y_mean, np.float64):
        if type(y_mean) != float and not isinstance(y_mean, np.float64):
            y_mean = y_mean.values[0]


        # y_std = train_y.std().values[0]
        # y_mean = train_y.mean().values[0]

        model_i["y"] = (model_i["y"] * y_std) + y_mean
        model_i["err"] = (model_i["err"] * y_std)

        self.model = model_i
        #__|

    def gp_model_catlearn(self,
        train_features,
        train_target,
        df_predict=None,
        gp_settings_dict={},
        opt_hyperparameters=False,
        ):
        """test_features
        """
        # | - gp_model_catlearn
        test_features = df_predict

        noise_default = 0.01  # Regularisation parameter.
        sigma_l_default = 0.8  # Length scale parameter.
        sigma_f_default = 0.2337970892240513  # Scaling parameter.
        alpha_default = 2.04987167  # Alpha parameter.

        noise = gp_settings_dict.get("noise", noise_default)
        sigma_l = gp_settings_dict.get("sigma_l", sigma_l_default)
        sigma_f = gp_settings_dict.get("sigma_f", sigma_f_default)
        alpha = gp_settings_dict.get("alpha", alpha_default)

        # | - Jose Optimized GP
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

        # | - Sandbox parameters

        # | - HIDE
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

            # | - Rational Quadratic Kernel
            # {
            #     'type': 'quadratic',
            #     'dimension': 'single',
            #     # 'dimension': 'features',
            #     'slope': sigma_l,
            #     'scaling': sigma_f,
            #     'degree': alpha,
            #     },
            #__|

            # | - Guassian Kernel (RBF)
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

            # | - Constant Kernel
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


        # print("train_features:", train_features.describe())
        # print("train_target:", train_target.describe())


        GP = GaussianProcess(
            kernel_list=kdict, regularization=noise, train_fp=train_features,
            train_target=train_target, optimize_hyperparameters=False,
            scale_data=False,
            )


        if opt_hyperparameters:
            GP.optimize_hyperparameters(
                global_opt=False,

                algomin='L-BFGS-B',  # The standard one ***********************

                # | - algomin
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


        # TEMP
        # print("test_features.describe():", test_features.describe())

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

    # | - FingerPrints *********************************************************
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
        # | - __init__

        # | - Setting Argument Instance Attributes
        self.df_pre = df_features_pre
        self.df_post = df_features_post

        self.pca_mode = pca_mode
        self.pca_comp = pca_comp
        self.pca_perc = pca_perc

        self.verbose = verbose
        #__|

        # | - Initializing Internal Instance Attributes
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
        # | - clean_data
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
        # | - clean_variance
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
        #__|

    def clean_skewness(self, df_train, df_test, labels):
        """
        """
        # | - clean_skewness
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
        # | - clean_infinite
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
        #__|

    def clean_standardize(self, df_train, df_test):
        """
        """
        # | - clean_standardize
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
        # | - pca_analysis
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

        shared_pca_attributes = dict(
            svd_solver="auto",
            whiten=True,
            )

        if pca_mode == "num_comp":

            num_data_points = df.shape[0]
            if num_data_points < pca_comp:
                pca_comp = num_data_points

            pca = PCA(
                n_components=pca_comp,
                **shared_pca_attributes)

        elif pca_mode == "perc":
            pca = PCA(
                n_components=pca_perc,
                **shared_pca_attributes)

        else:
            print("ISDJFIESIFJ NO GOODD")

        pca.fit(df)


        # | - Transforming the training data set
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
        #__|

        # | - Transforming the test data set
        df = df_test
        pca_features_cleaned = pca.transform(df)

        num_pca_comp = pca_features_cleaned.shape[-1]
        # print("num_pca_comp: ", num_pca_comp)

        df_pca = pd.DataFrame(
            pca_features_cleaned,
            columns=['PCA%i' % i for i in range(num_pca_comp)],
            index=df.index)

        df_test = df_pca
        #__|

        # out_dict["pca"] = pca


        # NEW | Setting pca instance as class attribute so that it can be called latter
        self.PCA = pca

        self.df_test = df_test
        self.df_train = df_train
        #__|

    #__| **********************************************************************



#  ██████  █████  ███    ██ ██████  ███████ ██████   █████   ██████ ███████
# ██      ██   ██ ████   ██ ██   ██ ██      ██   ██ ██   ██ ██      ██
# ██      ███████ ██ ██  ██ ██   ██ ███████ ██████  ███████ ██      █████
# ██      ██   ██ ██  ██ ██ ██   ██      ██ ██      ██   ██ ██      ██
#  ██████ ██   ██ ██   ████ ██████  ███████ ██      ██   ██  ██████ ███████

class CandidateSpace:
    """
    """

    # | - CandidateSpace *******************************************************
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
        # | - __init__

        # | - Setting Argument Instance Attributes
        self.Y_data = Y_data
        self.Y_key = Y_key
        self.FingerPrints = FingerPrints

        self.restrict_to_validated_space = restrict_to_validated_space
        #__|

        # | - Initializing Internal Instance Attributes

        #__|

        self.check_input_data()

        # #####################################################################
        Y_data = self.Y_data
        Y_key = self.Y_key

        self.Y_data_series = Y_data[[Y_key]]
        #__|

    def set_fingerprints(self, FingerPrints):
        """Set instance of FingerPrints class."""
        # | - set_fingerprints
        restrict_to_validated_space = self.restrict_to_validated_space
        # #####################################################################

        self.FingerPrints = FingerPrints

        if restrict_to_validated_space:
            tmp = 42
        #__|

    def create_mixed_candidate_space(self,
        computed_ids,
        ):
        """Create mixed candidate space from pre and post dataframes.

        Creates dataframe by pulling rows from the df_pre and df_post dfs.

        Resuling df will have the same rows as the all_ids list

        ids in computed_ids will be attempted to be pulled from df_post,
        otherwise the row will be pulled from df_pre
        """
        # | - create_mixed_candidate_space
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
        #__|

    def check_input_data(self):
        """
        """
        # | - check_input_data
        Y_data = self.Y_data

        if len(Y_data.index) != len(Y_data.index.unique()):
            print("Y_data has duplicate rows!!!")

        #__|

    #__| **********************************************************************
