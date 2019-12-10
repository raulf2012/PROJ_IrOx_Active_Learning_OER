
"""
"""

#| - Import Modules
import numpy as np
import pandas as pd

import gpflow

from sklearn.decomposition import PCA

from catlearn.regression.gaussian_process import GaussianProcess

from catlearn.preprocess.clean_data import (
    clean_infinite,
    clean_variance,
    clean_skewness)
from catlearn.preprocess.scaling import standardize
#__|


def gp_model_gpflow(
    train_x,
    train_y,
    df_predict=None,
    ):
    """
    input_dim=20
    # variance=
    # lengthscales=
    active_dims
    ARD=True
    """
    #| - gp_model_gpflow
    dim = len(list(train_x))
    # print(dim)
    k = gpflow.kernels.RBF(
        input_dim=int(dim),
        # input_dim=20,
        # variance=0.05,
        # active_dims=,
        ARD=False,
        lengthscales=0.5,
        )

    m = gpflow.models.GPR(
        train_x.values, np.vstack(train_y),
        kern=k)

    # m.likelihood.variance = np.sqrt(0.00003)
    m.likelihood.variance = np.sqrt(0.01)

    gpflow.train.ScipyOptimizer().minimize(m)

    # mean, var = m.predict_y(df_features_pca)
    mean, var = m.predict_y(df_predict)

    gp_out = {
        # "m": m,
        "prediction": mean.flatten(),
        "uncertainty": var.flatten(),
        # "uncertainty_with_reg": pred["uncertainty_with_reg"],
        }

    # df_i = pd.DataFrame(gp_out, index=df_features_pca.index.to_list())

    return(gp_out, m)
    # return(df_i)
    #__|


def gp_model_catlearn(
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
            },

        {
            'type': 'gaussian',
            'dimension': 'single',
            # 'dimension': 'features',
            'width': sigma_l / 10,
            'scaling': sigma_f / 10,
            },
        #__|

        #| - Constant Kernel
        {
            "type": "constant",
            # "operation": 0.2,
            # "features": ,
            "dimension": "single",
            # "dimension": "features",
            "const": 0.1,
            # "bound": ,
            },
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

    #| - __old__

    #| - Clean data
    # finite_numeric_data = clean_infinite(
    #     train_features,
    #     test=df_predict,
    #     targets=train_target)
    # train_features = finite_numeric_data['train']
    # test_features = finite_numeric_data['test']
    # train_target = finite_numeric_data['targets']
    #
    # data = clean_variance(train_features, test_features)
    # data = normalize(data['train'], data['test'])
    # train_features = data['train']
    # test_features = data['test']
    #__|

    # kernel = [{'type': 'gaussian', 'width': 3, 'scaling': 1.}]
    # GP = GaussianProcess(
    #     train_fp=train_features,
    #     train_target=train_target,
    #     kernel_list=kernel,
    #     regularization=3e-2,
    #     optimize_hyperparameters=True,
    #     scale_data=True,  # True is breaking code sometimes
    #     )
    #__|

    #__|



# #############################################################################
# #############################################################################
def gp_workflow(
    df_features_post=None,
    df_test=None,
    df_bulk_dft=None,
    df_bulk_dft_all=None,
    df_ids=None,

    gp_model=None,
    opt_hyperparameters=False,
    gp_params=None,

    y_train_key="energy_pa",

    run_gp=True,

    verbose=True,



    clean_variance_flag=True,
    clean_skewness_flag=True,
    clean_infinite_flag=True,
    standardize_data_flag=True,

    pca_comp=3,
    pca_perc=0.95,
    pca_mode="num_comp",

    index_order_list=None,
    ):
    """
    """
    #| - gp_workflow
    out_dict = {
        "model": None,
        "model_inst": None,
        "train_x": None,
        "train_y": None,
        "train_y_standard": None,
        "test_x": None,
        "pca": None,
        "df_train_pre_pca": None,
        "df_test_pre_pca": None,
        }


    #| - Standardize Data
    # Testing that the columns of the test and training data set have the columns
    # df_post_cols = list(df_features_post["voronoi"])
    df_post_cols = list(df_features_post)
    # df_pre_cols = list(df_test["voronoi"])
    df_pre_cols = list(df_test)

    mess_i = "Train and test data need to have the same columns"
    assert df_post_cols == df_pre_cols, mess_i
    col_labels = df_post_cols
    labels = col_labels

    # Test and train indices
    # train_index = df_features_post["voronoi"].index
    train_index = df_features_post.index
    # test_index = df_test["voronoi"].index
    test_index = df_test.index


    # train_data = df_features_post["voronoi"].values
    train_data = df_features_post.values
    # test_data = df_test["voronoi"].values
    test_data = df_test.values


    # #########################################################################
    #| - Clean Variance #######################################################
    if clean_variance_flag:
        if verbose:
            print("Cleaning variance:")
            print("train_data.shape:", train_data.shape)

        cleaned_data = clean_variance(
            train_data,
            test=test_data,
            labels=labels,
            # mask=None,
            )

        train_data = cleaned_data["train"]
        test_data = cleaned_data["test"]
        labels = cleaned_data["labels"]

        # print("ISJFSDIJFISDJF**$HT(WJWEGHJSF)")
        # print("train_data.shape:", train_data.shape)
        # print("len(train_index):", len(train_index))
        # print("ISJFSDIJFISDJF**$HT(WJWEGHJSF)")

        out_dict["data_clean__variance_train"] = pd.DataFrame(train_data, columns=labels, index=train_index)
        out_dict["data_clean__variance_test"] = pd.DataFrame(test_data, columns=labels, index=test_index)

        if verbose:
            print("train_data.shape:", train_data.shape); print("")
    #__|

    # #########################################################################
    #| - Clean Skewness #######################################################
    if clean_skewness_flag:
        if verbose:
            print("Cleaning skewness:")
            print("train_data.shape:", train_data.shape)
        cleaned_data = clean_skewness(
            train_data,
            test=test_data,
            labels=labels,
            # mask=None,
            # skewness=3.0,
            skewness=2.0,
            )

        train_data = cleaned_data["train"]
        test_data = cleaned_data["test"]
        labels = cleaned_data["labels"]

        out_dict["data_clean__skewness_train"] = pd.DataFrame(train_data, columns=labels, index=train_index)
        out_dict["data_clean__skewness_test"] = pd.DataFrame(test_data, columns=labels, index=test_index)

        if verbose:
            print("train_data.shape:", train_data.shape); print("")
    #__|

    # #########################################################################
    #| - Clean Infinite #######################################################
    if clean_infinite_flag:
        if verbose:
            print("Cleaning infinite:")
            print("train_data.shape:", train_data.shape)
        cleaned_data = clean_infinite(
            train_data,
            test=test_data,
            targets=None,
            labels=labels,
            mask=None,
            max_impute_fraction=0,
            strategy='mean',
            )

        train_data = cleaned_data["train"]
        test_data = cleaned_data["test"]
        labels = cleaned_data["labels"]

        out_dict["data_clean__infinite_train"] = pd.DataFrame(train_data, columns=labels, index=train_index)
        out_dict["data_clean__infinite_test"] = pd.DataFrame(test_data, columns=labels, index=test_index)

        if verbose:
            print("train_data.shape:", train_data.shape); print("")
    #__|

    # #########################################################################
    #| - Standardize ##########################################################
    if standardize_data_flag:
        cleaned_data = standardize(
            train_data,
            test_matrix=test_data,
            mean=None,
            std=None,
            local=True,  # COMBAK
            )

        train_data = cleaned_data["train"]
        test_data = cleaned_data["test"]

        out_dict["data_clean__std_train"] = pd.DataFrame(train_data, columns=labels, index=train_index)
        out_dict["data_clean__std_test"] = pd.DataFrame(test_data, columns=labels, index=test_index)

    #__|


    pd.DataFrame(train_data, columns=labels, index=train_index)


    # #########################################################################
    # Construct DataFrames ####################################################
    df_train = pd.DataFrame(
        train_data,
        columns=labels,
        index=train_index,
        )

    df_test = pd.DataFrame(
        test_data,
        columns=labels,
        index=test_index,
        )


    out_dict["df_train_pre_pca"] = df_train
    out_dict["df_test_pre_pca"] = df_test
    #__|

    #| - PCA Analysis
    ###########################################################################
    df = df_train

    if pca_mode == "num_comp":
        pca = PCA(
            n_components=pca_comp,
            svd_solver="full",
            # whiten=True,
            whiten=False,
            )
    elif pca_mode == "perc":
        pca = PCA(
            n_components=pca_perc,
            svd_solver="full",
            whiten=True,
            )
    else:
        print("ISDJFIESIFJ NO GOODD")


    pca.fit(df)

    # Transforming the training data set
    pca_features_cleaned = pca.transform(df)

    num_pca_comp = pca_features_cleaned.shape[-1]
    if verbose:
        print("num_pca_comp: ", num_pca_comp)

    df_pca = pd.DataFrame(
        pca_features_cleaned,
        columns=['PCA%i' % i for i in range(num_pca_comp)],
        index=df.index)

    df_train = df_pca


    # Transforming the test data set
    df = df_test
    pca_features_cleaned = pca.transform(df)

    num_pca_comp = pca_features_cleaned.shape[-1]
    # print("num_pca_comp: ", num_pca_comp)

    df_pca = pd.DataFrame(
        pca_features_cleaned,
        columns=['PCA%i' % i for i in range(num_pca_comp)],
        index=df.index)

    df_test = df_pca
    out_dict["pca"] = pca
    #__|


    #| - GP Model
    train_x = df_train
    train_y = df_bulk_dft[y_train_key]
    train_y_standard = (train_y - train_y.mean()) / train_y.std()

    # Pickling data ######################################################
    import os; import pickle
    directory = "out_data"
    if not os.path.exists(directory): os.makedirs(directory)
    with open(os.path.join(directory, "train_x_y.pickle"), "wb") as fle:
        pickle.dump({
            "train_x": train_x,
            "train_y": train_y,
            "train_y_standard": train_y_standard,
            },
            fle)

    out_dict["train_x"] = train_x
    out_dict["train_y"] = train_y
    out_dict["train_y_standard"] = train_y_standard
    out_dict["test_x"] = df_test
    #__|

    if run_gp:
        gp_model_out, m = gp_model(
            train_x,
            train_y_standard,
            df_predict=df_test,
            opt_hyperparameters=opt_hyperparameters,
            gp_settings_dict=gp_params,
            )

        #| - GP Model continued
        model_0 = pd.DataFrame(
            gp_model_out,
            index=df_test.index)

        df_tmp = pd.DataFrame(
            index=df_train[~df_train.index.duplicated(keep='first')].index)
        df_tmp["computed"] = True

        df_ids_0 = df_ids.set_index("unique_ids", inplace=False)
        model_i = pd.concat(
            [
                model_0,
                df_ids_0["id"],
                df_tmp,
                ],
            axis=1,
            sort=False)
        model_i = model_i.fillna(value={'computed': False})

        #  ####################################################################
        # Unstandardizing the output ##########################################
        model_i["prediction_unstandardized"] = \
            (model_i["prediction"] * train_y.std()) + train_y.mean()
        model_i["uncertainty_unstandardized"] = \
            (model_i["uncertainty"] * train_y.std())
        model_i["uncertainty_with_reg_unstandardized"] = \
            (model_i["uncertainty_with_reg"] * train_y.std())
        #__|

        # | - Add real energy data to model df
        if df_bulk_dft_all is not None:
            model_i = pd.concat(
                [model_i, df_bulk_dft_all[y_train_key]],
                axis=1, sort=True)
        elif df_bulk_dft is not None:
            model_i = pd.concat(
                [model_i, df_bulk_dft[y_train_key]],
                axis=1, sort=True)
        #__|

        # Sorting model df by energy and adding an x_axis column
        model_i = model_i.sort_values("prediction_unstandardized")
        model_i["x_axis_ind"] = list(range(len(model_i)))

        indices_sorted = model_i.index.sort_values().tolist()
        model_i = model_i.reindex(indices_sorted)

        out_dict["model"] = model_i
        out_dict["model_inst"] = m

    return(out_dict)
    #__|

# #############################################################################
# #############################################################################


def job_aquisition(
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


    pred = model_i["prediction_unstandardized"]
    uncert = model_i["uncertainty_unstandardized"]
    model_i["tmp"] = pred - uncert

    df_i = model_i.sort_values("tmp", ascending=True)
    df_i = df_i[df_i["computed"] == False]


    df_tmp = df_i[0:aqs_bin_size]
    df_not_comp_but_needed = df_tmp[df_tmp["actually_computed"] == False]


    df_tmp = df_i[df_i["actually_computed"] == True][0:aqs_bin_size]
    new_ids_to_compute = df_tmp.index.tolist()

    out_dict = {
        "new_ids_to_compute": new_ids_to_compute,
        "ids_needed_but_not_avail": df_not_comp_but_needed["id"].tolist(),
        }

    return(out_dict)
    #__|


def random_job_aquisition(
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


def test_al_conv(model_i):
    #| - test_al_conv
    al_converged = False

    model_i = model_i.sort_values("prediction_unstandardized")
    lowest_e_row = model_i.iloc[0]
    min_e = lowest_e_row["prediction_unstandardized"]
    min_index = lowest_e_row.name

    df_tmp = model_i[
        (
            model_i["prediction_unstandardized"] -
            model_i["uncertainty_unstandardized"]) < min_e].drop(
                min_index,
                errors="ignore",
                )

    # Systems that have a lower uncertainty tail lower than the lowest energy
    # structure and also are not computed
    df_tmp = df_tmp[df_tmp["computed"] == False]

    if len(df_tmp) == 0:
        al_converged = True

    return(al_converged)
    #__|






#| - __old__





    #| - __old__
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
    #__|


    #| - __old__
    # df_i = model_i.sort_values("prediction_unstandardized", ascending=True)
    # df_next_systems = df_i[df_i["computed"] == False][0:aqs_bin_size]
    # df_not_computed_but_needed = \
    #     df_next_systems[df_next_systems["actually_computed"] == False]
    #
    # new_ids_to_compute = df_next_systems.index.tolist()
    # df_not_computed_but_needed = \
    #     df_next_systems[df_next_systems["actually_computed"] == False]
    # return(df_i)
    # df_next_systems = df_i[df_i["computed"] == False]
    #
    # print(df_next_systems)
    # new_ids_to_compute = df_next_systems[0:aqs_bin_size].index.tolist()
    # computed_ids += df_next_systems.index.tolist()
    # print("Number of new acquisition points: ", len(df_next_systems[df_next_systems["actually_computed"] == True]))
    #__|

#
# def clean_data(df_features):
#     """
#     """
#     #| - clean_data
#     df_features_i = df_features
#     # df_features_i = df_i["features"]
#
#     df_features = df_features_i
#     df_features_cpy = df_features_i
#
#     train_features = df_features_cpy.values
#     train_labels = list(df_features_cpy)
#
#     #| - Clean variance
#     output = clean_variance(
#         train_features,
#         test=None,
#         labels=train_labels,
#         mask=None,
#         )
#     train_features = output["train"]
#     train_labels = output["labels"]
#     #__|
#
#     #| - Clean infinite
#     output = clean_infinite(
#         train_features,
#         test=None,
#         targets=None,
#         labels=train_labels,
#         mask=None,
#         max_impute_fraction=0,
#         strategy='mean',
#         )
#     train_features = output["train"]
#     train_labels = output["labels"]
#     #__|
#
#     #| - Clean skewness
#     # output = clean_skewness(
#     #     train_features,
#     #     test=None,
#     #     labels=train_labels,
#     #     mask=None,
#     #     skewness=3.,
#     #     )
#     # train_features = output["train"]
#     # train_labels = output["labels"]
#     #
#     # column_labels = output["labels"]
#     #__|
#
#     #| - Standardize Data
#     output = standardize(
#         train_features,
#         test_matrix=None,
#         mean=None,
#         std=None,
#         local=True,
#         )
#     #__|
#
#     df_features_cleaned = pd.DataFrame(
#         data=output["train"],
#         # columns=output["labels"],
#         )
#
#     df_features_cleaned = df_features_cleaned.set_index(df_features_i.index)
#     df_features_i = df_features_cleaned
#
#     return(df_features_i)
#
#     #| - __old__
#     # multi_index = pd.MultiIndex.from_tuples(
#     #     [tuple(i) for i in column_labels],
#     #     # names=("tmp1", "tmp2"),
#     #     )
#     # df_features_cleaned.columns = multi_index
#     # df_features_cleaned = df_features_cleaned.set_index(
#     #     df_features.index,
#     #     drop=True, append=False,
#     #     inplace=False, verify_integrity=False)
#     #__|
#
#     #__|

#__|
