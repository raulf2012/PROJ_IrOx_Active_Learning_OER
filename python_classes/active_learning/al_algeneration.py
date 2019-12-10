#!/usr/bin/env python

"""Module to TEMP TEMP.

Author: Raul A. Flores
"""

#| - IMPORT MODULES
import os
import copy

import numpy as np
import pandas as pd
#__|


#  █████  ██       ██████  ███████ ███    ██
# ██   ██ ██      ██       ██      ████   ██
# ███████ ██      ██   ███ █████   ██ ██  ██
# ██   ██ ██      ██    ██ ██      ██  ██ ██
# ██   ██ ███████  ██████  ███████ ██   ████

class ALGeneration:
    """
    """

    #| - ALGeneration *********************************************************
    _TEMP = "TEMP"


    def __init__(self,
        completed_ids=None,
        acquisition_bin=None,
        CandidateSpace=None,
        RegressionModel=None,
        DuplicateFinder=None,
        duplicate_analysis=True,
        index_acq_gen_dict=None,
        prev_acquisition=None,
        prev_duplicate_ids=None,
        duplicate_swap_dict=None,
        al_gen=None,
        acquisition_method="gp_ucb",
        verbose=True,
        ):
        """
        """
        #| - __init__

        #| - Setting Argument Instance Attributes
        # Copy completed_ids to this instance so it doesn't change
        # self.completed_ids = completed_ids
        completed_ids = copy.copy(completed_ids)
        self.completed_ids = completed_ids

        self.acquisition_bin = acquisition_bin
        self.CandidateSpace = CandidateSpace
        self.RegressionModel = RegressionModel
        self.DuplicateFinder = DuplicateFinder
        self.duplicate_analysis = duplicate_analysis
        self.index_acq_gen_dict = index_acq_gen_dict
        self.prev_acquisition = prev_acquisition
        self.prev_duplicate_ids = prev_duplicate_ids
        self.duplicate_swap_dict = duplicate_swap_dict
        self.al_gen = al_gen
        self.acquisition_method = acquisition_method
        self.verbose = verbose
        #__|

        #| - Initializing Internal Instance Attributes
        self.df_train = None
        self.df_test = None
        self.new_acquisition = None
        self.duplicate_swap_lists = []
        self.indices_that_are_duplicates = []
        # self.duplicate_swap_dict = dict()
        #__|

        # #####################################################################
        self.run_regression_model()

        # Acquisition #########################################################
        acquisition_method = self.acquisition_method
        self.new_acquisition = \
            self.acquisition(acquisition_method=acquisition_method)
        #__|

    def run_regression_model(self):
        """
        """
        #| - run_regression_model
        # #####################################################################
        RM = self.RegressionModel
        CandidateSpace = self.CandidateSpace
        completed_ids = self.completed_ids
        DuplicateFinder = self.DuplicateFinder
        get_df_train_test = self.get_df_train_test
        prev_duplicate_ids = self.prev_duplicate_ids
        duplicate_analysis = self.duplicate_analysis
        __run_duplicate_analysis__ = self.__run_duplicate_analysis__
        # #####################################################################

        # df_train, df_test = self.get_df_train_test()
        df_train, df_test = get_df_train_test()

        # train_x = df_train
        train_targets_all = CandidateSpace.Y_data_series
        train_targets = train_targets_all.loc[completed_ids]

        RM.set_df_train(df_train)
        RM.set_train_targets(train_targets)
        RM.set_df_test(df_test)

        RM.run_regression()


        model = pd.concat([
            CandidateSpace.Y_data_series,
            RM.model,
            ], axis=1, sort=False)

        self.model = model

        if DuplicateFinder is not None and duplicate_analysis:
            __run_duplicate_analysis__(prev_duplicate_ids=prev_duplicate_ids)
        else:
            tmp = 42

        duplicates = self.indices_that_are_duplicates

        indices = model.index.tolist()
        model["duplicate"] = \
            [True if i in duplicates else False for i in indices]

        self.model = model
        #__|

    def get_df_train_test(self):
        """
        """
        #| - get_df_train
        CS = self.CandidateSpace
        FP = CS.FingerPrints
        completed_ids = self.completed_ids
        # #####################################################################

        df_mixed = CS.create_mixed_candidate_space(completed_ids)

        Y_data_series = CS.Y_data_series
        Y_data_series_completed = Y_data_series.loc[completed_ids]



        # Pickling data ######################################################
        import os; import pickle
        # directory = "out_data"
        # if not os.path.exists(directory): os.makedirs(directory)
        # with open(os.path.join(os.environ["HOME"], "__temp__", "TEMP.pickle"), "wb") as fle:
        #     pickle.dump((Y_data_series_completed, df_mixed, completed_ids), fle)
        # #####################################################################



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

        return(df_train, df_test)
        #__|

    def __run_duplicate_analysis__(self,
        # duplicate_ids_prev=None,
        prev_duplicate_ids=None,
        ):
        """

        Args:
            duplicate_ids_prev: <None> or <list of indices>
                If provided, removes previously identified duplicate ids from
                consideration. Once a system has been identified as a duplicate
                then it will permenately stay that way
        """
        #| - __run_duplicate_analysis__
        # #####################################################################
        acquisition_bin = self.acquisition_bin
        model = self.model
        DuplicateFinder = self.DuplicateFinder
        index_acq_gen_dict = self.index_acq_gen_dict
        al_gen = self.al_gen
        # #####################################################################


        # #####################################################################
        #| - Apply 'gen_acquired' to model df
        def method(row_i, index_acq_gen_dict):
            index_i = row_i.name
            gen_i = index_acq_gen_dict.get(index_i, np.nan)
            return(gen_i)

        model["gen_acquired"] = model.apply(
            method, axis=1,
            args=(index_acq_gen_dict, ))
        #__|

        # #####################################################################
        #| - Preparing 'simil_dict_master'
        # Only consider duplicates in the set of structures that have been computed
        model_acq = model[model["acquired"] == True]
        filter_ids = model_acq.index.tolist()

        if prev_duplicate_ids is not None:
            filter_ids = [i for i in filter_ids if i not in prev_duplicate_ids]
        else:
            prev_duplicate_ids = []


        simil_dict_master = dict()
        for index_i in filter_ids:
            simil_dict = DuplicateFinder.i_all_similar(
                index_i, filter_ids=filter_ids)

            simil_dict_master[index_i] = simil_dict

        keys_to_delete = []
        for key, val in simil_dict_master.items():
            if val == dict() or val is None:
                keys_to_delete.append(key)

        for key in keys_to_delete:
            del simil_dict_master[key]
        #__|


        if len(simil_dict_master.keys()) == 0:
            self.indices_that_are_duplicates = prev_duplicate_ids

        else:
            keys = list(simil_dict_master.keys())

            tmp_list = \
                [np.array(list(i.keys())) for i in simil_dict_master.values()]
            all_ids_from_duplicate_analysis = \
                keys + list(np.hstack(tmp_list))
            all_ids_from_duplicate_analysis = \
                list(set(all_ids_from_duplicate_analysis))

            # #################################################################
            # Tracks ids that have already been identified as duplicates
            # Don't consider further, already being removed/treated

            # self.TEMP__df_tmp = df_tmp

            # TEMP
            TEMP_df_tmp_dict = dict()
            indices_that_are_duplicates = []
            duplicate_swap_lists = []
            for key, val in simil_dict_master.items():

                if key in indices_that_are_duplicates:
                    continue

                ids_of_duplicates = [key] + list(val.keys())
                ids_of_duplicates = \
                    [i for i in ids_of_duplicates if i not in indices_that_are_duplicates]

                # Skip loop if no duplicate ids are present
                if len(ids_of_duplicates) <= 1:
                    continue


                # ████████ ███████ ███    ███ ██████
                #    ██    ██      ████  ████ ██   ██
                #    ██    █████   ██ ████ ██ ██████
                #    ██    ██      ██  ██  ██ ██
                #    ██    ███████ ██      ██ ██
                # TEMP
                bool_0 = "momkzj9r84" in ids_of_duplicates
                bool_1 = "9yz2mt8hbh" in ids_of_duplicates
                TEMP_PRINT = False
                if bool_0 and bool_1:
                    TEMP_PRINT = True



                df_tmp = model.loc[ids_of_duplicates]
                df_tmp = df_tmp.sort_values("gen_acquired")
                TEMP_df_tmp_dict[key] = df_tmp

                #| - MISC checks
                assert df_tmp.shape[0] > 1, "Only one row in df_tmp"

                #__|



                earliest_acq_row = df_tmp.iloc[0]
                earlist_gen = earliest_acq_row["gen_acquired"]

                generations_acquired = df_tmp["gen_acquired"].tolist()


                #| - Figuring out what kind of df_tmp we have
                if len(list(set(generations_acquired))) == 1:
                    """
                    All duplicates were acquired from the same generation
                    Hopefully fromthe current one

                    Action: simply take the lowest energy system, no need to
                    replace an already existing system
                    """
                    tmp = 42
                    # display(df_tmp)
                    # print("All duplicates acquired at the same gen | OK")

                elif len(list(set(generations_acquired))) == 2:
                    """
                    Duplicates span 2 generations

                    Action:
                        * Check if lowest energy system lies in the earlier
                        generation, or the later one
                        * If the early gen wins out, then nothing special
                        necessary
                    """
                    df_early = df_tmp[df_tmp["gen_acquired"] == earlist_gen]
                    selected_row_early_gen = \
                        df_early.sort_values("y_real").iloc[0]
                    gen_acq_early = selected_row_early_gen["gen_acquired"]

                    most_stable_row = df_tmp.sort_values("y_real").iloc[0]
                    gen_acq_most_stable = most_stable_row["gen_acquired"]

                    if gen_acq_most_stable > gen_acq_early:

                        tmp_list = [
                            selected_row_early_gen.name,
                            most_stable_row.name,
                            ]

                        duplicate_swap_lists.append(tmp_list)

                elif len(list(set(generations_acquired))) > 2:
                    """
                    Duplicates span more than 2 generations

                    I suspect that this occurs when the following applies
                        A ~= B
                        B ~= C
                        but A != C

                    Action: Not sure, plan out
                    """
                    # df_early = df_tmp[df_tmp["gen_acquired"] == earlist_gen]
                    # selected_row_early_gen = \
                    #     df_early.sort_values("y_real").iloc[0]
                    # gen_acq_early = selected_row_early_gen["gen_acquired"]
                    #
                    # most_stable_row = df_tmp.sort_values("y_real").iloc[0]
                    # gen_acq_most_stable = most_stable_row["gen_acquired"]
                    #
                    # if gen_acq_most_stable > gen_acq_early:
                    #
                    #     tmp_list = [
                    #         selected_row_early_gen.name,
                    #         most_stable_row.name,
                    #         ]
                    #
                    #     duplicate_swap_lists.append(tmp_list)

                    tmp = 42
                #__|

                # IMPORTANT <--------------------------------------------------
                #| - Select desired systems from duplicates ###################
                # Select entry from array of duplicates that will be kept
                # Can either keep the one that was already there (earliest)
                # or we can always select the most stable one

                # selected_row = \
                selected_row_early_gen = \
                    df_tmp[df_tmp["gen_acquired"] == earlist_gen].sort_values("y_real").iloc[0]

                selected_row = df_tmp.sort_values("y_real").iloc[0]
                #__|


                # Adding the ids of all systems in df_tmp (other than the
                # 'selected' one), to the list of duplicate systems

                duplicate_indices_i = df_tmp.index.tolist()
                duplicate_indices_i.remove(selected_row.name)
                indices_that_are_duplicates.extend(duplicate_indices_i)







                # TEMP
                if TEMP_PRINT:
                    TEMP_DATA_DICT = {
                        "df_tmp": df_tmp,
                        "selected_row_early_gen": selected_row_early_gen,
                        "selected_row": selected_row,
                        }
                    self.TEMP_DATA_DICT = TEMP_DATA_DICT
                    print("TEMP_DATA_DICT being set")


            # TEMP
            self.TEMP_df_tmp_dict = TEMP_df_tmp_dict
            self.duplicate_swap_lists = duplicate_swap_lists

            indices_that_are_duplicates = list(set(indices_that_are_duplicates))

            self.indices_that_are_duplicates = \
                indices_that_are_duplicates + prev_duplicate_ids
        #__|

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
        #__|

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

        acquisition_index_ordered = model_acq.index

        return(acquisition_index_ordered)
        #__|

    def acquisition_random(self, model_i):
        """
        """
        #| - job_aquisition
        model_acq = copy.deepcopy(model_i)
        index_list = model_acq.index.tolist()

        np.random.shuffle(index_list)
        acquisition_index_out = index_list

        return(acquisition_index_out)
        #__|

    #__| **********************************************************************
