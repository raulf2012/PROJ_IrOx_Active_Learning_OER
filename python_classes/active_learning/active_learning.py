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

from IPython.display import display
#__|




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
        DuplicateFinder=None,  # Optional
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
            DuplicateFinder: <CCF or other class with required methods>
                Instance that can be used to identify duplicate structures
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
        self.DuplicateFinder = DuplicateFinder
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
        self.index_acq_gen_dict = dict()
        self.duplicate_ids = []
        self.duplicate_swap_dict = dict()
        self.performance__static_winners = dict()
        #__|

        # Initialize seed ids
        self.seed_ids = self.get_seed_ids()
        self.completed_ids.extend(self.seed_ids)

        self.__check_inputs__()
        #__|


    def run_AL(self):
        """
        """
        #| - run_AL

        #| - class attributes #################################################
        al_gen = self.al_gen
        verbose = self.verbose
        seed_ids = self.seed_ids
        acquisition_bin = self.acquisition_bin
        completed_ids = self.completed_ids
        CandidateSpace = self.CandidateSpace
        RegressionModel = self.RegressionModel
        DuplicateFinder = self.DuplicateFinder
        al_gen_dict = self.al_gen_dict
        duplicate_ids = self.duplicate_ids
        duplicate_swap_dict = self.duplicate_swap_dict

        stop_mode = self.stop_mode
        stop_num_generations = self.stop_num_generations

        index_acq_gen_dict = self.index_acq_gen_dict
        # __| #################################################################

        while not self.al_converged:
            print(str(self.al_gen).zfill(3), " | init  | ", 64 * "*")


            if self.al_gen == 0:
                index_acq_gen_dict_i = dict()
                for index_j in seed_ids:
                    index_acq_gen_dict_i[index_j] = int(self.al_gen)
                self.index_acq_gen_dict.update(index_acq_gen_dict_i)

                prev_acquisition = seed_ids
            else:
                prev_acquisition = al_gen_dict[self.al_gen - 1].new_acquisition


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

            if stop_mode == "num_generations":
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
                DuplicateFinder=DuplicateFinder,
                index_acq_gen_dict=index_acq_gen_dict,
                prev_acquisition=prev_acquisition,
                prev_duplicate_ids=self.duplicate_ids,
                duplicate_swap_dict=duplicate_swap_dict,
                al_gen=self.al_gen,
                verbose=verbose)

            al_gen_dict[self.al_gen] = ALGen_i
            self.completed_ids.extend(ALGen_i.new_acquisition)
            self.duplicate_ids.extend(ALGen_i.indices_that_are_duplicates)

            # TEMP
            print("'6r716sxr9t' in self.duplicate_ids | 2",
                "6r716sxr9t" in self.duplicate_ids)


            # ████████ ███████ ███    ███ ██████
            #    ██    ██      ████  ████ ██   ██
            #    ██    █████   ██ ████ ██ ██████
            #    ██    ██      ██  ██  ██ ██
            #    ██    ███████ ██      ██ ██


            self.duplicate_ids = list(set(self.duplicate_ids))

            # TEMP
            print("'6r716sxr9t' in self.duplicate_ids | 2",
                "6r716sxr9t" in self.duplicate_ids)



            self.duplicate_swap_dict[self.al_gen] = \
                ALGen_i.duplicate_swap_lists


            # #################################################################
            # Updating the 'index_acq_gen_dict' attribute
            # self.index_acq_gen_dict
            index_list = []
            index_acq_gen_dict_i = dict()
            for index_j in ALGen_i.new_acquisition:
                index_list.append(index_j)
                # +1 so that the seed calcs count as the 0th acquistiion
                index_acq_gen_dict_i[index_j] = int(self.al_gen + 1)
            self.index_acq_gen_dict.update(index_acq_gen_dict_i)

            mess = "Seems like an id was acquired in more than 1 generation?"
            assert len(index_list) == len(set(index_list)), mess
            # #################################################################
            # __|


            # Check various performance metrics for AL loop
            self.__evaluate_performance__()

            # Check that the completed_ids are within the available data
            for id_i in completed_ids:
                if id_i not in CandidateSpace.Y_data.index:
                    raise Exception("ISJIFDSJFISDJIfj")

            print(str(self.al_gen).zfill(3), " | final | ", 64 * "*"); print()
            self.al_gen += 1
        # __|

    def add_main_Y_to_model(self,
        model,
        plot_dft_instead_of_pred=True,
        prediction_key="y",
        uncertainty_key="err",
        ):
        """Construct main plotting Y column from model df

        Args:
            plot_dft_instead_of_pred: <bool>
                Use DFT energy if already acquired and available instead of
                predicted value from ML model
        """
        #| - add_main_Y_to_model

        def method(row_i):
            #| - method
            computed = row_i["acquired"]
            y_real = row_i["y_real"]

            actually_computed = False
            if not np.isnan(y_real):
                actually_computed = True
            # dft_energy = y_real

            predicted_energy = row_i[prediction_key]
            predicted_uncert = row_i[uncertainty_key]

            # #####################################################################
            new_column_values_dict = {
                "Y_main": None,
                "Y_uncer": None}

            # #####################################################################
            if computed and actually_computed:
                new_column_values_dict["Y_main"] = y_real
                new_column_values_dict["Y_uncer"] = 0.
            else:
                new_column_values_dict["Y_main"] = predicted_energy
                new_column_values_dict["Y_uncer"] = predicted_uncert

            # #####################################################################
            for key, value in new_column_values_dict.items():
                row_i[key] = value
            return(row_i)
            #__|

        if plot_dft_instead_of_pred:
            model = model.apply(method, axis=1)
        else:
            model["Y_main"] = model[prediction_key]
            model["Y_uncer"] = model[uncertainty_key]

        return(model)
        # __|



    def __evaluate_performance__(self,
        types=["static_winners"],
        ):
        """Evaluate performance metrics of AL loop in real-time.

        Args:
            types: <list of str>
                Types of performance anaysis to perform on AL, must have
                corresponding method defined
        """
        #| - __evaluate_performance__

        # #####################################################################
        _evaluate_performance__static_winners = \
            self._evaluate_performance__static_winners
        meth_static_winners = _evaluate_performance__static_winners
        # #####################################################################

        if "static_winners" in types:
            meth_static_winners()

        # __|

    def _evaluate_performance__static_winners(self):
        """Evaluate whether the identity and ordering of the highest
        performing systems are unchanging over several generations of the AL
        loop
        """
        #| - _evaluate_performance__

        #| - class attributes #################################################
        AL = self
        al_gen = self.al_gen
        verbose = self.verbose
        seed_ids = self.seed_ids
        acquisition_bin = self.acquisition_bin
        completed_ids = self.completed_ids
        CandidateSpace = self.CandidateSpace
        RegressionModel = self.RegressionModel
        DuplicateFinder = self.DuplicateFinder
        al_gen_dict = self.al_gen_dict

        stop_mode = self.stop_mode
        stop_num_generations = self.stop_num_generations

        index_acq_gen_dict = self.index_acq_gen_dict
        # __| #################################################################

        # #####################################################################
        mode = "lowest_N"  # 'lowest_N' or 'lowest_perc'

        N_ids = 10
        lowest_perc = 5

        # Number of consecutive generations that the Nth best systems must
        # remain static
        M_gens = 3
        # #####################################################################

        if mode == "lowest_perc":
            num_candidates = CandidateSpace.FingerPrints.df_pre.shape[0]
            N_ids = int(num_candidates * (lowest_perc * 0.01))

        gen_keys = list(AL.al_gen_dict.keys())

        if len(gen_keys) > M_gens:
            latest_M_keys = gen_keys[-(M_gens + 1):]
            last_gen_key = gen_keys[-1]

            al_gen_dict_subset_i = dict(zip(
                latest_M_keys,
                [AL.al_gen_dict.get(i, None) for i in latest_M_keys]))

            indices_list = []
            iterator = enumerate(al_gen_dict_subset_i.items())
            for i_cnt, (gen_i, AL_i) in iterator:
                model_i = AL_i.model

                model_i = AL.add_main_Y_to_model(
                    model_i, plot_dft_instead_of_pred=True)
                model_i = model_i[(model_i["duplicate"] == False)]
                model_i = model_i.sort_values("Y_main")

                indices_i = model_i.index.tolist()

                indices_list.append(indices_i)

                if i_cnt >= M_gens:
                    indices_i = indices_list[i_cnt][0:N_ids]
                    ids_static_list = []
                    for j in range(M_gens):
                        indices_j = indices_list[i_cnt - (j + 1)][0:N_ids]
                        ids_static = indices_j == indices_i
                        ids_static_list.append(ids_static)

                    ids_are_static = all(ids_static_list)

            self.performance__static_winners[last_gen_key] = ids_are_static
        # __|




    def __check_inputs__(self):
        """Check inputs to class."""
        #| - __check_inputs__
        # #####################################################################
        stop_mode = self.stop_mode
        stop_num_generations = self.stop_num_generations
        # #####################################################################

        if stop_mode == "num_generations":
            mess_i = "stop_mode='num_generations', \
                Must pass int to 'stop_num_generations'"
            assert type(stop_num_generations) == type(1), mess_i
        # __|

    def get_seed_ids(self):
        """Retrieve the ids for the initial seed calculations."""
        #| - get_seed_ids
        # #####################################################################
        mode = self.mode
        num_seed_calcs = self.num_seed_calcs
        CandidateSpace = self.CandidateSpace
        # TEMP
        seed = 20191025
        # #####################################################################

        ids_candidate_space = CandidateSpace.FingerPrints.df_pre.index.tolist()
        ids_candidate_space = np.sort(ids_candidate_space)

        # Randomize ids in candidate space
        np.random.seed(seed)
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

    #| - ALGeneration *********************************************************
    _TEMP = "TEMP"


    def __init__(self,
        completed_ids=None,
        acquisition_bin=None,
        CandidateSpace=None,
        RegressionModel=None,
        DuplicateFinder=None,
        index_acq_gen_dict=None,
        prev_acquisition=None,
        prev_duplicate_ids=None,
        duplicate_swap_dict=None,
        al_gen=None,
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
        self.index_acq_gen_dict = index_acq_gen_dict
        self.prev_acquisition = prev_acquisition
        self.prev_duplicate_ids = prev_duplicate_ids
        self.duplicate_swap_dict = duplicate_swap_dict
        self.al_gen = al_gen
        self.verbose = verbose
        #__|

        #| - Initializing Internal Instance Attributes
        self.df_train = None
        self.df_test = None
        self.new_acquisition = None
        self.duplicate_swap_lists = []
        # self.duplicate_swap_dict = dict()
        #__|

        # #####################################################################

        self.run_regression_model()
        self.new_acquisition = self.acquisition(acquisition_method="gp_ucb")
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

        if DuplicateFinder is not None:
            __run_duplicate_analysis__(
                prev_duplicate_ids=prev_duplicate_ids,
                )

            duplicates = self.indices_that_are_duplicates

            indices = model.index.tolist()
            model["duplicate"] = \
                [True if i in duplicates else False for i in indices]

            self.model = model
        # __|

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
        with open(os.path.join(os.environ["HOME"], "__temp__", "TEMP.pickle"), "wb") as fle:
            pickle.dump((Y_data_series_completed, df_mixed, completed_ids), fle)
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
        # __|


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
        # __|

        # #####################################################################
        #| - Preparing 'simil_dict_master'
        # Only consider duplicates in the set of structures that have been computed
        model_acq = model[model["acquired"] == True]
        filter_ids = model_acq.index.tolist()

        if prev_duplicate_ids is not None:
            filter_ids = [i for i in filter_ids if i not in prev_duplicate_ids]
        else:
            prev_duplicate_ids = []


        # TEMP
        print("'6r716sxr9t' in prev_duplicate_ids",
            "6r716sxr9t" in prev_duplicate_ids)


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
        # __|


        if len(simil_dict_master.keys()) == 0:
            self.indices_that_are_duplicates = []
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


                df_tmp = model.loc[ids_of_duplicates]
                df_tmp = df_tmp.sort_values("gen_acquired")
                TEMP_df_tmp_dict[key] = df_tmp

                # TEMP
                print("key:", key)
                print(40 * "-")
                # print("987gsdfg")
                # display(df_tmp)
                # self.TEMP__df_tmp = df_tmp

                #| - MISC checks
                assert df_tmp.shape[0] > 1, "Only one row in df_tmp"

                # __|



                # #############################################################
                # #############################################################
                # #############################################################
                # #############################################################
                # #############################################################

                TEMP_PRINT = False
                if "64cg6j9any" in df_tmp.index.tolist():
                    if "9yz2mt8hbh" in df_tmp.index.tolist():
                        # if al_gen ==
                        TEMP_PRINT = True
                        print(30 * "H")
                        display(df_tmp)
                        print("al_gen:", al_gen)
                        print(2 * "\n")

                if al_gen == 1 and "9yz2mt8hbh" in df_tmp.index.tolist():
                    TEMP_PRINT = True
                    print(30 * "H")
                    display(df_tmp)
                    print("al_gen:", al_gen)
                    # print("prev_duplicate_ids:", prev_duplicate_ids)
                    print(2 * "\n")

                if "6r716sxr9t" in df_tmp.index.tolist() and al_gen > 1:
                    display(df_tmp)


                # #############################################################
                # #############################################################
                # #############################################################
                # #############################################################
                # #############################################################



                earliest_acq_row = df_tmp.iloc[0]
                earlist_gen = earliest_acq_row["gen_acquired"]

                generations_acquired = df_tmp["gen_acquired"].tolist()

                if TEMP_PRINT:
                    print(len(list(set(generations_acquired))))





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
                # __|



                #| - __old__
                if len(list(set(generations_acquired))) == 1:
                    tmp = 42
                    # display(df_tmp)
                    # print("All duplicates acquired at the same gen | OK")
                else:
                    mess = "There shouldn't be more than one duplicate from previous generations"
                    num_early_gens = generations_acquired.count(
                        earliest_acq_row["gen_acquired"])
                    if num_early_gens > 1:
                        print(key, mess, key)

                    # assert num_early_gens == 1, mess

                # # Are there multiple early gen rows to choose from?
                # # Should only happen if multiple are acquired at once
                # multiple_early_gens_present = False
                # if len(list(set(generations_acquired))) == 1:
                #     print("multiple_early_gens_present")
                #     # display(df_tmp)
                #     # print("TEMP")
                #     multiple_early_gens_present = True
                #     # break
                # __|


                # IMPORTANT <--------------------------------------------------
                #| - Select desired systems from duplicates ###################
                # Select entry from array of duplicates that will be kept
                # Can either keep the one that was already there (earliest)
                # or we can always select the most stable one

                # selected_row = \
                selected_row_early_gen = \
                    df_tmp[df_tmp["gen_acquired"] == earlist_gen].sort_values("y_real").iloc[0]

                selected_row = df_tmp.sort_values("y_real").iloc[0]
                # __|



                indices_that_are_duplicates_i = df_tmp.index.tolist()
                indices_that_are_duplicates_i.remove(selected_row.name)

                indices_that_are_duplicates.extend(indices_that_are_duplicates_i)


            # TEMP
            print("'6r716sxr9t' in indices_that_are_duplicates",
                "6r716sxr9t" in indices_that_are_duplicates)




            # TEMP
            self.TEMP_df_tmp_dict = TEMP_df_tmp_dict
            self.duplicate_swap_lists = duplicate_swap_lists

            indices_that_are_duplicates = list(set(indices_that_are_duplicates))

            self.indices_that_are_duplicates = \
                indices_that_are_duplicates + prev_duplicate_ids


            # [i for i in all_ids_from_duplicate_analysis if i not in indices_that_are_duplicates]

                # | - OLD | Trying to replace value for lowest energy duplicate
                # lowest_y_row = df_tmp.sort_values("y_real").iloc[0]
                # TEMP
                # lowest_y_row = df_tmp.sort_values("y_real").iloc[1]
                # if earliest_acq_row.name != lowest_y_row.name:
                #     print(earliest_acq_row.name, lowest_y_row.name)
                # #     model.loc[lowest_y_row.name]
                # #     model.rename(
                # #         index={
                # #             lowest_y_row.name: earliest_acq_row.name + "_TEMP",
                # #             earliest_acq_row.name: lowest_y_row.name,
                # #             }, inplace=True)
                # #     model.rename(
                # #         index={
                # #             earliest_acq_row.name + "_TEMP": earliest_acq_row.name,
                # #             }, inplace=True)
                # __|

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
