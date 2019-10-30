#!/usr/bin/env python

"""Module to TEMP TEMP.

Author: Raul A. Flores
"""

#| - IMPORT MODULES
import os
# import copy
#
import pickle
#
import numpy as np
# import pandas as pd
#
# # SciKitLearn
# from sklearn.decomposition import PCA
#
# # Catlearn
# from catlearn.regression.gaussian_process import GaussianProcess
# from catlearn.preprocess.clean_data import (
#     clean_infinite,
#     clean_variance,
#     clean_skewness)
# from catlearn.preprocess.scaling import standardize
#
# from IPython.display import display
#__|

# TEMP
from al_algeneration import ALGeneration
# from active_learning import ALGeneration


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

            self.duplicate_ids = list(set(self.duplicate_ids))

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