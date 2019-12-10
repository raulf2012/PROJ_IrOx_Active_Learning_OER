#!/usr/bin/env python

"""Module to TEMP TEMP.

Author: Raul A. Flores
"""

#| - IMPORT MODULES
import os
import sys

#__|


class CCF:
    """
    """

    #| - CCF ******************************************************************
    _TEMP = "TEMP"


    def __init__(self,
        df_dij=None,
        d_thresh=None,
        ):
        """
        """
        #| - __init__

        #| - Setting Argument Instance Attributes
        self.df_dij = df_dij
        self.d_thresh = d_thresh
        #__|

        #| - Initializing Internal Instance Attributes

        #__|

        #__|

    def i_all_similar(self, index_i, filter_ids=None):
        """Check whether there are any systems similar to i.
        """
        #| - i_all_similar
        # self = CCF
        d_thresh = self.d_thresh
        df_dij = self.df_dij

        # index_i = "mkbrzh8kv5"
        # index_j = "folatese_05"
        # filter_ids = all_indices
        # filter_ids = None
        # #####################################################################

        index_in_df = index_i in df_dij.index
        if not index_in_df:
            return(None)


        if filter_ids is not None:
            filter_ids_inter = df_dij.index.intersection(filter_ids)
            df_dij = df_dij.loc[filter_ids_inter, filter_ids_inter]

        row_i = df_dij.loc[index_i]
        row_i = row_i.drop(labels=index_i)

        out_dict = row_i[row_i < d_thresh].to_dict()

        return(out_dict)
        #__|

    def i_j_similar(self, index_i, index_j):
        """Test whether systems i and j are similar to one another.
        """
        #| - i_j_similar
        # self = CCF
        d_thresh = self.d_thresh
        df_dij = self.df_dij

        # index_i = "budabebu_36"
        # index_j = "folatese_05"
        # #####################################################################

        dij = df_dij.loc[index_i, index_j]

        similar = False
        if dij < d_thresh:
            similar = True
        elif dij > d_thresh:
            similar = False
        else:
            assert False, "AHHHHHHHH!!!"

        return(similar)
        #__|

    #__| **********************************************************************
