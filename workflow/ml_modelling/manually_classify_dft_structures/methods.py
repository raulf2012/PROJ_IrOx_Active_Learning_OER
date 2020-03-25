# | - Import Modules
from pymatgen.analysis.local_env import (
    LocalStructOrderParams, get_neighbors_of_site_with_index)
#__|

def site_is_of_motif_type(
    struct,
    n,
    neighbors_list=None,

    approach="min_dist",
    delta=0.1,
    cutoff=10.0,
    thresh=None):
    """
    Returns the motif type of the site with index n in structure struct;
    currently featuring "tetrahedral", "octahedral", "bcc", and "cp"
    (close-packed: fcc and hcp) as well as "square pyramidal" and
    "trigonal bipyramidal".  If the site is not recognized,
    "unrecognized" is returned.  If a site should be assigned to two
    different motifs, "multiple assignments" is returned.

    Args:
        struct (Structure): input structure.
        n (int): index of site in Structure object for which motif type
                is to be determined.
        approach (str): type of neighbor-finding approach, where
              "min_dist" will use the MinimumDistanceNN class,
              "voronoi" the VoronoiNN class, "min_OKeeffe" the
              MinimumOKeeffe class, and "min_VIRE" the MinimumVIRENN class.
        delta (float): tolerance involved in neighbor finding.
        cutoff (float): (large) radius to find tentative neighbors.
        thresh (dict): thresholds for motif criteria (currently, required
                keys and their default values are "qtet": 0.5,
                "qoct": 0.5, "qbcc": 0.5, "q6": 0.4).

    Returns: motif type (str).
    """
    # | - site_is_of_motif_type
    # print("THIS IS MY CUSTOM METHOD, NOT NATIVE PYMATGEN")

    if thresh is None:
        thresh = {
            "qtet": 0.5, "qoct": 0.5, "qbcc": 0.5, "q6": 0.4,
            "qtribipyr": 0.8, "qsqpyr": 0.8}

    ops = LocalStructOrderParams([
        "cn", "tet", "oct", "bcc", "q6", "sq_pyr", "tri_bipyr"])

    # #########################################################################
    if neighbors_list is None:
        neighs_cent = get_neighbors_of_site_with_index(
            struct, n, approach=approach, delta=delta, cutoff=cutoff)
    else:
        neighs_cent = neighbors_list

    neighs_cent.append(struct.sites[n])

    opvals = ops.get_order_parameters(
        neighs_cent, len(neighs_cent) - 1, indices_neighs=[
            i for i in range(len(neighs_cent) - 1)])
    cn = int(opvals[0] + 0.5)
    motif_type = "unrecognized"
    nmotif = 0

    if cn == 4 and opvals[1] > thresh["qtet"]:
        motif_type = "tetrahedral"
        nmotif += 1
    if cn == 5 and opvals[5] > thresh["qsqpyr"]:
        motif_type = "square pyramidal"
        nmotif += 1
    if cn == 5 and opvals[6] > thresh["qtribipyr"]:
        motif_type = "trigonal bipyramidal"
        nmotif += 1
    if cn == 6 and opvals[2] > thresh["qoct"]:
        motif_type = "octahedral"
        nmotif += 1
    if cn == 8 and (opvals[3] > thresh["qbcc"] and opvals[1] < thresh["qtet"]):
        motif_type = "bcc"
        nmotif += 1
    if cn == 12 and (opvals[4] > thresh["q6"] and opvals[1] < thresh["q6"] and
                     opvals[2] < thresh["q6"] and opvals[3] < thresh["q6"]):
        motif_type = "cp"
        nmotif += 1

    if nmotif > 1:
        motif_type = "multiple assignments"

    return motif_type
    #__|
