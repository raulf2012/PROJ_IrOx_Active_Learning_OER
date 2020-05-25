# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# #############################################################################
import pickle; import os
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs",
    "out_data/AB3/gp_ucb_True",

    # "TEST_AL_2_devehowo.pickle",
    "TEST_AL_2_fugunefo.pickle",
    )
with open(path_i, "rb") as fle:
    AL = pickle.load(fle)
# #############################################################################

CS = AL.CandidateSpace
FP = CS.FingerPrints
PCA = FP.PCA

# +
print("PCA variance explained:", sum(PCA.explained_variance_))

print("num components:", PCA.n_components)
