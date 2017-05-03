INFORMATION ON TEST

GOAL: reconstruct one silce (FBP) of a volume using PySHT2 in cluster. Output volume in OUTPUT file.

SUCCESS: Yes

RUN TEST:
    > make_oar_pyhst2 phantom_al_mg_PET_thibaut_1dist_1_slice.par

    > oarsub -S ./phantom_al_mg_PET_thibaut_1dist_1_slice.oar
