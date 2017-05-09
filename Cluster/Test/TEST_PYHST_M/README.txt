INFORMATION ON TEST

GOAL: Use an octave script to reconstruct one silce (FBP) of a volume using PySHT2 in cluster. Output volume in OUTPUT file.

SUCCESS: Yes

RUN TEST:
    > make_oar_octave test 1 1

    > oarsub -S ./test.oar
