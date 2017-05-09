INFORMATION ON TEST

GOAL: Run a simple m file on cluster. This mfile displays 'hello_world' and writes in in OUTPUT file.

SUCCESS: Yes

RUN TEST:
    > make_oar_octave easy_test 0 1
    > oarsub -S ./easy_test.oar
