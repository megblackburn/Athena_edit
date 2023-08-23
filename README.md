# Athena_edit
Edited Athena directory with added problem for turbulent mixing layer.

## New files:
- all files in utils
- mls2.c = current working problem, mix_layer_kh.c = more specific problem (unfinished)
- tcm3.c = current townsend cooling function
- current input = tst/megKH/athinputmeg.kh

## Work needed:
- Temperature and energy equations need double checked
- check pressure calculation
- set correct variable values

## Work done:
- error with the MPIallreduce within the frame_shift function - fixed with adding #include <mpi.h>

