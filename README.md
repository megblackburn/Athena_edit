# Athena_edit
Athena fork with new problem for turbulent mixing layer.

New files:
- all files in utils
- mls.c is current working problem
- tcm3.c is the current townsend cooling function
- current input is tst/megKH/athinputmeg.kh


Work needed:
- error with the MPIallreduce within the frame_shift function - init_mesh.c file needs changed
- Temperature and energy equations need double checked
- definitions of U[k][j][i]. need checked
- upper and lower densities need added and assigned to U[k][j][i].d for different requirements - along with temperatures 
- check pressure calculation
- problem needs refined to specific TML requirements
- set correct variables values

