import os
import numpy as np

path='./'
with open(path+file, "r") as fid_read:
  with open(path+file, "w") as fid_write:
    for line in fid_read:
      if "AthenaArray" in line:
        line = "calloc_1d_array()
