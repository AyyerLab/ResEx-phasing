#!/usr/bin/python

# This function converts the intensity file to a symmetric form by intoducing three symmetry axis: X,Y and Z in the middle of the intensity 3D volume. To symmetrize we replace the intensity value by the average of its symemtric counterparts.

# Pass the filename with the path as an argument to the command line.

#Author: Moshir Harsh, FSC-CFEL1, DESY. btemoshir@gmail.com

import numpy as np
import sys
import os

file_name= sys.argv[1]
mat = np.fromfile(file_name,dtype='float32')
mat_reshape = np.reshape(mat,(301,301,301),order='C')

#Replace each term with its and its reflection term's average
mat_reshape = (mat_reshape + mat_reshape[::-1])/2
mat_reshape = (mat_reshape + mat_reshape[:,::-1])/2
mat_reshape = (mat_reshape + mat_reshape[:,:,::-1])/2

#Export the file
mat_reshape.tofile(os.path.splitext(file_name)[0]+ "-sym" + os.path.splitext(file_name)[1])

