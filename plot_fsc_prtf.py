import matplotlib.pyplot as mlp
import numpy as np

fsc  = np.loadtxt("fsc-0-0.dat",usecols = (0,2))
prtf = np.loadtxt("data/recon_avg/4et8_full-prtf.dat")

mlp.figure()
mlp.scatter(fsc[:,0],fsc[:,1])
mlp.title("FSC")
mlp.savefig("fsc_plot.pdf",format="pdf")
mlp.close()

mlp.figure()
mlp.scatter(prtf[:,0],prtf[:,1])
mlp.title("PRTF")
mlp.savefig("prtf_plot.pdf",format="pdf")
mlp.close()

del fsc
del prtf
