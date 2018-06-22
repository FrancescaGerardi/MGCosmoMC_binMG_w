from sklearn.gaussian_process import GaussianProcessRegressor
import sys
from sklearn.gaussian_process.kernels import RBF
import numpy as np
import matplotlib.pyplot as plt
import argparse



parser = argparse.ArgumentParser(description='Input parameters from CAMB')
parser.add_argument('--inired',metavar='zini',  type=float, nargs='+',
                   help='initial redshift')
parser.add_argument('--endred',metavar='zend',  type=float, nargs='+',
                   help='end redshift')
parser.add_argument('--ODEsteps',metavar='ODEsteps',  type=int, nargs='+',
                   help='number of steps for the ODE solver')
parser.add_argument('--redshifts',metavar='z',  type=float, nargs='*',default=[],
                   help='values of redshifts')
parser.add_argument('--eos',metavar='wde',  type=float, nargs='*',default=[],
                   help='equation of state')
parser.add_argument('--eos0',metavar='w0',  type=float, nargs='*',default=[],
                   help='equation of state at z=0')
parser.add_argument('--l',metavar='l',  type=float, nargs='+',
                   help='correlation length')
# parser.add_argument('--outfile', nargs='+', type=argparse.FileType('w'),default=sys.stdout)
parser.add_argument('--outfile', nargs='+', type=str ,default=sys.stdout)
args = parser.parse_args()

#print args.outfile
#print args.inired[0], args.endred[0], args.ODEsteps[0]
#print args.redshifts                                   valori in mezzo ai bin
#print args.eos
#print args.l[0]
#print args.outfile[0]

#Training points
inired = args.inired[0]
endred = args.endred[0]
ODEsteps = args.ODEsteps[0]
#z_edge = np.array(args.redshifts) #NH redshift at edge of each bin
wde = np.array(args.eos)
w0 = np.array(args.eos0)
l = args.l[0]
filename = args.outfile[0]

z = np.array(args.redshifts)#z_edge[:-1] + np.diff(z_edge)/2 #NH z is now the redshift in the middle of each bin
ini=[inired]
z=np.concatenate([ini,z])
#print  z
wde=np.concatenate([w0,wde])
#print  wde

nb=len(wde)
#defining the baseline -1
base = lambda x: -1+x-x

# Generation of the Gaussian Process
gp = GaussianProcessRegressor(kernel=RBF(l, (l, l)))

#Fit --> Training
g = gp.fit(z[:, np.newaxis], wde-base(z))

#Plotting points (if log use np.logspace)
z_sampling = np.linspace(inired, endred, ODEsteps)
#q_sampling = np.linspace(-1., -0.1, 1000)

#Predict points
w_pred, sigma = gp.predict(z_sampling[:, np.newaxis], return_std=True)
w_pred = w_pred + base(z_sampling)

#Plot the result: remove it from final verions
fig= plt.figure(figsize=(14,12))
plt.plot(z_sampling, w_pred, label = 'l=%s'%l)
plt.legend(fontsize=20)
plt.scatter(z, wde)
fig.savefig('test_figure.png')

# print to file
#f = open(filename,'w')
# print len(z_sampling)
#for i in range(0, ODEsteps):
#    print >>f, z_sampling[i], w_pred[i]
np.savetxt(filename, np.array([z_sampling, w_pred]).T, fmt="%15.8e")

#f.close()

exit()
