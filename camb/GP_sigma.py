from sklearn.gaussian_process import GaussianProcessRegressor
import sys
from sklearn.gaussian_process.kernels import RBF
import numpy as np
import matplotlib.pyplot as plt
import argparse



parser = argparse.ArgumentParser(description='Input parameters from CAMB')
parser.add_argument('--inia',metavar='aini',  type=float, nargs='+',
                   help='initial scale factor')
parser.add_argument('--enda',metavar='aend',  type=float, nargs='+',
                   help='end scale factor')
parser.add_argument('--ODEsteps',metavar='ODEsteps',  type=int, nargs='+',
                   help='number of steps for the ODE solver')
parser.add_argument('--scalefactors',metavar='a',  type=float, nargs='*',default=[],
                   help='values of redshifts')
parser.add_argument('--sigma',metavar='sigma',  type=float, nargs='*',default=[],
                   help='equation of state')
parser.add_argument('--sigman',metavar='sn',  type=float, nargs='*',default=[],
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
inia = args.inia[0]
enda = args.enda[0]
ODEsteps = args.ODEsteps[0]
#z_edge = np.array(args.redshifts) #NH redshift at edge of each bin
sig = np.array(args.sigma)
sn = np.array(args.sigman)
l = args.l[0]
filename = args.outfile[0]

nb=len(sig)
#print nb
a = np.array(args.scalefactors)#z_edge[:-1] + np.diff(z_edge)/2 #NH z is now the redshift in the middle of each bin
ini=[inia]
a=np.concatenate([ini,a])
#print  a
sig=np.concatenate([sn,sig])
#print  sig



#GP smoothed at last point
sma=[a[nb]-0.001]
a_gp=np.concatenate([a,sma])
sms=[sig[nb]]
sig_gp=np.concatenate([sig,sms])

#defining the baseline -1
base = lambda x: +1+x-x

# Generation of the Gaussian Process
gp = GaussianProcessRegressor(kernel=RBF(l, (l, l)))

#Fit --> Training
g = gp.fit(a_gp[:, np.newaxis], sig_gp-base(a_gp))

#Plotting points (if log use np.logspace)
a_sampling = np.linspace(inia, enda, ODEsteps)
#print a_sampling

#transforming a_sampling in z_sampling
z_sampling=np.zeros(ODEsteps)
for i in range (ODEsteps):
    z_sampling[i]= -1 + 1/a_sampling[i]
#print z_sampling
#Predict points
sig_pred, sigma = gp.predict(a_sampling[:, np.newaxis], return_std=True)
sig_pred = sig_pred + base(a_sampling)


#Plot the result: remove it from final verions
#fig= plt.figure(figsize=(14,12))
#plt.plot(a_sampling, sig_pred, label = 'l=%s'%l)
#plt.legend(fontsize=20)
#plt.scatter(a, sig)
#fig.savefig('test_figure_sigma_a.png')

# print to file
#f = open(filename,'mu')
# print len(z_sampling)
#for i in range(0, ODEsteps):
#    print >>f, z_sampling[i], mu_pred[i]
np.savetxt(filename, np.array([z_sampling, sig_pred]).T, fmt="%15.8e")

#f.close()

exit()
