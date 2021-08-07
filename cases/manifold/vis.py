import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from backwardc import floodsz_backwardc2
from backwardc import norm2angle
from backwardc import angle2norm
from backwardc import normalization1

p = np.arange(-np.pi,np.pi,np.pi/90)
t = np.arange(0,np.pi,np.pi/90)
P,T = np.meshgrid(p, t)
mesh = np.array([P,T])
pt = mesh.T.reshape(-1, 2)

n = len(pt)
norms = np.empty([n,3])
cs = np.empty([n,3])
pars = np.empty([n])
for i in range(len(pt)):
    norms[i,:] = angle2norm(pt[i,:])
    norms[i,:] = norms[i,:] / np.sum(np.abs(norms[i,:]))
    cs[i,:],pars[i]  = floodsz_backwardc2(norms[i,:],0.05)

xs = cs[:,0]
ys = cs[:,1]
zs = cs[:,2]
xs1 = xs[pars==1]; ys1=ys[pars==1]; zs1=zs[pars==1]; par1=pars[pars==1]
xs2 = xs[pars==2]; ys2=ys[pars==2]; zs2=zs[pars==2]; par2=pars[pars==2]
xs3 = xs[pars==3]; ys3=ys[pars==3]; zs3=zs[pars==3]; par3=pars[pars==3]
xs4 = xs[pars==4]; ys4=ys[pars==4]; zs4=zs[pars==4]; par4=pars[pars==4]
xs5 = xs[pars==5]; ys5=ys[pars==5]; zs5=zs[pars==5]; par5=pars[pars==5]

phs = pt[:,0]
ths = pt[:,1]
ph1 = phs[pars==1]; th1=ths[pars==1]
ph2 = phs[pars==2]; th2=ths[pars==2]
ph3 = phs[pars==3]; th3=ths[pars==3]
ph4 = phs[pars==4]; th4=ths[pars==4]
ph5 = phs[pars==5]; th5=ths[pars==5]

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter(xs,ys,zs,c=pars,cmap='viridis')

plt.show()
