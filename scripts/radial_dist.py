import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from os.path import join as pjoin
import sys
from numba import jit, prange

#import pandas as pd


def get_xyz(h5,t):
    pos = h5[str(t)+'/position/']
    return pos['x'][()],pos['y'][()],pos['z'][()]

def color_map(M):
    max_M = max(M)
    min_M = min(M)
    color_indx = (M[:] - min_M)/max_M

    return color_indx

def moving_avg(x,window):
    avg = np.zeros(len(x)-window)
    for i in range(len(x)-window):
        avg[i]=np.mean(x[i:i+window])
    return avg
def sort_xy_by_m(x,y,M):
    x1=[]
    y1=[]
    x2=[]
    y2=[]
    M1 = min(M)
    for i in range(len(x)):
        if M[i] == M1:
            x1.append(x[i])
            y1.append(y[i])
        else:
            x2.append(x[i])
            y2.append(y[i])
    return x1,y1,x2,y2

@jit(nopython=True)
def sort_r_by_m(r,M):
    M1 = float(min(M))
    r1 = []
    r2 = []
    for i in range(len(r)):
        if M[i] == M1:
            r1.append(r[i])
        else:
            r2.append(r[i])
    return r1,r2

@jit(nopython=True, parallel=True)
def avg_particle_dist(x,y,z):
    N = len(x)
    r = 0
    for i in prange(N):
        for j in range(N):
            if j!=i:
                xdiff = x[i]-x[j]
                ydiff = y[i]-y[j]
                zdiff = z[i]-z[j]
                r += np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
    return r/((N-1)*(N-1))

@jit(nopython=True, parallel=True)
def closest_neigbour(x,y,z):
    N = len(x)
    r_min = np.zeros(N)
    for i in prange(N):
        r=np.zeros(N)
        for j in range(N):
            if j!=i:
                xdiff = x[i]-x[j]
                ydiff = y[i]-y[j]
                zdiff = z[i]-z[j]
                r[j]= np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
        r = np.delete(r,[i])
        r_min[i] = r.min()
    return np.mean(r_min)


                


do_b = True

try:
    DIR =sys.argv[1]
except:
    print('ERROR: Provide the data directory path')

try:
    dim =int(sys.argv[2])
except:
    print('ERROR: Provide the data directory path')

file_name = "particle"#"saga12"#"2D1000p"#"particle"#"rhoNeutral" #"P"

h5 = h5py.File(pjoin(DIR,file_name+'.hdf5'),'r')

Nt   = h5.attrs['Nt']
dp   = h5.attrs['dp']
M = h5['particle/M'][()]

tmax = Nt*dp
x0,y0,z0 = get_xyz(h5,0)
xn,yn,zn = get_xyz(h5,tmax-5)


#M = h5['0/M/']

color_indx = color_map(M)
colors = cm.rainbow(color_indx)

print("processing positions and making plot...")
fig_pos, ax_pos = plt.subplots(1,2)
fig_pos.suptitle('Position of species 1 and species 2')
x01,y01,x02,y02 = sort_xy_by_m(x0,y0,M)
ax_pos[0].scatter(x01,y01,color='black',label='species 1',marker='.')
ax_pos[0].scatter(x02,y02,color='red',label='species 2',marker = '.')
ax_pos[0].set_xlabel('x')
ax_pos[0].set_ylabel('y')
ax_pos[0].set_title('initial')
xn1,yn1,xn2,yn2 = sort_xy_by_m(xn,yn,M)
ax_pos[1].scatter(xn1,yn1,color='black',marker='.')
ax_pos[1].scatter(xn2,yn2,color='red',marker = '.')
ax_pos[1].set_xlabel('x')
ax_pos[1].set_ylabel('y')
ax_pos[1].set_title('final')
fig_pos.legend()


#ax_pos[1].scatter(xn,yn,color=colors[color_indx])

if dim == 3:
    r0 = np.sqrt(x0[:]*x0[:] + y0[:]*y0[:] + z0[:]*z0[:])
    rn = np.sqrt(xn[:]*xn[:] + yn[:]*yn[:] + zn[:]*zn[:])
elif dim == 2:
    r0 = np.sqrt(x0[:]*x0[:] + y0[:]*y0[:])
    rn = np.sqrt(xn[:]*xn[:] + yn[:]*yn[:])

fig_r, ax_r = plt.subplots(1,2)
r01,r02 = sort_r_by_m(r0,M)

print("processing distributions and making plot...")

#Maybe scale according to r so the graph is "flat"

fig_r.suptitle('Distribution along r for species 1 and 2')
ax_r[0].hist(r01,histtype = 'step', color='black',label='species 1',bins=12)
ax_r[0].hist(r02,histtype = 'step', color='red', label='species 2',bins=12)
ax_r[0].axvline(np.mean(r01),color='black')
ax_r[0].axvline(np.mean(r02),color='red')
ax_r[0].set_xlabel('r')
ax_r[0].set_ylabel('freq')
rn1,rn2 = sort_r_by_m(rn,M)
ax_r[1].hist(rn1,histtype = 'step', color='black',bins=12)
ax_r[1].hist(rn2,histtype = 'step', color='red',bins=12)
ax_r[1].axvline(np.mean(rn1),color='black')
ax_r[1].axvline(np.mean(rn2),color='red')
ax_r[1].set_xlabel('r')
ax_r[1].set_ylabel('freq')
fig_r.legend()

if do_b:
    print("processing b and making plot...")
    b=[]

    b_step = 50
    for t in range(0,tmax-b_step,b_step):
        x,y,z = get_xyz(h5,t)
        if dim == 2:
            r = np.sqrt(x[:]*x[:] + y[:]*y[:])
        elif dim == 3:
            r = np.sqrt(x[:]*x[:] + y[:]*y[:] + z[:]*z[:])
        b1,b2 = sort_r_by_m(r,M)
        b.append(np.mean(b2)/np.mean(b1))

    fig_b,ax_b = plt.subplots()
    fig_b.suptitle('Ratio of average r for species 2 over species 1')
    ax_b.plot(b, color='black', label='raw')
    #ax_b.plot(moving_avg(b,len(b)//40), color='red', label='mean')
    ax_b.axhline(1,color='blue')
    ax_b.set_xlabel('time')
    ax_b.set_ylabel('b2/b1')
    fig_b.legend()

print("processing avg particle distance and making plot...")
r_avg=[]
r_min=[]

b_step = 50
for t in range(0,tmax-b_step,b_step):
    x,y,z = get_xyz(h5,t)
    r_avg.append(avg_particle_dist(x,y,z))
    r_min.append(closest_neigbour(x,y,z))


fig_ravg,ax_ravg = plt.subplots(1,2)
fig_ravg.suptitle('avg interparticle distance')
ax_ravg[0].plot(r_avg, color='black', label='raw')
ax_ravg[1].plot(r_min, color='black', label='raw')
#ax_b.plot(moving_avg(b,len(b)//40), color='red', label='mean')
ax_ravg[0].set_xlabel('time')
ax_ravg[0].set_ylabel('$\hat{r}$')
ax_ravg[1].set_xlabel('time')
ax_ravg[1].set_ylabel('$\hat{r_{min}}$')
plt.tight_layout()
fig_ravg.legend()

plt.show()