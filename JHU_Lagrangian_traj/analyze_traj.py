#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 13:13:50 2018

@author: ron
"""

import numpy as np
import matplotlib.pyplot as plt


import matplotlib
cmap = matplotlib.cm.get_cmap('viridis')

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

tr = np.load('trajectories_t=4_steps=150.npz')

pos, time, vel = tr['x'], tr['t'], tr['u']
dt = time[1]-time[0]
N = pos.shape[1]
t_eta = 0.045
eta = 0.00287
puffs = np.array(zip(range(0,8000,400),range(400,8400,400)))


def get_r_lists(pos, vel, N_start, N_end, gamma_th):
    eps = 0.092 # <- dissipation
    eta, L = 0.00287, 1.376
    r_lst = []
    for i in range(N_start,N_end):
        for j in range(i+1,N_end):
            r = np.linalg.norm(pos[:,i,:] - pos[:,j,:], axis=1)
            
            #if r[0] < 50*eta: continue
            #r0 = (pos[0,i,:] - pos[0,j,:])/r[0]
            #v0 = np.dot( vel[i,:] - vel[j,:] , r0 )
            
            v0 = (r[1] - r[0])/dt
            
            gamma = v0**2 / (r[0] * eps)**(0.666)
            if gamma < gamma_th:
                r_lst.append(r)

    print len(r_lst)
    r_lst = np.array(r_lst)
    dr_lst = r_lst[:,:] - r_lst[:,0:1]
    dr2_lst = (r_lst[:,:]/r_lst[:,0:1] - 1)**2 
    
    return r_lst, dr_lst, dr2_lst




# averaging over puffs:
    
g = [0.01,0.1,1.0,10]
curves = []
p = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19] # < -- 15 is a BAD PUFF :/
for e,gamma in enumerate(g):
    dr2_lst_lst = []
    for seq in puffs[p,:]:
        r_lst, dr_lst, dr2_lst = get_r_lists(pos, vel, seq[0], seq[-1], gamma)
        dr2_lst_lst.append(np.mean(dr2_lst,axis=0))
    av_r2 = np.mean(dr2_lst_lst,axis=0)
    #ax.loglog(time, av_r2,lw=1,c=cmap(1.0*e/len(g)), label = r'$\gamma=%.2f$'%gamma)
    curves.append(av_r2)


curves = np.array(curves) / (time/t_eta)**3 
np.savez('r2_DNS_curves.npz', data = curves, time_ax = time/t_eta)


fig, ax = plt.subplots()
TL = 1.99

ax.loglog([TL/t_eta,TL/t_eta], [0.3e-4,0.3e-3],'k--',lw=1)
ax.text(TL/t_eta*0.90, 0.35e-3,r'$T_L$',size=18)
#ax.annotate(r'$T_L$', xy=(TL/t_eta, 0.5e-4), xytext=(TL/t_eta, 0.5e-3),
#            arrowprops=dict(facecolor='black',width=2))

for e,c in enumerate(curves):
    ax.loglog(time/t_eta, c / (time/t_eta)**3, lw=2, c=cmap(1.0*e/len(g)),
              label = r'$\gamma=%.2f$'%(g[e]))
ax.legend()
ax.set_xlabel(r'$\tau / \tau_\eta$',size=18)
ax.set_ylabel(r'$\langle  \Delta r^2 / r_0^2 \rangle \cdot (\tau / \tau_\eta)^{-3}$',size=18)
ax.grid()
plt.tight_layout()
fig.savefig('normed_r2_dissipation.pdf')






# ======================================
# ================ mixing between puffs:
# ======================================
N = 400
indexes = []
for i in range(N):
    test = 1
    while test:
        rnd = np.random.randint(0,8000)
        test = rnd in indexes
    indexes.append(rnd)

rand_pos = pos[:,indexes,:]



g = [0.01,0.1,1.0,10]
curves = []
for e,gamma in enumerate(g):
    r_lst, dr_lst, dr2_lst = get_r_lists(rand_pos, vel, 0, N, gamma)
    av_r2 = np.mean(dr2_lst_lst,axis=0)
    curves.append(av_r2)


fig, ax = plt.subplots()
TL = 1.99

ax.loglog([TL/t_eta,TL/t_eta], [0.5e-4,1e-3],'k--')

for e,c in enumerate(curves):
    ax.loglog(time/t_eta, c / (time/t_eta)**3, lw=2, c=cmap(1.0*e/len(g)),
              label = r'$\gamma=%.2f$'%(g[e]))
ax.legend()
ax.set_xlabel(r'$\tau / \tau_\eta$',size=18)
ax.set_ylabel(r'$\langle  \Delta r^2 / r_0^2 \rangle \cdot (\tau / \tau_\eta)^{-3}$',size=18)
ax.grid()







# ==============================
# get gamma distributions:
#==============================



def get_gamma_list(pos,N_start, N_end):
    eps = 0.092 # <- dissipation
    gamma_lst = []
    r0_lst = []
    for i in range(N_start,N_end):
        for j in range(i+1,N_end):
            r = np.linalg.norm(pos[:,i,:] - pos[:,j,:], axis=1)
            v0 = (r[1]-r[0])/dt
            gamma = v0**2 / (r[0] * eps)**(0.666)
            gamma_lst.append(gamma)
            r0_lst.append(r[0])
    return gamma_lst, r0_lst


g_lst, r0_lst = [] , []
for seq in range(0,8000,400):
    g_, r0_ = get_gamma_list(pos, seq, seq + 400)
    for g in g_: g_lst.append(g)
    for r0 in r0_: r0_lst.append(r0)

res = zip(g_lst, r0_lst) 
    
eta = 0.00287
fig, ax = plt.subplots()
shp = ['o','d','^','v']

for e,r0 in enumerate([100.0,35.0,20.0,10.0]):
    temp = []
    for i in res:
        if i[1]<r0*eta:
            temp.append(i[0])
    print len(temp)
    h = np.histogram(temp,bins=np.logspace(-2,1.4772,num=25), normed=True)
    x,y = 0.5*(h[1][:-1] + h[1][1:]) , h[0]
    ax.plot(x,y, shp[e]+'-',color=colors[e], label = r'$r_0 < %d \eta $'%(r0))
    #ax.semilogx(x[::5],y[::5], label = r'$r_0 < %d \eta $'%(r0))
    
ax.legend()
ax.grid()
ax.set_xlabel(r'$\gamma$')
ax.set_ylabel(r'$PDF(\gamma | r_0)$')
ax.set_xlim([-0.5,15.5])
ax.set_ylim([1e-2,4e-1])
fig.savefig('gamma_distribution.pdf')





















