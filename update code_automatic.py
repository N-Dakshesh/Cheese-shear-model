'''
Adding equilibration steps
'''
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import disk as disk

F_tensor = 0


#Only depends on d, distance between particles.
def Mie(d, e = 1, sigma = 1, n=6, m=3):
    factor = (n/(n-m))*(n/m)**(m/(n-m))
    main = (sigma/d)**n - (sigma/d)**m
    return factor*e*main

#Derivative
def Mie_diff(d, e = 1, sigma = 1, n = 6, m = 3):
    factor = (n/(n-m))*(n/m)**(m/(n-m))
    main = -n*(sigma/d)**(n+1) + m*(sigma/d)**(m+1)
    return factor*e*main

#Input x, positions of points
def distances(x):
    dx = x[:, 0] - x[:, 0][:, None]
    dy = x[:, 1] - x[:, 1][:, None]
    dz = x[:, 2] - x[:, 2][:, None]
    #r = np.sqrt(dx*dx + dy*dy + dz*dz)
    #return r, dx, dy, dz

    r = pdist(x)
    r = np.nan_to_num(r)

    U = np.sum(Mie(x))
    #r[r>length/2] = 0
    r1 = squareform(r)
    return r1, dx, dy, dz, U #Also return U, potential


'''
def __distances(X):
    # --> Number of particles.
    n = len(X)
    
    # --> Gram matrix.
    d2 = -2 * X @ X.T
    
    # --> Squared pairwise distances.
    diag = -0.5 * np.einsum('ii->i', d2)
    d2 += diag + diag[:, None]
    
    # --> Prevent division by 0.
    np.einsum('ii->i', d2)[...] = 1
    
    # --> Net forces.
    F = np.nansum( (X[:, None, :] - X) * d2[..., None]**-1.5, axis=0)
    
    return d2
'''


def norm(v):
    #input vector
    return v/np.linalg.norm(v)

def unit_dir(dx, dy, dz):
    #directions = np.dstack((dx, dy))
    #directions = np.dstack((directions, dz))
    al1 = np.array([dx, dy, dz])
    al2 = al1.transpose(1, 2, 0)
    directions = al2

    #To normalise
    #Replace nans
    unit_directions = directions/np.linalg.norm(directions, axis=2)[:,:,None]
    unit_directions = np.nan_to_num(unit_directions)
    return unit_directions



#slow don't use
def _stress(directions, F_tensor_per_particle):
    #Need to sum F per particle*distance to that particle.
    #3x3 array for 3x3 dimensions?
    #dirs = unit_dir*r
    sigma = np.zeros((3, 3))
    for i in range(0, len(directions)):
        for j in range(0, len(directions)):
            for k in range(0, 3):
                for l in range(0, 3):
                    #Need to get correct line from direction. k and l refer to dimension, i and j refer to particles.
                    sigma[k, l] = directions[i, j, k]*F_tensor_per_particle[i, j, l]
    return sigma
'''
def stress(r, unit_dir, F_tensor_per_particle):
    dirs = r*unit_dir
    sigma = np.einsum('ijk,ijl->kl', dirs, F_tensor_per_particle)
    return sigma
'''






def image_copy(x):
    '''First N rows (where N is number of particles) are original positions! Save these.'''
    all_x = x


    shifts = np.stack(list(product([-1, 0, 1], repeat = 3)))
    #print(shifts)

    for i in range(0, len(shifts)):
        x1 = np.copy(x)
        for j in range(0, 3):

            check = shifts[i] == np.array([0,0,0])
            if check.all()==True:
                int(1)
            else:
                x1[:, j] = x1[:, j] + 2*shifts[i,j]

        check = shifts[i] == np.array([0,0,0])
        if check.all()==True:
            int(1)
        else:
            all_x = np.vstack((all_x, x1))


    return all_x

#With this

def F_per_particle(all_x):
    #Get 'image' positions for forces
    #x_map = np.memmap('x.mymemmap', dtype = 'float64', shape = np.shape(x))
    r = distances(all_x)
    dis = r[0]
    dis[ dis < 0.7] = 0
    unit_directions = unit_dir(r[1], r[2], r[3]) #unit vectors of forces
    directions = np.dstack((r[1], r[2]))
    directions = np.dstack((directions, r[3])) #Need non unit directions
    #unit_directions[abs(unit_directions) < 1] = 0

    F_mag = Mie_diff(dis)
    F_tensor_per_particle = F_mag[:,:,None] * unit_directions
    F_tensor_per_particle = np.nan_to_num(F_tensor_per_particle)

    #sigma = np.einsum('ijk,ijl->kl', directions, F_tensor_per_particle) #This gives stress matrix. Need distances*forces
    return F_tensor_per_particle, directions, r[4]

def stress(F_tensor_per_particle, all_v, directions):
    term2 = 0.5*np.einsum('ijk,ijl->kl', directions, F_tensor_per_particle)
    term1 = np.einsum('ia,ib->ab',all_v, all_v)
    term1 = 0
    sigma = term1 + term2
    sigma = sigma/(dims*2*width*height)
    return sigma


def F_update(F_tensor_per_particle):
    F_tensor = F_tensor_per_particle.sum(axis = 1)
    return F_tensor



def shift(x):
    '''input stepped array. Shift across'''
    #Need to add the initial shift in
    x[x > dims + 3*dims] -= 2*dims
    x[x < -dims + 3*dims] += 2*dims
    if (x < -dims + 3*dims).any():
        print('shift')
    if (x > dims + 3*dims).any():
        print('shift')
    return x

def step_equil(all_x, all_v, gamma, dt = 0.1, F_tensor = F_tensor, i = 0): #This doesn't strain
    '''all_x is x with copied unit cells'''
    #Pass through step function with argument
    dv = F_tensor*dt*dt

    all_v = all_v + dv
    dr = dv*dt
    dr[:,0] = dr[:,0]
    all_x = all_x + dr

    return all_x, all_v

def step_strain(all_x, all_v, gamma, dt = 0.1, F_tensor = F_tensor, i = 0): #This does strain
    '''Need to step again after equilibration'''
    #Pass through step function with argument
    g = gamma
    dv = F_tensor*dt*dt
    dv[:, 0] = dv[:,0] + g*dt*all_x[:,0]
    all_v = all_v + dv
    dr = dv*dt
    dr[:,0] = dr[:,0]
    all_x = all_x + dr

    return all_x, all_v

def many_step_equil(steps, x, step_func = step_equil, gamma = 0, dt = 0.1):
    #time x is x points evolving with time
    time_x = x
    #x positions at the current point in time
    new_x = x
    new_v = np.zeros((27*len(x), 3))
    #new_v = np.vstack((v, new_v)) #Need for every particle due to ghost PBC. Combine
    all_v = np.zeros((len(x), 3)) #Final v's only matter for inside central box
    sigma = np.zeros((3, 3))
    all_sigma = np.zeros((3, 3))
    a = -1 #Makes sure you get number of steps under strain you want
    i = 0 #Iteration number
    all_U = np.array([])

    #Need to update F each time!
    while a < steps:
        try:
            copy_x = image_copy(new_x) #create periodic cell arrangement
            new_F_per, directions, new_U = F_per_particle(copy_x)#Find forces for each particle, and get direction vectors out
            b=0
            all_U= np.append(all_U, new_U)
            sigma_before = np.diag(sigma).sum() #Get stress before step

            sigma = stress(new_F_per, new_v, directions)
            #print(sigma)
            new_F = F_update(new_F_per) #Sum across particles
            #print('hi')
            copy_x_update, new_v = step_func(copy_x, new_v, gamma, dt, new_F, i) #Step forwards
            new_x = copy_x_update[0:np.shape(x)[0]] # delete particles that have moved outside bound
            #if i ==84:
            #    int(1)
            new_x = shift(new_x)
            time_x = np.dstack((time_x, new_x))
    
            sigma_now = np.diag(sigma).sum() #get stress after step
    
    
    
            all_sigma = np.dstack((all_sigma, sigma))
    
            all_v = np.dstack((all_v, new_v[0:np.shape(x)[0]]))
            print(i, sigma_now)
            i = i+1
    
            if step_func==step_strain:
                a = a+1
                #print('hi')
    
            if abs(sigma_now - sigma_before) < 1e-5: #Equilibrium if sigma is constant
                print('equilibrium')
                b = i #When you reached equilibrium
                step_func = step_strain
        except KeyboardInterrupt:
            return time_x, all_v, all_sigma, copy_x_update, new_v, all_U



    all_sigma = all_sigma[:, :, 1:] #Remove initial 0's.
    return time_x, all_v, all_sigma, copy_x_update, new_v, all_U #Return positions of initial particles at every point every t
                                                          #Return velocities of initial particles at every point every t
                                                          #Return stress tensor at every point every t
                                                          #Return full ghost positions
                                                          #Return full ghost velocities. These are needed if you want to start the system from a specific point.



#%%
import random
r=0.1
length = 3
width = 3
height = 3
grid = disk.Grid(r, length, width, height)

rand = (random.uniform(0, length), random.uniform(0, width), random.uniform(0, height))
#Seed creates repeatable points
data = grid.poisson(rand, 30)
x = np.array(data)

x = x + length
dims = length/2

#%%
gamma = 0.01
dt = 0.1
results, v, sigma, ghost_final_x, ghost_final_v, all_U = many_step_equil(1000, x, step_equil, gamma, dt)
#To get diagonals
sigma1 = np.einsum('ijk->k', sigma)
#%%
kinetic = v*v
kinetic = kinetic.sum(axis = 1)
kinetic = kinetic.sum(axis = 0)
kinetic = kinetic*0.5
energy = kinetic[1:] + all_U


#%%
t = np.arange(0, dt*len(sigma1), 0.1)
plt.plot(t, sigma1)
plt.xlabel('time, reduced units')
plt.ylabel('Stress, reduced units')
plt.title('Euler, stress - strain')
#plt.xlim(49.5,50)


plt.show()
#%%
plt.figure(figsize = [12, 8])
plt.plot(t[1:], energy)
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('Energy against time for Euler')

#%%
plt.plot(0.1*t, sigma1)
plt.xlabel('Strain, reduced units')
plt.ylabel('Stress, reduced units')
#plt.xlim(0.3, 0.5)
plt.show()

#%%
np.savetxt('Final_results4.csv', sigma1)
#np.savetxt('after_equilibration6', results[:,:,-1])
#np.savetxt('after_equilibrationv6', v[:,:,-1])

