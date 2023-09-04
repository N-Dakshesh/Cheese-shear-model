'''
Adding equilibration steps
Fixing 2
Different velocity verlet
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
        x = shift(x)
    if (x > dims + 3*dims).any():
        print('shift')
        x = shift(x)
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
    dv = F_tensor*dt
    dv[:, 0] = dv[:,0] + g*all_x[:,0]
    all_v = all_v + dv
    dr = dv*dt
    dr[:,0] = dr[:,0]
    all_x = all_x + dr

    return all_x, all_v


def position_update(x,v,F,dt, gamma):
    x_new = x + dt*v
    return x_new

def position_update_stress(x,v,F,dt, gamma):
    x_new = x + dt*v
    x_new[:,0] = x_new[:,0] + gamma*x_new[:,0]*dt
    return x_new

def velocity_update(v,F_new,dt):
    v_new = v + 0.5*dt*(F_new)
    return v_new

def force_n_stress(copy_x, all_v):
    new_F_per, directions, new_U = F_per_particle(copy_x)
    sigma = stress(new_F_per, all_v, directions)
    new_F = F_update(new_F_per)
    return new_F, new_U, sigma

def many_step_equil(steps, x, step_func = position_update, gamma = 0, dt = 0.1):
    equil = False #Equilibrium condition
    #time x is x points evolving with time
    time_x = x
    #x positions at the current point in time
    new_x = x
    new_v = np.zeros((27*len(x), 3))
    #new_v = np.vstack((v, new_v)) #Need for every particle due to ghost PBC. Combine
    all_v = np.zeros((len(x), 3)) #Final v's only matter for inside central box
    sigma = np.zeros((3, 3))
    all_sigma = np.zeros((3, 3))
    all_U = np.array([])
    a = -1 #Makes sure you get number of steps under strain you want
    i = 0 #Iteration number
    new_F = 0


    #Need to update F each time!
    while a < steps:
        try:
            new_v = velocity_update(new_v, new_F, dt)
            copy_x = image_copy(new_x) #create periodic cell arrangement
            copy_x = step_func(copy_x, new_v, new_F, dt, gamma) #Updating positions using most recent force
            old_F = new_F #Shift forces. new force is now old F
            sigma_before = np.diag(sigma).sum() #Get stress before step

            new_F, new_U, sigma = force_n_stress(copy_x, new_v) #Updating forces. We now have new and old force
            new_v = velocity_update(new_v, new_F, dt) #Updating velocities
            all_U = np.append(all_U, new_U)


            new_x = copy_x[0:np.shape(x)[0]] # delete particles that have moved outside bound
            new_x = shift(new_x) #Apply boundary conditions
            time_x = np.dstack((time_x, new_x))

            b = 0
            if i==0:
                sigma_before_by2 = 0
                sigma_before_by3 = 0
            else:
                sigma_before_by2 = all_sigma[:,:,-2]
                sigma_before_by2 = np.diag(sigma_before_by2).sum()
                sigma_before_by3 = 0

            if i>1:
                sigma_before_by3 = all_sigma[:,:,-3]
                sigma_before_by3 = np.diag(sigma_before_by3).sum()

            sigma_now = np.diag(sigma).sum() #get stress after step
    
    
    
            all_sigma = np.dstack((all_sigma, sigma))
    
            all_v = np.dstack((all_v, new_v[0:np.shape(x)[0]])) #Full list of velocities
            print(i, sigma_now)
            i = i+1
    
            if equil==True:
                a = a+1
                #print('hi')
    
            if abs(sigma_now - sigma_before) < 1e-6: #Equilibrium if sigma is constant
                print('equilibrium')
                b = i #When you reached equilibrium
                if equil == False:
                    equil = True
                    step_func = position_update_stress
                    #new_v = new_v + gamma*image_copy(new_x)
                    dt = dt/10
            # elif abs(sigma_now - sigma_before_by2) < 1e-6:
            #     print('equilibrium_osc')
            #     b = i #When you reached equilibrium
            #     if equil == False:
            #         equil = True
            #         step_func = position_update_stress
            #         #new_v = new_v + gamma*image_copy(new_x)
            # elif abs(sigma_now - sigma_before_by3) < 1e-7:
            #     print('equilibrium_osc')
            #     b = i #When you reached equilibrium
            #     if equil == False:
            #         equil = True
            #         step_func = position_update_stress
            #         #new_v = new_v + gamma*image_copy(new_x)

        except KeyboardInterrupt:
            kinetic = all_v*all_v
            kinetic = kinetic.sum(axis = 1)
            kinetic = kinetic.sum(axis = 0)
            kinetic = kinetic*0.5
            energy = kinetic[1:] + all_U
            return time_x,  all_sigma,  b, energy


    kinetic = all_v*all_v
    kinetic = kinetic.sum(axis = 1)
    kinetic = kinetic.sum(axis = 0)
    kinetic = kinetic*0.5
    z = len(all_U)
    energy = kinetic[:z] + all_U
    all_sigma = all_sigma[:, :, 1:] #Remove initial 0's.
    return time_x, all_sigma, b, energy #Return positions of initial particles at every point every t
                                                          #Return velocities of initial particles at every point every t
                                                          #Return stress tensor at every point every t
                                                          #Return full ghost positions
                                                          #Return full ghost velocities. These are needed if you want to start the system from a specific point.



#%%
import random
r=3
length = 7.2
width = 7.2
height = 7.2
grid = disk.Grid(r, length, width, height)

rand = (random.uniform(0, length), random.uniform(0, width), random.uniform(0, height))
#Seed creates repeatable points
data = grid.poisson(rand, 30)
x = np.array(data)

x = x + length
dims = length/2

#%%
gamma = 0.1
dt = 0.1
results, sigma, equil_step, energy = many_step_equil(50, x, position_update, gamma, dt)
#To get diagonals
sigma1 = np.einsum('ijk->k', sigma)



#%%
t = np.arange(0, 1*len(sigma1), 1)
plt.plot(t, abs(sigma1))
plt.xlabel('time, reduced units')
plt.ylabel('Stress, reduced units')

#plt.xlim(49.5,50)


plt.show()
#%%
plt.plot(t[1:], energy)
plt.show()

#%%
np.savetxt('Verlet_results1.csv', sigma1)
#np.savetxt('after_equilibration6', results[:,:,-1])
#np.savetxt('after_equilibrationv6', v[:,:,-1])

