#%%
import numpy as np
import matplotlib.pyplot as plt

def Mie(e, sigma, d, n=6, m=3):
    factor = (n/(n-m))*(n/m)**(m/(n-m))
    main = (sigma/d)**n - (sigma/d)**m
    return factor*e*main

def diff(x, fx):
    Delta_fx = fx[1:]-fx[0:-1]
    dfdx = Delta_fx/np.diff(x)# Note that the returned array is 1 element shorter than the original array
    #print(dfdx[7916])
    return dfdx

def square_diff(x, fx):
    #Cuts off first and last point
    a = diff(x, fx)
    b = diff(x[0:-1], a)

    return b

k = 1.38e-23
def modulus(U, N, T, phi_cp, d, d_equil, j=11, k = k):
    A = phi_cp*j/(5*np.pi*d_equil*d_equil)
    #Need to cut off first point to work with square diff
    term1 = 4*diff(d, U)
    term2 = d_equil*square_diff(d, U)
    G = N*k*T + A*(term1[0:-1] + term2)
    return G

#Returns index number of closest value to element.
def close(array, value):
    distance = array - value
    distance = np.abs(distance)
    return np.argmin(distance)

#%%
#Calculate R_hs
import scipy.integrate as si
T = 298

def hs(d, Req = 5e-9):
    k = 1.38e-23
    T = 298
    v = -Mie(3*k*T, Req*2**(2/3), d)
    v = v/(k*T)
    exp = np.exp(v)
    f = 1 - exp
    return f

Req = 5e-9
Rhs = si.quad(hs, 0, 2*4e-9, args = Req)
print(Rhs)

#%%

ds = np.linspace(1e-10, 0.8e-8, 1000)
res = hs(ds)
plt.plot(ds, hs(ds))
#%%
#Solving for R_eq depending on R_hs (i.e. rootfinding)

def root_func(R_eq, Rhs):
    #d = np.linspace(1e-12, 1e-7, 10000)
    f = 2*Rhs - si.quad(hs, 0, 2*Req, args = (R_eq))[0]
    return f



import scipy.optimize as op
roots = op.root(root_func, x0 = 1, args = (10e-9))
roots = roots.copy()
R_eq = roots['x']
print(R_eq)

def R_eq_calc(Rhs):
    roots = op.root(root_func, x0 = 1e-9, args = (Rhs))
    roots = roots.copy()
    R_eq = roots['x']
    return R_eq/2
#%%
#d = np.linspace(1e-12, 1e-7, 10000)
R_eq_calc(5e-9)

#%%
#Vary R_hs
T = 298
k = 1.38e-23
e = 3*k*T
#R_eq = 10e-9
#sigma = 2*R_eq/(2**(1/3))

#sigma = 9e-9
d = np.linspace(0.1e-9, 30e-9, 30000)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = 0.6
phi = 0.7
Rhs = np.linspace(1e-10, 10e-9, 1000)

N_assemb = N/6

G = np.array([])

#Sigma variation
for i in Rhs:




    #Minimum energy
    #implementing y/n close packing
    if phi < 0.6:
        R_eq = abs(R_eq_calc(i))
        #R_eq = 4e-9
        d_equil = 2*R_eq
    else:
        R_eq = abs(R_eq_calc(i))
        #R_eq = 4e-9
        d_equil = 2*R_eq*(phi_cp/phi)**(1/3)
        #print(d_equil)

    print(d_equil)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]
    print(d_min)

    U = Mie(e=e, sigma=(R_eq)*2**(2/3), d=d)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
    G = np.append(G, G_d)

plt.figure(figsize = [15, 10])
plt.plot(Rhs, G, label = 'model (2)')
#plt.scatter(0.5, G[8], c = 'red', label = "phi = 0.5. Values below aren't suitable for equation")

plt.grid()
plt.minorticks_on()
plt.xlabel('Volume Fraction')
plt.ylabel('Shear Modulus')
#plt.ylim(0, 3e4)
plt.xlabel('Rhs')
plt.ylabel('G, Shear Modulus (Pa)')
plt.title('Shear Modulus against Volume Fraction, Prediction')
plt.legend()
#plt.show()

phi_dep = diff(Rhs, G)
print('Rhs dependence', phi_dep[500]) #need 5nm

#%%
#Vary sigma
T = 298
k = 1.38e-23
e = 3*k*T
#R_eq = 10e-9
#sigma = 2*R_eq/(2**(1/3))

#sigma = 9e-9
d = np.linspace(0.01e-9, 3e-9, 500)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = 0.6
phi = 0.8


G = np.array([])

#Sigma variation
for i in range(5, 13, 1):



    #Changing sigma
    #Minimum energy

    R_eq = (i*1e-9)*2**(1/3)
    d_equil = 2*R_eq*(phi_cp/phi)**(1/3)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]

    U = Mie(e=e, sigma=i*1e-9, d=d)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
    G = np.append(G, G_d)

plt.figure(figsize = [8, 6])
sigma = np.arange(5, 13, 1)
plt.plot(sigma*1e-9, G,'.')
plt.xlabel('sigma')
plt.grid()
plt.ylabel('Shear Modulus, Pa')
plt.title('Relationship between shear modulus and sigma')

#%%
U = Mie(e=e, sigma=9e-9, d=d, n = 4, m = 2)
plt.plot(d, U)
plt.scatter(d_min, U[min_arg])
plt.xlim(6e-9, 30e-9)
plt.ylim(-3*e, 10*e)
plt.grid()
plt.minorticks_on()
plt.hlines(0, 0, 20, color = 'black')
plt.xlabel('Centre to centre separation')
plt.ylabel('Potential Energy')
plt.show()

#R_eq = np.linspace(3.5e-9, 7.5e-9, 50)

#%%
#Vary phi
T = 298
k = 1.38e-23
e = 3*k*T
#R_eq = 10e-9
#sigma = 2*R_eq/(2**(1/3))

#sigma = 9e-9
d = np.linspace(0.1e-9, 50e-9, 50000)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = 0.6
phi = np.arange(0.1, 1, 0.0005)
m = 9
n = 10




N_assemb = N/6

G = np.array([])

#Sigma variation
for i in phi:




    #Minimum energy
    #implementing y/n close packing
    if i < 0.6:
        R_eq = 3.8223e-9
        d_equil = 2*R_eq
    else:
        R_eq = 3.8223e-9
        d_equil = 2*R_eq*(phi_cp/i)**(1/3)
        #print(d_equil)

    #print(d_equil)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]
    #print(d_min)

    U = Mie(e=e, sigma=2*(3.5e-9)*(m/n)**(1/(n-m)), d=d, m = m, n = n)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
    G = np.append(G, G_d)


plt.figure(figsize = [15, 10])
plt.plot(phi, G, label = 'n = 10, m = 9')
#plt.scatter(0.5, G[8], c = 'red', label = "phi = 0.5. Values below aren't suitable for equation")

# T = 298
# k = 1.38e-23
# e = 3*k*T
# #R_eq = 10e-9
# #sigma = 2*R_eq/(2**(1/3))

# #sigma = 9e-9
# d = np.linspace(0.1e-9, 50e-9, 50000)
# bulk_density = 1200
# molar_mass = 150000
# N = 1200/(molar_mass*1.67e-27)
# phi_cp = 0.6
# phi = np.arange(0.1, 1, 0.0005)
# m = 3
# n = 6




# N_assemb = N/6

# G = np.array([])

# #Sigma variation
# for i in phi:




#     #Minimum energy
#     #implementing y/n close packing
#     if i < 0.6:
#         R_eq = 4e-9
#         d_equil = 2*R_eq
#     else:
#         R_eq = 4e-9
#         d_equil = 2*R_eq*(phi_cp/i)**(1/3)
#         #print(d_equil)

#     print(d_equil)
#     min_arg = close(d[:-2], d_equil)
#     d_min = d[min_arg]
#     print(d_min)

#     U = Mie(e=e, sigma=2*(3.5e-9)*(m/n)**(1/(n-m)), d=d, m = m, n = n)
#     #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

#     G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
#     G = np.append(G, G_d)


# plt.plot(phi, G, label = 'n = 6')


plt.grid()
plt.minorticks_on()
plt.xlabel('Volume Fraction')
plt.ylabel('Shear Modulus')
#plt.ylim(0, 3e4)
plt.xlabel(r'$\phi$, Volume Fraction')
plt.ylabel('G, Shear Modulus (Pa)')
plt.title('Shear Modulus against Volume Fraction, Prediction')
plt.legend()
#plt.show()

phi_dep = diff(phi, G)
print('phi dependence', phi_dep[1200]) #need 0.7

vol_frac = np.array([0.33348, 0.348, 0.348, 0.348, 0.5623, 0.43211])
caprino_frac = np.array([0.33348])
#feta_frac = np.array([0.348, 0.348, 0.348])
gouda_frac = np.array([0.5623])
mozz_frac = np.array([0.43211])
other_frac = np.array([0.48911, 0.46602, 0.49063, 0.59279, 0.64763, 0.35753, 0.37541, 0.37421])
ched_frac = 0.72
pec_frac = 0.67
par_frac = 0.77
more_frac = np.array([0.736, 0.629, 0.707])

capm = np.array([17.86756])*1000
#fetam = np.array([781.5, 503.5, 425.5])*1000/6
goudam = np.array([79.9])*1000
mozzm = np.array([26.1111])*1000
otherm_elastic = np.array([29.1, 27.6, 34.3, 30.1, 69.3, 9.3, 13.6, 10.5])
otherm_plastic = np.array([9.9, 10.3, 12.3, 11.2, 32.2, 3.3, 3.9, 3.0])
otherm = np.sqrt(otherm_elastic**2 + otherm_plastic**2)*1000
cheddm = 300*1000
pecm = (2840/3)*1000
parm = 706e3
morem = np.array([308, 326.6, 216])*1000


cape = np.array([1.414])*1000
#fetae = np.array([129.29, 79.52, 42.81])*1000/6
goudae = np.array([6.1])*1000
mozze = np.array([3.402])*1000
othere = 0.05*otherm
chedde = 10*1000
pece = 0.05*pecm
pare = 50
moree = np.array([25, 90, 0.05*216])*1000

#modulus = np.array([17.86756, 781.5, 503.5, 425.5, 79.9, 26.1111])
#error = np.array([1.414, 129.29, 79.52, 42.81, 6.1, 3.402])

#plt.figure(figsize = [12, 8])
plt.errorbar(caprino_frac, capm, yerr = cape, fmt = 'x')
#plt.errorbar(feta_frac, fetam, fetae, fmt = 'x')
plt.errorbar(gouda_frac, goudam, goudae, fmt = 'x')
plt.errorbar(mozz_frac, mozzm, mozze, fmt = 'x')
plt.errorbar(other_frac, otherm, othere, fmt = 'x')
plt.errorbar(ched_frac, cheddm, chedde, fmt = 'x')
plt.errorbar(pec_frac, pecm, pece, fmt = 'x')
plt.errorbar(par_frac, parm, pare, fmt = 'x')
plt.errorbar(more_frac, morem, moree, fmt = 'x')


plt.minorticks_on()
plt.xlabel('Volume Fraction')
plt.ylabel('Shear Modulus')
#plt.ylim(0, 3e4)
plt.xlabel(r'$\phi$, Volume Fraction')
plt.ylabel('G, Shear Modulus (Pa)')
plt.title('Shear Modulus against Volume Fraction, Comparing to $\phi$ > $\phi_{CP}$')
plt.legend()
#plt.xlim(0, 0.8)
#%%

#Vary phi
T = 298
k = 1.38e-23
e = 3*k*T
#R_eq = 10e-9
#sigma = 2*R_eq/(2**(1/3))

#sigma = 9e-9
d = np.linspace(0.1e-9, 50e-9, 50000)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = 0.6
phi = np.arange(0.1, 1, 0.0005)
m = 3
n = 6




N_assemb = N/6

G = np.array([])

#Sigma variation
for i in phi:




    #Minimum energy
    #implementing y/n close packing
    if i < 0.6:
        R_eq = 3.8223e-9
        d_equil = 2*R_eq
    else:
        R_eq = 3.8223e-9
        d_equil = 2*R_eq*(phi_cp/i)**(1/3)
        #print(d_equil)

    #print(d_equil)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]
    #print(d_min)

    U = Mie(e=e, sigma=2*(3.5e-9)*(m/n)**(1/(n-m)), d=d, m = m, n = n)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
    G = np.append(G, G_d)


plt.figure(figsize = [15, 10])
plt.plot(phi, G, label = 'n = 10, m = 9')
#plt.scatter(0.5, G[8], c = 'red', label = "phi = 0.5. Values below aren't suitable for equation")

# T = 298
# k = 1.38e-23
# e = 3*k*T
# #R_eq = 10e-9
# #sigma = 2*R_eq/(2**(1/3))

# #sigma = 9e-9
# d = np.linspace(0.1e-9, 50e-9, 50000)
# bulk_density = 1200
# molar_mass = 150000
# N = 1200/(molar_mass*1.67e-27)
# phi_cp = 0.6
# phi = np.arange(0.1, 1, 0.0005)
# m = 3
# n = 6




# N_assemb = N/6

# G = np.array([])

# #Sigma variation
# for i in phi:




#     #Minimum energy
#     #implementing y/n close packing
#     if i < 0.6:
#         R_eq = 4e-9
#         d_equil = 2*R_eq
#     else:
#         R_eq = 4e-9
#         d_equil = 2*R_eq*(phi_cp/i)**(1/3)
#         #print(d_equil)

#     print(d_equil)
#     min_arg = close(d[:-2], d_equil)
#     d_min = d[min_arg]
#     print(d_min)

#     U = Mie(e=e, sigma=2*(3.5e-9)*(m/n)**(1/(n-m)), d=d, m = m, n = n)
#     #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

#     G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
#     G = np.append(G, G_d)


# plt.plot(phi, G, label = 'n = 6')


plt.grid()
plt.minorticks_on()
plt.xlabel('Volume Fraction')
plt.ylabel('Shear Modulus')
#plt.ylim(0, 3e4)
plt.xlabel(r'$\phi$, Volume Fraction')
plt.ylabel('G, Shear Modulus (Pa)')
plt.title('Shear Modulus against Volume Fraction, Prediction')
plt.legend()
#plt.show()

phi_dep = diff(phi, G)
print('phi dependence', phi_dep[1200]) #need 0.7

vol_frac = np.array([0.33348, 0.348, 0.348, 0.348, 0.5623, 0.43211])
caprino_frac = np.array([0.33348])
#feta_frac = np.array([0.348, 0.348, 0.348])
gouda_frac = np.array([0.5623])
mozz_frac = np.array([0.43211])
other_frac = np.array([0.48911, 0.46602, 0.49063, 0.59279, 0.64763, 0.35753, 0.37541, 0.37421])
ched_frac = 0.72
pec_frac = 0.67
par_frac = 0.77
more_frac = np.array([0.736, 0.629, 0.707])

capm = np.array([17.86756])*1000
#fetam = np.array([781.5, 503.5, 425.5])*1000/6
goudam = np.array([79.9])*1000
mozzm = np.array([26.1111])*1000
otherm_elastic = np.array([29.1, 27.6, 34.3, 30.1, 69.3, 9.3, 13.6, 10.5])
otherm_plastic = np.array([9.9, 10.3, 12.3, 11.2, 32.2, 3.3, 3.9, 3.0])
otherm = np.sqrt(otherm_elastic**2 + otherm_plastic**2)*1000
cheddm = 300*1000
pecm = (2840/3)*1000
parm = 706e3
morem = np.array([308, 326.6, 216])*1000


cape = np.array([1.414])*1000
#fetae = np.array([129.29, 79.52, 42.81])*1000/6
goudae = np.array([6.1])*1000
mozze = np.array([3.402])*1000
othere = 0.05*otherm
chedde = 10*1000
pece = 0.05*pecm
pare = 50
moree = np.array([25, 90, 0.05*216])*1000

#modulus = np.array([17.86756, 781.5, 503.5, 425.5, 79.9, 26.1111])
#error = np.array([1.414, 129.29, 79.52, 42.81, 6.1, 3.402])

#plt.figure(figsize = [12, 8])
plt.errorbar(caprino_frac, capm, yerr = cape, fmt = 'x')
#plt.errorbar(feta_frac, fetam, fetae, fmt = 'x')
plt.errorbar(gouda_frac, goudam, goudae, fmt = 'x')
plt.errorbar(mozz_frac, mozzm, mozze, fmt = 'x')
plt.errorbar(other_frac, otherm, othere, fmt = 'x')
plt.errorbar(ched_frac, cheddm, chedde, fmt = 'x')
plt.errorbar(pec_frac, pecm, pece, fmt = 'x')
plt.errorbar(par_frac, parm, pare, fmt = 'x')
plt.errorbar(more_frac, morem, moree, fmt = 'x')


plt.minorticks_on()
plt.xlabel('Volume Fraction')
plt.ylabel('Shear Modulus')
#plt.ylim(0, 3e4)
plt.xlabel(r'$\phi$, Volume Fraction')
plt.ylabel('G, Shear Modulus (Pa)')
plt.title('Shear Modulus against Volume Fraction, Comparing to $\phi$ > $\phi_{CP}$')
plt.legend()
#plt.xlim(0, 0.8)
#%%
#Maybe try fitting?
#Vary phi

def phi_func(phi, n, m):
    n = np.round(n)
    m = np.round(m)

    T = 298
    k = 1.38e-23
    e = 3*k*T
    #R_eq = 10e-9
    #sigma = 2*R_eq/(2**(1/3))

    #sigma = 9e-9
    d = np.linspace(0.1e-9, 50e-9, 5000)
    bulk_density = 1200
    molar_mass = 150000
    N = 1200/(molar_mass*1.67e-27)
    phi_cp = 0.6
    #phi = np.arange(0.1, 1, 0.005)





    N_assemb = N/6

    G = np.array([])

    #Sigma variation
    for i in phi:




        #Minimum energy
        #implementing y/n close packing
        if i < 0.6:
            R_eq = 4e-9
            d_equil = 2*R_eq
        else:
            R_eq = 4e-9
            d_equil = 2*R_eq*(phi_cp/i)**(1/3)
            #print(d_equil)
    
        print(d_equil)
        min_arg = close(d[:-2], d_equil)
        d_min = d[min_arg]
        print(d_min)
    
        U = Mie(e=e, sigma=2*(3.5e-9)*(m/n)**(1/(n-m)), d=d, m = m, n = n)
        #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.
    
        G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
        G = np.append(G, G_d)

    return G

from scipy.optimize import curve_fit

frac_data =np.array([0.33348, 0.5623, 0.43211, 0.48911, 0.46602, 0.49063, 0.59279, 0.64763, 0.35753, 0.37541, 0.37421, 0.72, 0.67, 0.77, 0.736, 0.629, 0.707])
G_data = np.array([17.86756, 79.9, 26.1111])

otherm_elastic = np.array([29.1, 27.6, 34.3, 30.1, 69.3, 9.3, 13.6, 10.5])
otherm_plastic = np.array([9.9, 10.3, 12.3, 11.2, 32.2, 3.3, 3.9, 3.0])
otherm = np.sqrt(otherm_elastic**2 + otherm_plastic**2)
G_data = np.append(G_data, otherm)
G_data = np.append(G_data, [300, 2840/3, 706, 308, 326.6, 216])
G_data = G_data*1000

fit, cov = curve_fit(phi_func, frac_data, G_data, p0 = [6,3] )
print(fit)

plt.figure(figsize = [15, 10])
plt.plot(phi, G, label = 'model (2)')
#plt.scatter(0.5, G[8], c = 'red', label = "phi = 0.5. Values below aren't suitable for equation")

plt.grid()
plt.minorticks_on()
plt.xlabel('Volume Fraction')
plt.ylabel('Shear Modulus')
#plt.ylim(0, 3e4)
plt.xlabel(r'$\phi$, Volume Fraction')
plt.ylabel('G, Shear Modulus (Pa)')
plt.title('Shear Modulus against Volume Fraction, Prediction')
plt.legend()
#plt.show()

phi_dep = diff(phi, G)
print('phi dependence', phi_dep[1200]) #need 0.7

vol_frac = np.array([0.33348, 0.348, 0.348, 0.348, 0.5623, 0.43211])
caprino_frac = np.array([0.33348])
#feta_frac = np.array([0.348, 0.348, 0.348])
gouda_frac = np.array([0.5623])
mozz_frac = np.array([0.43211])
other_frac = np.array([0.48911, 0.46602, 0.49063, 0.59279, 0.64763, 0.35753, 0.37541, 0.37421])
ched_frac = 0.72
pec_frac = 0.67
par_frac = 0.77
more_frac = np.array([0.736, 0.629, 0.707])

capm = np.array([17.86756])*1000
#fetam = np.array([781.5, 503.5, 425.5])*1000/6
goudam = np.array([79.9])*1000
mozzm = np.array([26.1111])*1000
otherm_elastic = np.array([29.1, 27.6, 34.3, 30.1, 69.3, 9.3, 13.6, 10.5])
otherm_plastic = np.array([9.9, 10.3, 12.3, 11.2, 32.2, 3.3, 3.9, 3.0])
otherm = np.sqrt(otherm_elastic**2 + otherm_plastic**2)*1000
cheddm = 300*1000
pecm = (2840/3)*1000
parm = 706e3
morem = np.array([308, 326.6, 216])*1000


cape = np.array([1.414])*1000
#fetae = np.array([129.29, 79.52, 42.81])*1000/6
goudae = np.array([6.1])*1000
mozze = np.array([3.402])*1000
othere = 0.05*otherm
chedde = 10*1000
pece = 0.05*pecm
pare = 50
moree = np.array([25, 90, 0.05*216])*1000

#modulus = np.array([17.86756, 781.5, 503.5, 425.5, 79.9, 26.1111])
#error = np.array([1.414, 129.29, 79.52, 42.81, 6.1, 3.402])

#plt.figure(figsize = [12, 8])
plt.errorbar(caprino_frac, capm, yerr = cape, fmt = 'x')
#plt.errorbar(feta_frac, fetam, fetae, fmt = 'x')
plt.errorbar(gouda_frac, goudam, goudae, fmt = 'x')
plt.errorbar(mozz_frac, mozzm, mozze, fmt = 'x')
plt.errorbar(other_frac, otherm, othere, fmt = 'x')
plt.errorbar(ched_frac, cheddm, chedde, fmt = 'x')
plt.errorbar(pec_frac, pecm, pece, fmt = 'x')
plt.errorbar(par_frac, parm, pare, fmt = 'x')
plt.errorbar(more_frac, morem, moree, fmt = 'x')


plt.minorticks_on()
plt.xlabel('Volume Fraction')
plt.ylabel('Shear Modulus')
#plt.ylim(0, 3e4)
plt.xlabel(r'$\phi$, Volume Fraction')
plt.ylabel('G, Shear Modulus (Pa)')
plt.title('Shear Modulus against Volume Fraction, Comparing to $\phi$ > $\phi_{CP}$')
plt.legend()
plt.xlim(0, 0.8)
plt.show()



#%%
#Vary phi_cp!
T = 298
k = 1.38e-23
e = 3*k*T
#R_eq = 10e-9
#sigma = 2*R_eq/(2**(1/3))

#sigma = 9e-9
d = np.linspace(0.1e-9, 50e-9, 500)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = np.arange(0, 1, 0.0001)
phi = 0.7

N_assemb = N/6

G = np.array([])

#Sigma variation
for i in phi_cp:
    #Changing sigma
    #Minimum energy
    #implementing y/n close packing
    if i < 0.6:
        R_eq = 10e-9
        d_equil = 2*R_eq
    else:
        R_eq = 10e-9
        d_equil = 2*R_eq*(i/phi)**(1/3)

    print(d_equil)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]

    U = Mie(e=e, sigma=7e-9, d=d)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, i, d, d_equil)[min_arg]
    G = np.append(G, G_d)


phi_dep = diff(phi_cp, G)
print('phicp dependence', phi_dep[6000]) #need 0.7
#%%
#Sensitivity analysis
#Set phi = 0.7
T = 298
k = 1.38e-23
e = 3*k*T
#R_eq = 10e-9
#sigma = 2*R_eq/(2**(1/3))

#sigma = 9e-9
d = np.linspace(0.1e-9, 50e-9, 5000)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = 0.6
phi = np.array([0.7])

N_assemb = N/6

G = np.array([])

#Sigma variation
for i in phi:



    #Changing sigma
    #Minimum energy
    #implementing y/n close packing
    if i < 0.6:
        R_eq = 4e-9
        d_equil = 2*R_eq
    else:
        R_eq = 4e-9
        d_equil = 2*R_eq*(phi_cp/i)**(1/3)

    print(d_equil)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]

    U = Mie(e=e, sigma=(3.5e-9)*2**(2/3), d=d)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
    G = np.append(G, G_d)

print('j dependence', (G - N*k*T)/11)
#print('phi_cp dependence', (G - N*k*T)/0.6)
print('e dependence', (G - N*k*T)/(3*k*T))
#%%
#Vary R_eq
T = 298
k = 1.38e-23
e = 3*k*T
#R_eq = 10e-9
#sigma = 2*R_eq/(2**(1/3))

#sigma = 9e-9
d = np.linspace(0.1e-9, 50e-9, 500)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = 0.6
R = np.arange(1e-9, 6e-9, 0.5e-11)
phi = 0.7
N_assemb = N/6

G = np.array([])

#Sigma variation
for i in R:
    #Changing sigma
    #Minimum energy
    #implementing y/n close packing
    if phi < 0.6:
        R_eq = i
        d_equil = 2*R_eq
    else:
        R_eq = i
        d_equil = 2*R_eq*(phi_cp/0.7)**(1/3)

    print(d_equil)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]

    U = Mie(e=e, sigma=i*2**(2/3), d=d)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
    G = np.append(G, G_d)

plt.plot(R, G)
R_dep = diff(R, G)
print('R dependence', R_dep[300]) #need 0.7
#%%
#Vary phi
T = 298
k = 1.38e-23
e = 3*k*T
#R_eq = 10e-9
#sigma = 2*R_eq/(2**(1/3))

#sigma = 9e-9
d = np.linspace(0.1e-9, 50e-9, 500)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = 0.6
phi = np.arange(0.1, 1, 0.0005)

N_assemb = N/6

G = np.array([])

#Sigma variation
for i in phi:



    #Changing sigma
    #Minimum energy
    #implementing y/n close packing
    if i < 0.6:
        R_eq = 10e-9
        d_equil = 2*R_eq
    else:
        R_eq = 10e-9
        d_equil = 2*R_eq*(phi_cp/i)**(1/3)

    print(d_equil)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]

    U = Mie(e=e, sigma=7e-9, d=d, m = 1)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
    G = np.append(G, G_d)

plt.figure(figsize = [12, 8])
plt.plot(phi, G,'.')
#plt.scatter(0.5, G[8], c = 'red', label = "phi = 0.5. Values below aren't suitable for equation")
plt.grid()
plt.minorticks_on()
plt.xlabel('Volume Fraction')
plt.ylabel('Shear Modulus')
#plt.ylim(0, 3e4)
plt.xlabel(r'$\phi$, Volume Fraction')
plt.ylabel('G, Shear Modulus (Pa)')
plt.title('Shear Modulus against Volume Fraction, Prediction')
plt.legend()
plt.show()

phi_dep = diff(phi, G)
print('phi dependence', phi_dep[1200]) #need 0.7
#%%
#n dependence? Differentiate G.
#Sensitivity analysis
#Set phi = 0.7
T = 298
k = 1.38e-23
e = 3*k*T
#R_eq = 10e-9
#sigma = 2*R_eq/(2**(1/3))

#sigma = 9e-9
d = np.linspace(0.1e-9, 50e-9, 500)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = 0.6
n = np.arange(1, 10, 0.001)
phi = 0.4
N_assemb = N/6

G = np.array([])

#Sigma variation
for i in n:



    #Changing sigma
    #Minimum energy
    #implementing y/n close packing
    if i < 0.6:
        R_eq = 4e-9
        d_equil = 2*R_eq
    else:
        R_eq = 4e-9
        d_equil = 2*R_eq*(phi_cp/phi)**(1/3)

    print(d_equil)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]

    U = Mie(e=e, sigma=2*R_eq*(3/i)**(1/(i-3)), d=d, n = i)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
    G = np.append(G, G_d)

n_dep = diff(n, G)
print('n dependence', n_dep[5000])
plt.plot(n, G)
#%%

#m dependence? Differentiate G.
#Sensitivity analysis
#Set phi = 0.7
T = 298
k = 1.38e-23
e = 3*k*T
#R_eq = 10e-9
#sigma = 2*R_eq/(2**(1/3))

#sigma = 9e-9
d = np.linspace(0.1e-9, 50e-9, 5000)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = 0.6
m = np.arange(1, 10, 0.001)
phi = 0.4

N_assemb = N/6

G = np.array([])

#Sigma variation
for i in m:



    #Changing sigma
    #Minimum energy
    #implementing y/n close packing
    if i < 0.6:
        R_eq = 4e-9
        d_equil = 2*R_eq
    else:
        R_eq = 4e-9
        d_equil = 2*R_eq*(phi_cp/phi)**(1/3)

    print(d_equil)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]

    U = Mie(e=e, sigma=2*R_eq*(i/6)**(1/(6-i)), d=d, m = i)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
    G = np.append(G, G_d)

m_dep = diff(m, G)
print('m dependence', m_dep[2000])
plt.plot(m, G)

#%%
#Vary phi_cp
T = 298
k = 1.38e-23
e = 3*k*T


d = np.linspace(0.1e-9, 30e-9, 500)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = np.arange(0.1, 1, 0.05)
phi = 0.8


G = np.array([])

#Sigma variation
for i in phi_cp:



    #Changing sigma
    #Minimum energy

    R_eq = (9e-9)*2**(1/3)
    d_equil = 2*R_eq*(i/phi)**(1/3)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]

    U = Mie(e=e, sigma=9e-9, d=d)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, i, d, d_equil)[min_arg]
    G = np.append(G, G_d)


plt.plot(phi_cp, G,'.')
#plt.scatter(0.5, G[8], c = 'red', label = "phi = 0.5. Values below aren't suitable for equation")
plt.grid()
plt.minorticks_on()

plt.xlabel('Volume Fraction close packing')
plt.ylabel('Shear Modulus')
plt.legend()
plt.show()
#%%
#Vary N
T = 298
k = 1.38e-23
e = 3*k*T


d = np.linspace(0.1e-9, 30e-9, 500)
bulk_density = 1200
molar_mass = 150000
N = np.arange(0.1*1200/(molar_mass*1.67e-27),0.5*1500/(molar_mass*1.67e-27), 1e23)
phi_cp = 0.6
R_HS = 4.4e-9



G = np.array([])

#Sigma variation
for i in N:

    phi = (1/3)*i*4*np.pi*R_HS**3

    #Changing sigma
    #Minimum energy

    R_eq = (9e-9)*2**(1/3)
    d_equil = 2*R_eq*(phi_cp/phi)**(1/3)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]

    U = Mie(e=e, sigma=9e-9, d=d)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, i, T, phi_cp, d, d_equil)[min_arg]
    G = np.append(G, G_d)


plt.plot(N, G,'.')
#plt.scatter(0.5, G[8], c = 'red', label = "phi = 0.5. Values below aren't suitable for equation")
plt.grid()
plt.minorticks_on()

plt.xlabel('Number Density of assemblies')
plt.ylabel('Shear Modulus')
plt.legend()
plt.show()
#%%
#Vary R_HS
T = 298
k = 1.38e-23
e = 3*k*T


d = np.linspace(0.1e-9, 30e-9, 500)
bulk_density = 1200
molar_mass = 150000
N = (1/6)*1200/(molar_mass*1.67e-27)
phi_cp = 0.6
R_hs = np.arange(3.5e-9, 7e-9, 0.5e-9)


G = np.array([])

#Sigma variation
for i in R_hs:

    phi = (1/3)*N*4*np.pi*i**3


    #Changing sigma
    #Minimum energy

    R_eq = 10e-9
    d_equil = 2*R_eq*(phi_cp/phi)**(1/3)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]

    U = Mie(e=e, sigma=9e-9, d=d)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
    G = np.append(G, G_d)


plt.plot(R_hs, G,'.')
#plt.scatter(0.5, G[8], c = 'red', label = "phi = 0.5. Values below aren't suitable for equation")
plt.grid()
plt.minorticks_on()

plt.xlabel('Hard sphere radius')
plt.ylabel('Shear Modulus')
plt.legend()
plt.show()
#%%
#Vary m
T = 298
k = 1.38e-23
e = 3*k*T


d = np.linspace(0.1e-9, 30e-9, 500)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = 0.6
phi = 0.7
m = np.arange(1, 13)

G = np.array([])

#Sigma variation
for i in m:



    #Changing sigma
    #Minimum energy

    R_eq = 4e-9
    d_equil = 2*R_eq*(phi_cp/phi)**(1/3)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]

    U = Mie(e=e, sigma=2*R_eq*(i/6)**(1/(6-i)), d=d, m = i)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
    G = np.append(G, G_d)


plt.plot(m, G,'.')
#plt.scatter(0.5, G[8], c = 'red', label = "phi = 0.5. Values below aren't suitable for equation")
plt.grid()
plt.minorticks_on()

plt.xlabel('m, attraction')
plt.ylabel('Shear Modulus')
plt.legend()
plt.show()

#%%
#Vary n
#Vary R_HS
T = 298
k = 1.38e-23
e = 3*k*T


d = np.linspace(0.1e-9, 30e-9, 500)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = 0.6
phi = 0.8
n = np.arange(1, 13)

G = np.array([])

#Sigma variation
for i in n:



    #Changing sigma
    #Minimum energy

    R_eq = (9e-9)*2**(1/3)
    d_equil = 2*R_eq*(phi_cp/phi)**(1/3)
    min_arg = close(d[:-2], d_equil)
    d_min = d[min_arg]

    U = Mie(e=e, sigma=9e-9, d=d, n = i)
    #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

    G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]
    G = np.append(G, G_d)


plt.plot(n, G,'.')
#plt.scatter(0.5, G[8], c = 'red', label = "phi = 0.5. Values below aren't suitable for equation")
plt.grid()
plt.minorticks_on()

plt.xlabel('n, attraction')
plt.ylabel('Shear Modulus')
plt.legend()
plt.show()
#%%
#Vary m and n
T = 298
k = 1.38e-23
e = 3*k*T


d = np.linspace(0.1e-9, 30e-9, 500)
bulk_density = 1200
molar_mass = 150000
N = 1200/(molar_mass*1.67e-27)
phi_cp = 0.6
phi = 0.4
n = np.arange(1, 15)
m = np.arange(1, 15)
T = 298
k = 1.38e-23
e = 3*k*T
#R_eq = 10e-9
#sigma = 2*R_eq/(2**(1/3))


all_data = np.zeros(14)
#Sigma variation
for i in m:
    some_data = np.array([])
    for j in n:



       #Changing sigma
       #Minimum energy

        R_eq = 4e-9
        d_equil = 2*4e-9 #if phi < 0.6
        #d_equil = 2*R_eq*(phi_cp/phi)**(1/3)
        min_arg = close(d[:-2], d_equil)
        d_min = d[min_arg]

        U = Mie(e=e, sigma=2*(3.5e-9)*(i/j)**(1/(j-i)), d=d, n =j, m = i)
        #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

        G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]

        #data = np.array([G_d])
        some_data = np.append(some_data, G_d)

    all_data = np.vstack([all_data, some_data])



all_data = all_data[1:]
#all_data = np.nan_to_num(all_data)

# from mpl_toolkits import mplot3d
 
# # creating an empty canvas
# fig = plt.figure()

 
# # defining the axes with the projection
# # as 3D so as to plot 3D graphs
# ax = plt.axes(projection="3d")
 

# ax.scatter3D(all_data[:,0], all_data[:,1], all_data[:,2], c=all_data[:,2], cmap='cividis');
# plt.xlabel('n')
# plt.ylabel('m')

# Showing the above plot
plt.show()

X, Y = np.meshgrid(m, n)

#levels = np.arange(600, 700, 1)

plt.figure(figsize=(20,10))
plt.contourf(X, Y, all_data, cmap='inferno')
plt.colorbar()
plt.xlabel('m')
plt.ylabel('n')
plt.title('Value of G at $\phi$ = 0.4 for different n, m')
#%%
#The 2 variable loop
def G_loop(sigma, phi_cp, phi, e, d, d_equil, n, m, N, T):
        R_eq = (sigma)*2**(1/3)
        #d_equil = 2*R_eq*(phi_cp/phi)**(1/3)
        min_arg = close(d[:-2], d_equil)
        d_min = d[min_arg]

        U = Mie(e=e, sigma=9e-9, d=d, n =j, m = i)

        #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

        G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]

        data = np.array([i, j, G_d])

        return data
#%%
#
def G_loop_simple(sigma, phi_cp, phi, e, d, d_equil, n, m, N, T):
        R_eq = (sigma)*2**(1/3)
        #d_equil = 2*R_eq*(phi_cp/phi)**(1/3)
        min_arg = close(d[:-2], d_equil)
        d_min = d[min_arg]

        U = Mie(e=e, sigma=9e-9, d=d, n = n, m = m)

        #Need to pass through all d to differentiate, but actual value for G is at min - look at 2R_eq.

        G_d = modulus(U, N, T, phi_cp, d, d_equil)[min_arg]


        return G_d
#%%
#Vary m and d*
T = 298
k = 1.38e-23
e = 3*k*T


d = np.linspace(0.1e-9, 30e-9, 500)
bulk_density = 1200
molar_mass = 150000
N = (1/6)*1200/(molar_mass*1.67e-27)
phi_cp = 0.6
phi = 0.8
n = np.arange(1, 13)
m = 3
sigma = 9e-9
d_equil = np.arange(9e-9, 20e-9, 1e-9)

global all_data
all_data = np.array([0,0,0])


for i in n:
    for j in d_equil:
       data = G_loop(sigma, phi_cp, phi, e, d, j, i, m, N, T)
       all_data = np.vstack([all_data, data])


all_data = all_data[1:]

from mpl_toolkits import mplot3d
 
# creating an empty canvas
fig = plt.figure(figsize = (10, 9))
ax = plt.axes(projection="3d")
 

ax.scatter3D(all_data[:,0], all_data[:,1], all_data[:,2], c=all_data[:,2], cmap='cividis');
ax.set_xlabel('n')
ax.set_ylabel('d*')
ax.set_zlabel('G')
#%%
#Vary m and n
T = 298
k = 1.38e-23
e = 3*k*T


d = np.linspace(0.1e-9, 30e-9, 500)
bulk_density = 1200
molar_mass = 150000
N = (1/6)*1200/(molar_mass*1.67e-27)
phi_cp = 0.6
phi = 0.8
n = np.arange(4, 9)
m = 3

d_equil = np.arange(3e-9, 20e-9, 1e-9)





for i in n:
    G = np.array([])
    for j in d_equil:
        data = G_loop_simple(d_equil*(phi/2*phi_cp)**(1/3), phi_cp, phi, e, d, j, i, m, N, T)
        G = np.append(G, data)

    plt.plot(d_equil, G,'.')
    #plt.scatter(0.5, G[8], c = 'red', label = "phi = 0.5. Values below aren't suitable for equation")
    plt.grid()
    plt.minorticks_on()

    plt.xlabel('d_equil, n = {}'.format(i))
    plt.ylabel('Shear Modulus')
    plt.legend()
    #plt.ylim(0.001 + 3.28335e3, 0.011 +  3.28335e3)
    plt.show()
