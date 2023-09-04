#%%
import numpy as np
import matplotlib.pyplot as plt
#%%

data = np.loadtxt('final_results4.csv')
#%%

dt = 0.1
time = np.arange(0, dt*len(data), dt)

length = 7

strain = 0.01*time
strain = strain/length
plt.plot(time, data)
#%%

#Step 823, volume fraction 0.8, gives stress strain
plt.plot(time[828:], data[828:],'.')

#%%
#Need linear stress region. 5 points.
plt.plot(time[827:833], data[827:833],'.')


#%%
#Shift


fit, cov = np.polyfit(strain[828:834], data[828:834], deg=1, cov = True)

parameters = np.poly1d(fit)
print(parameters)

plt.figure(figsize =[12,8])

plt.rc('font', size=10) #controls default text size
plt.rc('axes', titlesize=20) #fontsize of the title
plt.rc('axes', labelsize=20) #fontsize of the x and y labels
plt.rc('xtick', labelsize=15) #fontsize of the x tick labels
plt.rc('ytick', labelsize=15) #fontsize of the y tick labels
plt.rc('legend', fontsize=20) #fontsize of the legend

plt.plot(strain[828:834], parameters(strain[828:834]), label = 'fit')
plt.plot(strain[828:834], data[828:834], label = 'simulation')
plt.xlabel('strain')
plt.ylabel('stress, reduced units')
plt.legend()

plt.show()



print(np.sqrt(cov[0,0]))
#%%
kb = 1.38e-23
T = 300
zero_u = 5e-9
factor = kb*T/(zero_u**3)

mod = parameters[1]*kb*T/(zero_u**3)
print(mod, np.sqrt(cov[0,0])*factor)


#%%


def plot_fit(fname, first_slice, second_slice, length, dt, step_in_slice = 1, skip_t = False):
    #fname, then put in slices between where you want to fit i.e. where you applied strain, then put in the length
    data = np.loadtxt(fname)
    #dt = 0.01
    time = np.arange(0, dt*len(data), dt)
    #time = time[0:-1]

    if skip_t == True:
        time = time[0:-1]

    strain = 0.01*time
    strain = strain/length
    plt.plot(time, data)
    fit, cov = np.polyfit(strain[first_slice:second_slice:step_in_slice], data[first_slice:second_slice:step_in_slice], deg=1, cov = True)

    parameters = np.poly1d(fit)
    print(parameters)

    plt.figure(figsize =[12,8])

    plt.rc('font', size=10) #controls default text size
    plt.rc('axes', titlesize=20) #fontsize of the title
    plt.rc('axes', labelsize=20) #fontsize of the x and y labels
    plt.rc('xtick', labelsize=15) #fontsize of the x tick labels
    plt.rc('ytick', labelsize=15) #fontsize of the y tick labels
    plt.rc('legend', fontsize=20) #fontsize of the legend

    plt.plot(strain[first_slice:second_slice], parameters(strain[first_slice:second_slice]), label = 'fit')
    plt.plot(strain[first_slice:second_slice], data[first_slice:second_slice], label = 'simulation')
    plt.xlabel('strain')
    plt.ylabel('stress, reduced units')
    plt.legend()

    plt.show()



    print(np.sqrt(cov[0,0]))

    kb = 1.38e-23
    T = 300
    zero_u = 5e-9
    factor = kb*T/(zero_u**3)

    mod = parameters[1]*kb*T/(zero_u**3)
    print(mod, np.sqrt(cov[0,0])*factor)
#%%
data = np.loadtxt('final_results5.csv')
#%%
#Other strain data at 0.01
dt = 0.01
time = np.arange(0, dt*len(data), dt)

length = 7

strain = 0.01*time
strain = strain
plt.plot(time, data)
#%%

#Step 823, volume fraction 0.8, gives stress strain
plt.plot(time[1336:1346], data[1336:1346],'.')

#%%
#Can try fit to that
#Note dt has been set to 0.01
plot_fit('final_results5.csv', 1336,1346, length = 7.2, dt = 0.01)
#%%
#0.69
data = np.loadtxt('final_results6.5.csv')

#Other strain data at 0.01
dt = 0.01
time = np.arange(0, dt*len(data), dt)

length = 7

strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.xlabel('time, reduced units')
plt.ylabel('Stress, reduced units')
plt.title('Stress against time, simulation')
plt.show()
#%%


plt.plot(time[1562:1572], data[1562:1572],'.')
#%%
plot_fit('final_results6.5.csv', 1562,1572, length = 7.4, dt = 0.1)
#%%
data = np.loadtxt('final_results9.csv')

#Other strain data at 0.01
dt = 0.01
time = np.arange(0, dt*len(data), dt)

length = 7

strain = 0.01*time
strain = strain
plt.plot(time, data)
#%%

#Step 823, volume fraction 0.8, gives stress strain
plt.plot(time[825:], data[825:],'.')
#%%
plot_fit('final_results9.csv', 825, -1, length = 6.5, dt = 0.01)
#%%
#0.5
data = np.loadtxt('final_results7.csv')

#Other strain data at 0.01
dt = 0.01
time = np.arange(0, dt*len(data), dt)

length = 8

strain = 0.01*time
strain = strain
plt.plot(time, data)
#%%

#Step 823, volume fraction 0.8, gives stress strain
plt.plot(time[756:], data[756:],'.')
#%%
plot_fit('final_results7.csv',756, -1, length = 8, dt = 0.01)
#%%
#0.6 repeat
data = np.loadtxt('final_results10.csv')

#Other strain data at 0.01
dt = 0.01
time = np.arange(0, dt*len(data), dt)

length = 7.4

strain = 0.01*time
strain = strain
plt.plot(time, data)
#%%

#Step 823, volume fraction 0.8, gives stress strain
plt.plot(time[1311:], data[1311:],'.')
#%%
plot_fit('final_results10.csv',1311, -1, length = 7.4, dt = 0.01)
#%%
#0.4 repeat
data = np.loadtxt('final_results8.csv')

#Other strain data at 0.01
dt = 0.01
time = np.arange(0, dt*len(data), dt)

length = 8

strain = 0.01*time
strain = strain
plt.plot(time, data)
#%%


plt.plot(time[1775:], data[1775:],'.')
#%%
plot_fit('final_results8.csv',1775, -1, length = 8.5, dt = 0.01)
#%%
#0.4 repeat
data = np.loadtxt('final_results11.csv')

#Other strain data at 0.01
dt = 0.01
time = np.arange(0, dt*len(data), dt)

length = 8

strain = 0.01*time
strain = strain
plt.plot(time, data)
#%%


plt.plot(time[1186:1300], data[1186:1300],'.')
#%%
plot_fit('final_results11.csv',1186, 1300, length = 7.1, dt = 0.01)
#%%
#0.4 repeat
data = np.loadtxt('final_results12.csv')

#Other strain data at 0.01
dt = 0.01
time = np.arange(0, dt*len(data), dt)

length = 7

strain = 0.01*time
strain = strain
plt.plot(time, data)
#%%


plt.plot(time[2328:], data[2328:],'.')
#%%
plot_fit('final_results12.csv',2328, -1, length = 7.2, dt = 0.01)
#%%
data = np.loadtxt('final_results4.csv')

#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)

length = 7

strain = 0.01*time
strain = strain
plt.rc('font', size=12) #controls default text size
plt.rc('axes', titlesize=12) #fontsize of the title
plt.rc('axes', labelsize=12) #fontsize of the x and y labels
plt.rc('xtick', labelsize=12) #fontsize of the x tick labels
plt.rc('ytick', labelsize=12) #fontsize of the y tick labels
plt.rc('legend', fontsize=12) #fontsize of the legend
plt.plot(time, data)
plt.xlabel('Time, reduced units')
plt.ylabel('Stress, reduced units')
plt.title('Stress strain simulation with Euler, energy breakdown')
plt.show()
#%%


plt.plot(time[1176:1181], data[1176:1181],'.')
#%%
#dt is now 0.1
plot_fit('final_results13.csv',1176,1181, length = 7.0, dt = 0.1)

#%%
data = np.loadtxt('final_results14.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)

length = 7

strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()
#%%


plt.plot(time[1684:1689], data[1684:1689],'.')
#%%
#dt is now 0.1
plot_fit('final_results14.csv',1684,1690, length = 7.3, dt = 0.1)
#%%
data = np.loadtxt('final_results16.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)

length = 7

strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()
#%%


plt.plot(time[3288:3293], data[3288:3293],'.')
#%%
#dt is now 0.1
plot_fit('final_results16.csv',3288, 3293, length = 7.3, dt = 0.1)

#%%
data = np.loadtxt('final_results17.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)

length = 7

strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()
#%%


plt.plot(time[1393:1400], data[1393:1400],'.')
#%%
#dt is now 0.1
plot_fit('final_results17.csv',1393, 1400, length = 7.3, dt = 0.1)
#%%
data = np.loadtxt('final_results18.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)

length = 6.9

strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()
#%%


plt.plot(time[1508:1518:2], data[1508:1518:2],'.')
#%%

#new_time = time[1508:1518:2]
#np.savetxt('Cutdata0', X)

#%%
#dt is now 0.1
plot_fit('final_results18.csv',1508, 1518, length = 6.9, dt = 0.1, step_in_slice=2)
#%%

data = np.loadtxt('final_results21.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)

length = 7.5

strain = 0.01*time
strain = strain
plt.plot(time[:-1], data)
plt.show()
#%%


plt.plot(time[1339:1350:2], data[1339:1350:2],'.')

#%%
#dt is now 0.1
plot_fit('final_results21.csv',1339, 1350, length = 7.5, dt = 0.1, step_in_slice=2)

#%%

data = np.loadtxt('final_results22.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)

length = 7.5

strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()
#%%


plt.plot(time[676:686:2], data[676:686:2],'.')

#%%
#dt is now 0.1
plot_fit('final_results22.csv', 676, 686, length = 7.5, dt = 0.1, step_in_slice=2)
#%%

data = np.loadtxt('final_results23.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)

length = 7.5

strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()

#%%


plt.plot(time[1689:1693], data[1689:1693],'.')

#%%
#dt is now 0.1
plot_fit('final_results23.csv', 1689, 1693, length = 6.6, dt = 0.1)

#%%

data = np.loadtxt('final_results23.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)

length = 7.5

strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()

#%%


plt.plot(time[1689:1693], data[1689:1693],'.')

#%%
#dt is now 0.1
plot_fit('final_results23.csv', 1689, 1693, length = 6.6, dt = 0.1)

#%%

data = np.loadtxt('final_results25.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)

length = 6.8

strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()

#%%


plt.plot(time[1648:1656:2], data[1648:1656:2],'.')

#%%
#dt is now 0.1
plot_fit('final_results25.csv', 1648, 1656, length = 6.8, dt = 0.1, step_in_slice=2)
#%%
data = np.loadtxt('final_results26.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)

length = 6.6

strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()

#%%


plt.plot(time[790:794], data[790:794],'.')

#%%
#dt is now 0.1
plot_fit('final_results26.csv', 790, 794, length = 6.6, dt = 0.1, step_in_slice=1)
#%%
data = np.loadtxt('final_results27.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)

length = 6.6

strain = 0.01*time
strain = strain
plt.plot(time[0:-1], data)
plt.show()

#%%


plt.plot(time[1405:1418], data[1405:1418],'.')

#%%
#dt is now 0.1
plot_fit('final_results27.csv', 1405, 1418, length = 6.6, dt = 0.1, step_in_slice=2, skip_t = True)

#%%
data = np.loadtxt('final_results28.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)
time = time


strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()

#%%


plt.plot(time[1097:1101], data[1097:1101],'.')

#%%
#dt is now 0.1
plot_fit('final_results28.csv', 1097, 1101, length = 6.6, dt = 0.1)

#%%
data = np.loadtxt('final_results29.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)
time = time


strain = 0.01*time
strain = strain
plt.plot(time[0:-1], data)
plt.show()

#%%


plt.plot(time[1435:1445], data[1435:1445],'.')

#%%
#dt is now 0.1
plot_fit('final_results29.csv', 1435, 1445, length = 7.7, dt = 0.1, step_in_slice = 2, skip_t = True)

#%%
data = np.loadtxt('final_results30.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)
time = time


strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()

#%%


plt.plot(time[1381:1387], data[1381:1387],'.')

#%%
#dt is now 0.1
plot_fit('final_results30.csv', 1381, 1387, length = 7.7, dt = 0.1, step_in_slice = 2, skip_t = False)

#%%
data = np.loadtxt('final_results31.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)
time = time


strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()

#%%


plt.plot(time[1617:1624], data[1617:1624],'.')

#%%
#dt is now 0.1
plot_fit('final_results31.csv', 1617, 1624, length = 8, dt = 0.1, step_in_slice = 2, skip_t = False)
#%%
data = np.loadtxt('final_results33.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)
time = time


strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()

#%%


plt.plot(time[977:982], data[977:982],'.')

#%%
#dt is now 0.1
plot_fit('final_results33.csv', 977, 982, length = 9, dt = 0.1, step_in_slice = 1, skip_t = False)

#%%
data = np.loadtxt('final_results34.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)
time = time


strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()

#%%


plt.plot(time[1115:1121], data[1115:1121],'.')

#%%
#dt is now 0.1
plot_fit('final_results34.csv', 1115, 1121, length = 9, dt = 0.1, step_in_slice = 1, skip_t = False)
#%%

#%%
data = np.loadtxt('final_results35.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)
time = time


strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()

#%%


plt.plot(time[1115:1121], data[1115:1121],'.')

#%%
#dt is now 0.1
plot_fit('final_results34.csv', 1115, 1121, length = 9, dt = 0.1, step_in_slice = 1, skip_t = False)
#%%
data = np.loadtxt('final_results36.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)
time = time


strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()

#%%


plt.plot(time[1058:1070], data[1058:1070],'.')

#%%
#dt is now 0.1
plot_fit('final_results36.csv', 1058, 1070, length = 10, dt = 0.1, step_in_slice = 2, skip_t = False)

#%%
data = np.loadtxt('final_results37.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)
time = time


strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()

#%%


plt.plot(time[2560:2568], data[2560:2568],'.')

#%%
#dt is now 0.1
plot_fit('final_results37.csv', 2560, 2568, length = 10, dt = 0.1, step_in_slice = 1, skip_t = False)
#%%
data = np.loadtxt('final_results38.csv')
plt.figure()
#Other strain data at 0.01
dt = 0.1
time = np.arange(0, dt*len(data), dt)
time = time


strain = 0.01*time
strain = strain
plt.plot(time, data)
plt.show()

#%%


plt.plot(time[2485:2495], data[2485:2495],'.')

#%%
#dt is now 0.1
plot_fit('final_results38.csv', 2485, 2495, length = 10, dt = 0.1, step_in_slice = 1, skip_t = False)



#%%
#Don't use really big values.


x = np.array([0.8, 0.93, 0.51, 0.41, 0.67])
y = np.array([15687345.518498838, 21254539.762724902, 3542257.3889302723, 2135931.542880334, 10706994.507800948])
yerr = np.array([1925162.5699594666, 135365.08775883666, 1229.9887187284721, 992.772689827703, 1568284.5103769186])

plt.figure(figsize = [12,8])
plt.errorbar(x, y, yerr = yerr, fmt = 'x')
plt.yscale('log')
plt.ylim(1e6, 1e8)
plt.xlabel('Volume Fraction')
plt.ylabel("Young's Modulus (Pa)")
plt.grid()
plt.tight_layout()
plt.minorticks_on()
plt.show()

#%%
#Don't use really big values.
#Everything

x = np.array([0.8,0.69,0.62,0.6,0.93,0.5,0.41,0.667,0.72,0.67,0.55,0.57,0.75,0.75,0.57,0.565,0.86,0.746,0.83,0.81,0.87,0.53,0.5,0.466, 0.23, 0.23, 0.32, 0.33])
y = np.array([15687345.518498838,3368930.876375354,17748952.96679625,43576887.90293616,21254539.762724902,3542257.3889302723,2135931.542880334,54625787.66535763,45375579.853711955,10706994.507800948,26343464.274849147,18906195.903194763,21220746,15850547.177105863,18810203.025445092,32374830.123484746,18500648.98572915,20676542.12783383,14239467.972406685,9027952.782596895,16818849.038759362,7457919.276435808,23542396.215326816,8930638.015828501,5801914.641785345, 5632961.767778324,20312326.789456934,5632961.767778324])/3

yerr = np.array([1925162.5699594666,

 615.722269198783,

 6047539.130952778,

 342894.040363888,

 135365.08775883666,

 1229.9887187284721,

 992.772689827703,

 407800.30524966045,

 504802.5587395366,

 1568284.5103769186,

 2141653.124694977,

 2139781.234095671,

 4913783.64781045,

 1568206.0435185134,

 1312997.7610932211,

 2551013.707705302,

 2572670.352133408,

 2525679.7276102635,

 2734599.9770028056,

 1170101.5965196749,

 2497647.0459891185,

 841186.6803182963,

 3924906.188237604,

 1475284.5346049394,

340981.96569802787,346723.12212560105,1108314.8005955187,1737259.3829056236
])/3

plt.figure(figsize = [12,8])
plt.errorbar(x, y, yerr = yerr, fmt = 'x')
plt.yscale('log')
plt.ylim(1e6, 3e7)
plt.xlabel('Volume Fraction')
plt.ylabel("Shear Modulus (Pa)")
plt.title('Simulation of Modulus against volume fraction')
plt.grid()
plt.tight_layout()
plt.minorticks_on()
plt.show()
