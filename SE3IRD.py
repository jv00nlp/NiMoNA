

"""

 The MuensterVirus SIR network - SE3IRD - Model


"""

import numpy as np
import matplotlib.pyplot as plt
import csv

# DEFINE TWO NECESSARY FUNCTIONS

def rk4_step(rhs, x, function_parameters, h):
    k1 = rhs(x, *function_parameters)
    k2 = rhs(x + k1*h/2., *function_parameters)
    k3 = rhs(x + k2*h/2., *function_parameters)
    k4 = rhs(x + k3*h, *function_parameters)
    return x + h/6.*(k1 + 2*(k2 + k3) + k4)

def SIR(U, alpha, beta, gamma, age_distr, leth, rec, matrix):
    S, E, Iy, Im, Io, R, D = U
    N = np.zeros(len(matrix))
    
    for i in range(0,len(matrix)):
        N[i] = Sn[0,i] + En[0,i] + Iyn[0,i] + Imn[0,i] + Ion[0,i] + Rn[0,i]
   
    dS = - alpha*S*(Iy+Im+Io)/N
    dE = alpha*S*(Iy+Im+Io)/N - gamma*E
    dIy = gamma*age_dist[0]*E - beta*Iy
    dIm = gamma*age_dist[1]*E - beta*Im
    dIo = gamma*age_dist[2]*E - beta*Io
    dR =  beta * (Iy*rec[0] + Im*rec[1] + Io*rec[2])
    dD = beta * (Iy*leth[0] + Im*leth[1] + Io*leth[2])
    

  
    
    for n in range(0,len(matrix)):
            for m in range(0,len(matrix)):
                
                if m != n:
                
                    
                    dS[n] = dS[n] + matrix[n,m]/N[m] * S[m] - matrix[m,n]/N[n] * S[n]
                    dE[n] = dE[n] + matrix[n,m]/N[m] * E[m] - matrix[m,n]/N[n] * E[n]
                    dIy[n] = dIy[n] + matrix[n,m]/N[m] * Iy[m] - matrix[m,n]/N[n] * Iy[n]
                    dIm[n] = dIm[n] + matrix[n,m]/N[m] * Im[m] - matrix[m,n]/N[n] * Im[n]
                    dIo[n] = dIo[n] + matrix[n,m]/N[m] * Io[m] - matrix[m,n]/N[n] * Io[n]
                    dR[n] = dR[n] + matrix[n,m]/N[m] * R[m] - matrix[m,n]/N[n] * R[n]
                    dD[n] = dD[n] + matrix[n,m]/N[m] * D[m] - matrix[m,n]/N[n] * D[n]
                    
                  
                
    return np.array([dS, dE, dIy, dIm, dIo, dR, dD])



# PARAMETERS
    
dt = 0.1        # timestep
Tend = 150      # Number of days
Nt = int(Tend/dt)
ts = np.linspace(0, Tend, Nt)

# PANDEMIC PARAMETERS

alpha = 0.35    # infection rate    (from S to E)
beta = 0.035    # removal rate      (from I to R/D)
gamma = 0.3     # transmission rate (from E to I)

age_dist = [0.18, 0.62, 0.2] # 18% of the population are <20y. old, 62% are 20-66y. old, 20% are 67 or older

leth = [0.001, 0.02, 0.09]  # lethality for 3 age groups
rec = [1 - leth[0], 1 - leth[1], 1 - leth[2]]   # recovery rate for 2 age groups, calculated via lethality

 


# LOAD EXTERN DATA

adj_matrix  =   np.loadtxt('AdjacencyMuenster.csv', delimiter = ',', skiprows = 1, usecols = (1,2,3,4,5,6,7,8,9,10,11))
pop         =   np.loadtxt('Populations2.csv', delimiter = ',', usecols = (1)) 

M = adj_matrix + np.transpose(adj_matrix)  # add transposed matrix to make sure population is constant

Sn = np.zeros((len(ts),len(M)))
Sn[0,:] = pop                   # put in population numbers from ext. data
En = np.zeros((len(ts),len(M)))
En[0,0] = 20                    # let virus start in Muenster, 20 young people get it
Iyn = np.zeros((len(ts),len(M)))                    
Imn = np.zeros((len(ts),len(M)))
Ion = np.zeros((len(ts),len(M)))
Sn[0,0] = Sn[0,0] - En[0,0]     # in order to make that total start population in Münster is still according to 'Populations.csv'
Rn = np.zeros((len(ts),len(M)))
Dn = np.zeros((len(ts),len(M)))


U = [Sn[0,:], En[0,:], Iyn[0,:], Imn[0,:], Ion[0,:], Rn[0,:], Dn[0,:]]

# the following lines are only for the last plot in the right corner, where all the cities are added

sumS        =   np.zeros(np.shape(ts))   # array for total numbers of all cities
sumE        =   np.zeros(np.shape(ts))
sumIy       =   np.zeros(np.shape(ts))
sumIm       =   np.zeros(np.shape(ts))
sumIo       =   np.zeros(np.shape(ts))
sumR        =   np.zeros(np.shape(ts))
sumD        =   np.zeros(np.shape(ts))

sumS[0]     =   sum(Sn[0,:])            # generate first entry from ext. data
sumIy[0]    =   sum(Iyn[0,:])
sumIm[0]    =   sum(Imn[0,:])
sumIo[0]    =   sum(Ion[0,:])
sumR[0]     =   sum(Rn[0,:])
sumD[0]     =   sum(Dn[0,:])



for i in range(1,Nt):
    U = rk4_step(SIR, U, [alpha, beta, gamma, age_dist, leth, rec, M], dt)
    Sn[i], En[i], Iyn[i], Imn[i], Ion[i], Rn[i,:], Dn[i,:] = U                   # Sn[i] funktioniert hier wie Sn[i,:] 
    sumS[i] = sum(Sn[i,:])
    sumE[i] = sum(En[i,:])
    sumIy[i] = sum(Iyn[i,:])
    sumIm[i] = sum(Imn[i,:])
    sumIo[i] = sum(Ion[i,:])
    sumR[i] = sum(Rn[i,:])
    sumD[i] = sum(Dn[i,:])
#    
#    all_I[:,i] = Iyn[i,:] + Imn[i,:] + Ion[i:]
#    print(all_I[i])
#    print('ooo')
#    print(Ion[i])
#    sumall_I[i] = sum(all_I[i,:])
    
    
    
name_list = [[]] * len(M)               # to read first row of populations2.csv, to use for plot titles
c = 0
with open('Populations2.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    for row in readCSV:
        name_list[c] = row[0]
        c = c+1



# PLOT THE SUBPLOTS

fig1 = plt.figure(1)
plt.clf()



for c in range(0,len(M)):
    plt.subplot(3, 4, c+1)
    plt.plot(ts, Sn[:,c], 'C0-')    # plots Susceptible  
    plt.plot(ts, En[:,c], 'C9-', linewidth = 1)   # plots Exposed
    plt.plot(ts, Iyn[:,c], 'C3--', linewidth = 0.5)   # plots young Infected
    plt.plot(ts, Imn[:,c], 'C3--', linewidth = 0.6)   # plots middle-aged Infected
    plt.plot(ts, Ion[:,c], 'C3--', linewidth = 0.7)   # plots old Infected
    plt.plot(ts, Iyn[:,c]+Imn[:,c]+Ion[:,c], 'C3-') # plots ALL Infected
    plt.plot(ts, Rn[:,c], 'C2-')    # plots Recovered
    plt.plot(ts, Dn[:,c], 'C4-')    # plots dead people
    plt.title(name_list[c], fontsize=9)
    
for c in range(0,8):        # remove x label for first two rows
    plt.subplot(3, 4, c+1)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=True,      # ticks along the bottom edge are on
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
    
    
# plot sum of all regions
plt.subplot(3,4,12)
plt.plot(ts, sumS,'C0-')
plt.plot(ts, sumE, 'C9-', linewidth = 1)
plt.plot(ts, sumIy, 'C3--', linewidth = 0.5)
plt.plot(ts, sumIm, 'C3--', linewidth = 0.6)
plt.plot(ts, sumIo, 'C3--', linewidth = 0.7)
plt.plot(ts, sumIy+sumIm+sumIo, 'C3-')  # plots total sum of all Infected
plt.plot(ts, sumR, 'C2-')
plt.plot(ts, sumD, 'C4-')

plt.annotate(str(round(max(sumD)))+' dead!',(Tend,max(sumD)), xytext=(Tend+10, 0), arrowprops=dict(arrowstyle= '-|>')) # highlight total death count

plt.title('total numbers', fontsize=10)
plt.xlabel('t in days')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

plt.suptitle('SIR network with infection rate ' r'$\alpha='+ 
             str(alpha)+'$' '  & removal rate ' r'$\beta =' + str(beta)+'$' ' & transmission rate ' r'$\gamma='+ str(gamma)+'$', fontsize=19, fontweight='bold')

plt.show()

print('Overall ' + str(round(max(sumD))) + ' people have died after ' + str(Tend) + ' days.' )




