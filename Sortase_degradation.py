#now we model enzyme degrading -> reaction constant decreases over time
#Try different values for c's

# import modules and libraries
import numpy as np
from scipy.special import factorial
import matplotlib.pyplot as plt
%matplotlib inline


ds = 1 #sortase is working at 100% efficiency

# reaction constants
c1 = 1*ds
c2 = 0.3*ds #this is reverse reaction 

# initial molecular population numbers of S1 and S2
A = 10
B = 30
C = 100
num_species = 3
# specify the maximum time, tmax, and maximum number of reactions, nrmax
# choose nrmax sufficiently high (for higher tmax, choose higher tmax)

tmax = 1
nrmax = 200 #number of reactions
# Initialise variables to store time and molecule numbers
current_t = 0
t_count = 0
largenum = 2*nrmax #twice the resolution of number of reactions
store_t = np.zeros(largenum) #initialize an array where we have filled it with zeros.  The number of zeros is double the number of steps we have
store_mols = np.zeros((largenum, num_species)) #make an array of arrays.  The num_species=2 means each of the nested arrays is [0, 0]. The largenum means we have length of array as 200
store_r = np.zeros(largenum)

mu = 0
next_r =0

while current_t < tmax:
    
    # store current system
    store_t[t_count] = current_t #this is keeping track of the frame that we are on for each segment of time
    store_mols[t_count,0] = A #for the matrix of matrix of zeros, update the value for the molecule A
    store_mols[t_count,1] = B
    store_mols[t_count,2] = C
    store_r[t_count] = next_r
    
    # calculate current reaction propensities
    a1 = A*B*c1
    a2 = C*c2  #whichever of a1 or a2 is greater will determine which of the two states (whether X1 decays first or X2 decays first so X1 levels rise then fall)
    
    ds = 1*((tmax - current_t + 1)/(tmax*(current_t+1))) #try different models of sortase enzyme degradation
    
    a = np.array([a1, a2])
    a0 = sum(a)
    
    #calculate the new reaction constants:
    c1 = 1*ds
    c2 = 0.3*ds #this is reverse reaction 
    
    #to debug and see why a0 is reaching zero, do print
    print(A)
    print(B)
    print(C)
    print(a)
    
     # generate two random numbers
    r = np.random.random(2)
    r1 = r[0]
    r2 = r[1]
    

    # choose next time increment 
    dt = (1/a0)*np.log(1/r1) #time interval afyter which next chemical reaction happens.  remember smaller the propensity of a reaction, a, the longer you have to wait
    #print("T: ", T)
    
    # choose next reaction 
    mu = 0
    N = r2*a0 - a[mu]

    while N > 0:
        print(N)
        mu = mu + 1
        N = N - a[mu]

    next_r = mu
    
    # define the time of next reaction
    current_t += dt

    # add one to the number of reactions count
    t_count += 1

    # update the system according to the reaction that was chosen
    if next_r == 0:
        A = A - 1
        B = B - 1
        C = C + 1
    
    if next_r == 1:
        C = C - 1
        A = A + 1
        B = B + 1


# plot of simulation
fig, ax = plt.subplots()
ax.plot(store_t, store_mols[:,0], '*', label='A')
ax.plot(store_t, store_mols[:,1], '*', label='B')
ax.plot(store_t, store_mols[:,2], '*', label='C')
legend = ax.legend(loc='right', shadow=True)
fig.suptitle('Equal Reaction Rates, A=200, B=100, C=0', fontsize=12)
ax.set_xlabel('Time (s)', fontsize=10)
ax.set_ylabel('Molecule Count', fontsize='medium')
plt.show()

# histogram of reactions chosen by algorithm
bins = np.array([-0.5,0.5,1.5])
plt.hist(store_r, bins)
