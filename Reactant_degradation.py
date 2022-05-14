# import modules and libraries
import numpy as np
from scipy.special import factorial
import matplotlib.pyplot as plt
%matplotlib inline

#try 3 reactions with layout above to model the degradation of one of the reactants
#Try different values for c's
# reaction constants
c1 = 1
c2 = 1 #this is reverse reaction 
c3 = 10 #when this is high, A decays  super fast
# initial molecular population numbers of S1 and S2
A = 100
B = 100
C = 0
num_species = 3
# specify the maximum time, tmax, and maximum number of reactions, nrmax
# choose nrmax sufficiently high (for higher tmax, choose higher tmax)

tmax = 500
nrmax = 500 #number of reactions
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
    a3 = A*c3
    a = np.array([a1, a2, a3])
    a0 = sum(a)
    
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
    dt = (1/a0)*np.log(1/r1) #time interval after which next chemical reaction happens.  remember smaller the propensity of a reaction, a, the longer you have to wait
    #print("T: ", T)
    
    # choose next reaction 
    mu = 0
    N = r2*a0 - a[mu]

    while N > 0:
        print(N)
        mu = mu + 1
        #if mu >= len(a):
        #    mu = int(len(a)) - 1
        #xcom = int(a[mu])
        N = N - a[mu] #xcom

    next_r = mu
    
    # define the time of next reaction
    current_t += dt

    # add one to the number of reactions count
    t_count += 1

    # update the system according to the reaction that was chosen
    if next_r == 0:
        print("reaction 1")
        A = A - 1
        B = B - 1
        C = C + 1
    
    if next_r == 1:
        print("reaction 2")
        C = C - 1
        A = A + 1
        B = B + 1
        
    if next_r == 2:
        print("reaction 3")
        A = A - 1
        
    if next_r > 2:
        print("error mu too high")
        

# plot of simulation
fig, ax = plt.subplots()
ax.plot(store_t, store_mols[:,0], 'g.', markersize=12, label='A')
ax.plot(store_t, store_mols[:,1], 'r.', markersize=3, label='B')
ax.plot(store_t, store_mols[:,2], 'k.', markersize=3, label='C')
legend = ax.legend(loc='upper right', shadow=True)
fig.suptitle('Reactant A Decay, A=100, B=100, C=0', fontsize=12)
ax.set_xlabel('Time', fontsize=10)
ax.set_ylabel('Molecule Count', fontsize='medium')
plt.show()



# histogram of reactions chosen by algorithm
print(store_r)
bins = np.array([-0.5,0.5,1.5, 2.5])
plt.hist(store_r, bins)
