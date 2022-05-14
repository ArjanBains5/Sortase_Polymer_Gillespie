# import modules and libraries
import numpy as np
from scipy.special import factorial
import matplotlib.pyplot as plt
%matplotlib inline


# FOR PROJECT
#try big reaction system  with layout above WITH ENZYME DYING OVER TIME

# RUN SIMULATION
# R1: A + B -> C
# R2: C -> A + B
# R3: A -> 0 a.k.a. A degrades
# R4: C + B -> D
# R5: D -> C + B
# R6: D + B -> E
# R7: E -> D + B
# R8: E + B -> F
# R9: F -> E + B
# R10: F + B -> G
# R11: G -> F + B

# reaction constants; include ds term for death of sortase over time
ds = 1 #sortase is working at 100% efficiency

c1 = 1*ds
c2 = 0.3*ds #this is reverse reaction 
c3 = 0.2 #when this is high, A decays  super fast
c4 = 1*ds
c5 = 0.3*ds
c6 = 1*ds
c7 = 0.3*ds
c8 = 1*ds
c9 = 0.3*ds
c10 = 1*ds
c11 = 0.3*ds

# initial molecular population numbers of S1 and S2
A = 300
B = 300
C = 150
D = 0
E = 0
F = 0
G = 0

num_species = 7
# specify the maximum time, tmax, and maximum number of reactions, nrmax
# choose nrmax sufficiently high (for higher tmax, choose higher tmax)

tmax = 500
nrmax = 20000 #number of reactions
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
    store_mols[t_count,3] = D
    store_mols[t_count,4] = E
    store_mols[t_count,5] = F
    store_mols[t_count,6] = G
    store_r[t_count] = next_r
    
    # calculate current reaction propensities
    a1 = A*B*c1
    a2 = C*c2  #whichever of a1 or a2 is greater will determine which of the two states (whether X1 decays first or X2 decays first so X1 levels rise then fall)
    a3 = A*c3
    a4 = C*B*c4
    a5 = D*c5
    a6 = D*B*c6
    a7 = E*c7
    a8 = E*B*c8
    a9 = F*c9
    a10 = F*B*c10
    a11 = G*c11
    
    ds = 1*((tmax - current_t + 1)/(tmax*(current_t+1))) #try different models of sortase enzyme degradation
    
    a = np.array([a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11])
    a0 = sum(a)
    
    #calculate the new reaction constants:
    c1 = 1*ds
    c2 = 0.3*ds #this is reverse reaction 
    c3 = 0.2 #when this is high, A decays  super fast
    c4 = 1*ds
    c5 = 0.3*ds
    c6 = 1*ds
    c7 = 0.3*ds
    c8 = 1*ds
    c9 = 0.3*ds
    c10 = 1*ds
    c11 = 0.3*ds
    
    #to debug and see why a0 is reaching zero, do print
    print(A)
    print(B)
    print(C)
    print(D)
    print(E)
    print(F)
    print(G)
    print(a)
    print(ds)
    print("c1 is:", c1)
    
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
        
    if next_r == 3:
        print("reaction 4")
        C = C - 1
        B = B - 1
        D = D + 1
        # R4: C + B -> D
        
    if next_r == 4:
        print("reaction 5")
        C = C + 1
        B = B + 1
        D = D - 1
        # R5: D -> C + B
        
    if next_r == 5:
        print("reaction 6")
        E = E + 1
        B = B - 1
        D = D - 1
# R6: D + B -> E

    if next_r == 6:
        print("reaction 7")
        E = E - 1
        B = B + 1
        D = D + 1
# R7: E -> D + B

    if next_r == 7:
        print("reaction 8")
        E = E - 1
        B = B - 1
        F = F + 1
# R8: E + B -> F

    if next_r == 8:
        print("reaction 9")
        F = F - 1
        B = B + 1
        E = E + 1
# R9: F -> E + B

    if next_r == 9:
        print("reaction 10")
        F = F - 1
        B = B - 1
        G = G + 1
# R10: F + B -> G

    if next_r == 10:
        print("reaction 11")
        G = G - 1
        B = B + 1
        F = F + 1
# R11: G -> F + B
        
    if next_r > 11:
        print("error mu too high")
   
  
# plot of simulation
fig, ax = plt.subplots()
ax.plot(store_t, store_mols[:,0], '.', label='A')
ax.plot(store_t, store_mols[:,1], '.', label='B')
ax.plot(store_t, store_mols[:,2], '.', label='C')
ax.plot(store_t, store_mols[:,3], '.', label='D')
ax.plot(store_t, store_mols[:,4], '.', label='E')
ax.plot(store_t, store_mols[:,5], '.', label='F')
ax.plot(store_t, store_mols[:,6], '.', label='G')
fig.suptitle('Sortase Decay Polymer Creation, A=300, B=300, C=150', fontsize=12)
ax.set_xlabel('Time', fontsize=10)
ax.set_ylabel('Molecule Count', fontsize='medium')
legend = ax.legend(loc='right', shadow=True)
plt.show()


# histogram of reactions chosen by algorithm
import matplotlib.pyplot as plt

print(store_r)
bins = np.array([-0.5,0.5,1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5])
#plt.xlabel
plt.ylabel("Number of Times Reaction Occurred")
plt.hist(store_r, bins, color='r')
