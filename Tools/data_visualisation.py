import matplotlib.pyplot as plt
import numpy as np




def Bin_Analysis( array):
    '''
    The analasys takes the last 40% of the data where from visual observation
    we see that it has equibrilated. And split in 40 blocks of which 30 are used
    for binning statistics and then are left unused to decorelate the bins.
    @In : The complete array
    @out: Average, Stdv (Standart deviation)
    '''

    if type(array) != np.ndarray:   # If it's not a numpy array
        array = np.array(array)     # Make it a Numpy array ;)

    averages = []
    N = len(array)  
    for n in range(6*N/10, (N+1)-N/25, N/25): #Remember that in Python 2.7 int/int=int 
        average_n = np.average(array[n:(n+3*N/100)])
        averages.append(average_n)

    Average = np.average(averages)
    Stdv = 2*np.std(averages)/np.sqrt(len(averages)) # 95% certainty interval

    #print Average, Stdv
    return int(np.round(Average)), int(np.round(Stdv))

                            
    

# Read in Energies and Momenta

Epot = []   # Polarization
TEner = []   # Total Energy
time = []   # Time

with open("ener.dat", 'r') as data:
    for line in data:
        words = line.split()
        try:
            words = [float(words[i]) for i in range(len(words))]
            time.append(words[0])
            TEner.append(words[4])
            Epot.append(words[5])
            
        except:
            pass

# Plot Polarization
x_list = time
plt.plot(x_list, Epot, label="Polarization: %i (%i))"%(Bin_Analysis(Epot)))
plt.legend()

plt.title('Visualisation of Polarization, values given are average of last 40% with 95% certainty interval')
plt.xlabel('time (reduced)')
plt.ylabel('Polarization (reduced)')

plt.savefig('pol.png', bbox_inches='tight')
#plt.show()
plt.clf()
#'''
# Plot Energy
x_list = time
plt.plot(x_list, TEner, label="Total Energy: %i (%i))"%(Bin_Analysis(Epot)))
plt.legend()

plt.title('Visualisation of Total Energy, values given are average of last 40% with 95% certainty interval')
plt.xlabel('Time (reduced)')
plt.ylabel('Energy (reduced)')

plt.savefig('energy.png', bbox_inches='tight')
#plt.show()
plt.clf()

#'''
