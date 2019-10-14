import numpy as np
import matplotlib.pyplot as plt

# spectral density files

'''
file1= open("drude_spectral_density_2.txt")
file2= open("superohmic_spectral_density_3.txt")
file3 = open("fractionalfunction_spectral_density_1.txt" )
file4 =open("lognormal_spectral_density_4.txt")

labels = ['Drude', 'Super-Ohmic', 'Fractional', 'Lognormal']

labels = ['Drudė', 'Super-ohinis', 'Trupmeninis', 'Lognormalusis']
'''
file1 = open("fractionalfunction_spectral_density_1.txt" )

labels =['Trupmeninis']
files = [file1]

def frequencylist(data):
    freqlist =[]
    for i in range(len(files)):
            pn = int(data[i][0])
            initfreq = data[i][1]
            step = data[i][2]
            freq= np.zeros(pn)
            for j in range(pn):
                freq[j] =initfreq+step*j            
            freqlist.append(freq)
    return freqlist

def read(file):
    data = []
    for i in range(len(file)):
        X = []
        for line in file[i].readlines():
            if line.startswith("#"):
                continue
            y = [float(value) for value in line.split()]
            X.append(y)
        X =np.array(X)
        data.append(X)
        file[i].close()
    return(data)       

def pictureprint(freqstorage, storage):
    fig, ax = plt.subplots()
    for i in range(len(files)):
        plt.ylim(0,45)
        plt.xlim(0,800)  
        ax.plot(freqstorage[i],storage[i][3:],lw=2, label=labels[i])
    #ax.legend(title = r'$\alpha$ =', loc="upper right", shadow=False)
    ax.legend(loc="upper right", shadow=False, fontsize=12)
    plt.ylabel('C"($\omega), cm^{-1}$', fontsize=12)
    plt.xlabel('$\omega , cm^{-1}$', fontsize=12)
    for i in range(len(files)):
        additionalfig = fig.add_axes([0.3, 0.45, 0.3, 0.35])
        additionalfig.plot(freqstorage[i],storage[i][3:],lw=2, label=labels[i])
        additionalfig.set_xlim(0,10)
        additionalfig.set_ylim(0,10)
    return()

storage = read(files)
freqstorage = frequencylist(storage)
pictureprint(freqstorage, storage)
