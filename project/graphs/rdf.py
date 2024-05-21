import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

for i in range(1,6):
    df = pd.read_csv("tmp"+str(i)+".rdf")
    r = df.c_2 / 1000
    rad = df.c_1 * 10
    # r = r*2*np.pi*rad
    
    plt.title(r'RDF for '+str(i)+'th region', fontsize=14)
    plt.xlabel(r'Radius, mm', fontsize=14)
    plt.ylabel(r'g(r)', fontsize=14)
    plt.grid(True)
    plt.plot(rad, r, color="firebrick")
    plt.savefig("rdf"+str(i)+".png")
    plt.clf()
