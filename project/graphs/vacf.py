import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

for i in range(1,6):
    df = pd.read_csv("tmp"+str(i)+".vacf")
    t = df.TimeStep
    t = (t - 500000)/10000
    r = df.c_vacf

    
    plt.title(r'VACF for '+str(i)+'th region', fontsize=14)
    plt.xlabel(r'Time, seconds', fontsize=14)
    plt.ylabel(r'VACF, (cm/s)$^2$', fontsize=14)
    plt.grid(True)
    plt.plot(t, r, color="firebrick")
    plt.savefig("vacf"+str(i)+".png")
    plt.clf()