import tellurium
import matplotlib
matplotlib.use('TkAgg')
import os, glob, numpy, pandas
import matplotlib.pyplot as plt
import seaborn

plt.figure()

for i in range(25):
    A0, B0 = numpy.random.uniform(0, 100, 2)
    ant = f"""
    model new_model()
        R1: A => B; k1*A;
        R2: B => A; k2*B;
    
        A = 50;//{A0};
        B = {B0};
    
        k1 = 3;
        k2 = 4;
    end
    """
    mod = tellurium.loada(ant)
    data = mod.simulate(0, 5, 11)
    plt.plot(data[:, 0], data[:, 1], 'r', label='A')
    plt.plot(data[:, 0], data[:, 2], 'g', label='B')

plt.show()












