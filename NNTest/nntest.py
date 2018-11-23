import os, glob, pandas, numpy
import tensorflow as tf
from tensorflow import keras
from sklearn.model_selection import train_test_split

## simulate multivariate data sets

def simulate_synthetic_data(repeats=50, features=10):
    size = (repeats, features)  ## 50 features and 10 instances of input
    uniform = numpy.random.uniform(0, 1, size=size)
    normal = numpy.random.normal(0, 1, size=size)
    exp = numpy.random.exponential(size=size)

    uniform = pandas.DataFrame(uniform)
    normal = pandas.DataFrame(normal)
    exp = pandas.DataFrame(exp)

    uniform['label'] = ['uniform'] * uniform.shape[0]
    normal['label'] = ['normal'] * normal.shape[0]
    exp['label'] = ['exp'] * exp.shape[0]

    uniform.set_index('label', append=True, inplace=True)
    normal.set_index('label', append=True, inplace=True)
    exp.set_index('label', append=True, inplace=True)

    uniform.columns = ["Feature{}".format(i) for i in uniform.columns]
    normal.columns = ["Feature{}".format(i) for i in normal.columns]
    exp.columns = ["Feature{}".format(i) for i in exp.columns]

    df = pandas.concat([uniform, normal, exp])
    df = df.reset_index(drop=True, level=0)
    df['ID'] = range(df.shape[0])
    df = df.set_index('ID', append=True)
    df.index = df.index.swaplevel(0, 1)
    return df


df = simulate_synthetic_data(repeats=50, features=10)

train_test_split()












