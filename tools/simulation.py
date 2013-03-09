from rpy import *
import numpy as np
def rand_model(number_genes=10000,record_length=12,type='white'):
    """generate the white noise or AR(1) simulation datasets"""
    #time_points = '\t'.join(map(str,range(0,45,4)))
    datasets = {}

    # generate random time series of white noise and AR(1)
    r.set_seed(314159265)
    for i in range(number_genes):
        if type == 'white':
            data = r.rnorm(record_length)
        else:
            set_default_mode(NO_CONVERSION)
            ar1 = r.arima_sim(r.list(ar=0.5),n=record_length)
            set_default_mode(BASIC_CONVERSION)
            data = ar1.as_py()
        datasets[type+'_'+str(i)] = data
    return datasets

def cosine_model_fixed(number_genes=5000,record_length=12,type='simple'):
    """generate simple or complicated cosine curve model"""
    datasets = {}
    amp = 3.5   # signal/noise ratio = 3.5:1
    mea = 0     # mean level of cosine model
    pha = 0     # phase of cosine model
    per = 24    # period of the cosine model
    r.set_seed(1)
    x = np.arange(0,45,4)
    for i in range(number_genes):
        if type == 'simple':
            data = mea + amp*cos(2*np.pi/per * x + pha) + r.rnorm(record_length)
        else:
            mea = 500
            amp = 140
            sn = 3.5    # signal-noise ratio 
            k = 0.01    # the half life of amplitude and mean, it takes for the maximum value to decrease by a factor of 2. It is computed as 0.693/k
            data = mea*np.exp(-k*x) + amp*np.exp(-k*x)*cos(2*np.pi/per*x + pha) + r.rnorm(record_length,sd=amp/sn)
        datasets[type+'_'+str(i)] = data
    return datasets

def cosine_model_varied(record_length=12,type='sta'):
    """generate stationary or nonstationary cosine curve model"""
    datasets = {}
    r.set_seed(314159265)
    x = np.arange(0,45,4)
    for per in arange(20,28,0.1):
        for pha in range(25):
            for sn in range(1,6,1):
                if type == 'sta':
                    amp = 2     # signal/noise ratio = 1:1
                    mea = 0     # mean level of cosine model
                    data = mea +sn*amp*cos(2*np.pi*(x/per - 1/25*pha)) + r.rnorm(record_length)
                else:
                    mea = 500
                    amp = 100
                    k = 0.01    # the half life of amplitude and mean, it takes for the maximum value to decrease by a factor of 2. It is computed as 0.693/k
                    data = mea*np.exp(-k*x) + sn*amp*np.exp(-k*x)*cos(2*np.pi*(x/per - 1/25*pha)) + r.rnorm(record_length,sd=50)
                datasets['%s_%s_%s_%s'%(type,per,pha,sn)] = data
    return datasets

def cosine_model_multiple(number_genes=2500,record_length=12,type='sta'):
    """generate composed two cosine stationary models"""
    datasets = {}
    r.set_seed(314159265)
    x = np.arange(0,45,4)
    for i in range(number_genes):
        for step in range(1,5,1):
            amp = 3     # signal/noise ratio = 3:1
            mea = 0     # mean level of cosine model
            per1 = 24 + step
            per2 = 24 - step
            data = mea + amp*cos(2*np.pi*(x/per1)) + amp*cos(2*np.pi*(x/per2)) + r.rnorm(record_length)
            datasets['%s%d_%s_%s'%(type,i,per1,per2)] = data
    return datasets

if __name__ == '__main__':
    fou = open('multiple_cosine.txt','w')
    data = cosine_model_multiple()
    for name in data:
        value = '\t'.join(map(str,data[name]))
        print >>fou, name+'\t'+value
