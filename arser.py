'''
ARSER: Analysis circadian expression data by harmonic regression 
based on autoregressive spectral estimation
Usage: pyhton arser.py inputfile outputfile
Author: Rendong Yang
Email: cauyrd@gmail.com
Version: 2.0
'''
import lib.utility as HR
import rpy2.robjects as robjects
from pylab import detrend_linear
import numpy as np
import sys
r = robjects.r
spec_ar = r['spec.ar']
p_adjust = r['p.adjust']

class Arser:
    '''
    Doc: Class for harmonic regression with MESE period identification
    Input: 
          x values -> time points
          y values -> expression value
    Output:
          self.x -> raw x values
          self.y -> raw y values
          self.mean -> mean value for raw y values
          self.dt_y -> y values detrended
          self.delta -> sampling interval
          self.estimate_period -> period identified by MESE
          self.amplitude -> amplitude for single cosine model
          self.phase -> phase for single cosine model
          self.R2 -> R square of regression curve
          self.R2adj -> adjusted R square of regression curve
          self.coef_var -> (standard deviation)/mean
          self.pvalue -> F test for testing significant regression model
    '''

    def __init__(self, x, y):
        '''
            initialized Arser instance object, detrendedy values, mean, period are calculated
        '''
        self.x = np.array(x)
        self.y = np.array(y)
        self.mean = np.mean(self.y)
        self.dt_y = detrend_linear(self.y)
        self.delta = self.x[1] - self.x[0]
        self.R2 = -1
        self.R2adj = -1
        self.coef_var = -1
        self.pvalue = -1
        self.phase = []
        self.estimate_period = []
        self.amplitude = []
        

    def get_period(self, is_filter=True, ar_method='mle'):
        '''
        estimate possiable cycling period of time series
        '''
        num_freq_mese = 500
        set_order = 24/self.delta
        if (set_order == len(self.x)):
            set_order = len(self.x)/2
        try:
            filter_y = HR.savitzky_golay(self.dt_y)
        except:
            filter_y = HR.savitzky_golay(self.dt_y, kernel=5, order=2)
        if is_filter:
            try:
                mese = spec_ar(robjects.FloatVector(filter_y.tolist()), n_freq=num_freq_mese, plot=False, method=ar_method, order=set_order)
            except:
                print 'spec_ar running error at line 70'
                sys.exit(1)
        else:
            try:
                mese = spec_ar(robjects.FloatVector(self.dt_y.tolist()), n_freq=num_freq_mese, plot=False, method=ar_method, order=set_order)
            except:
                print 'spec_ar running error at line 76'
                sys.exit(1)
 
        # search all the locial peaks of maximum entropy spectrum
        peaks_loc = []      # the locition for each peak in mese spectrum
        for i in range(1, num_freq_mese-1):
            if mese.rx2('spec')[i] > mese.rx2('spec')[i+1] and mese.rx2('spec')[i] > mese.rx2('spec')[i-1]:
                peaks_loc.append((mese.rx2('spec')[i], i))
        peaks_loc.sort(reverse=True)    # sort frequency by spectrum value
        try:
            periods = [1/mese.rx2('freq')[item[1]]*self.delta for item in peaks_loc]
        except:
            periods = []
        return periods

    def harmonic_regression(self, period):    
        '''
        general harmonic regression
        dt_y = mesor + sigma( A*cos(2*pi/T*x) + B*sin(2*pi/T*x) ) + error
        '''
        x = np.array([])
        x = x.reshape(len(self.x), 0)
        x_varnm_names = []
        for T in period:
            cosx = np.cos(2*np.pi/T*self.x)
            sinx = np.sin(2*np.pi/T*self.x)
            x = np.c_[x, cosx, sinx]
            x_varnm_names += ['cos%.1f' % T,'sin%.1f' % T]
        model = HR.ols(self.dt_y, x, y_varnm = 'y', x_varnm = x_varnm_names)
        return model

    def evaluate(self, T_start=20, T_end=28, T_default=24):    
        '''
        evaluate the best model for each time series
        '''
        is_filter = [True, False]
        ar_method = ['yule-walker', 'mle', 'burg']
        best_model = {'AIC':1e6, 'ols_obj':None, 'period':None, 'filter':None, 'ar_method':''}
        for p1 in is_filter:
            for p2 in ar_method:
                
                # choose best model's period from 'mle','yw','burg'
                period = filter((lambda x:x>=T_start and x<=T_end), self.get_period(is_filter=p1, ar_method=p2))
                if not period:
                    p2 = 'default'
                    period = [T_default]
                m = self.harmonic_regression(period)

                # model selection by aic
                aic = m.ll()[1]
                if aic <= best_model['AIC']:
                    best_model['AIC'] = aic
                    best_model['ols_obj'] = m
                    best_model['period'] = period
                    best_model['filter'] = p1
                    best_model['ar_method'] = p2

        # record the best model parameters
        self.estimate_period = best_model['period']
        self.amplitude = []
        self.phase = []
        m = best_model['ols_obj']
        for i in range(len(self.estimate_period)):
            phi = np.angle(np.complex(m.b[2*i+1], -m.b[2*i+2]))
            
            # for float point number can not compare with 0, so use <1e-6 as nonpositive real number
            self.phase.append(np.fabs(phi)/(2*np.pi)*self.estimate_period[i] if phi<=1e-6 else self.estimate_period[i] - phi/(2*np.pi)*self.estimate_period[i]) 
            self.amplitude.append(np.sqrt(m.b[2*i+1]**2 + m.b[2*i+2]**2))
        self.R2 = m.R2
        self.R2adj = m.R2adj
        self.pvalue = m.Fpv
        self.coef_var = np.std(self.y)/self.mean
        return {'period':self.estimate_period, 'amplitude':self.amplitude, 'phase':self.phase, 'R2':self.R2, 'R2adj':self.R2adj, 'coef_var':self.coef_var, 'pvalue':self.pvalue, 'filter':best_model['filter'], 'ar_method':best_model['ar_method']}


if __name__ == '__main__':
    import time
    import os
    import random
    run_start = time.time()
    random.seed(run_start)
    if len(sys.argv)==1:
		print 'version 2.0'
        print 'Usage: python arser.py inputfile outputfile start(optional) end(optional) default_period(optional)'
        print 'start: default 20'
        print 'end: default 28'
        print 'default_period: 24'
        sys.exit(0)
    fin = open(sys.argv[1])
    tempfile = "temp" + str(random.random())
    fou = open(tempfile,'w')
    start = 20
    end = 28
    T_default = 24
    if len(sys.argv)>3:
        start = float(sys.argv[3])
        end = float(sys.argv[4])
        T_default = float(sys.argv[5])

    # get time_points at header of input file
    time_points = map(float, fin.readline().split()[1:])
    print >> fou, 'probe\tfilter_type\tar_method\tperiod_number\tperiod\tamplitude\tphase\tmean\tR_square\tR2_adjust\tcoef_var\tpvalue'
    
    # estimate each probe using general  harmonic regression
    pvalues = []
    for line in fin:
        probe  = line.split()
        y_value = map(float, probe[1:])
        arser = Arser(time_points, y_value)
        d = arser.evaluate(start, end, T_default)    # search period in frequency domain range:[start,end]
        opt_line = [probe[0]]	# probe name column
        opt_line.append(int(d['filter']))	# model type column
        opt_line.append(d['ar_method'])    # arse method column
        opt_line.append(len(d['period']))	# period_number column
        opt_line.append(','.join(map(str, d['period'])))    # period colume	
        opt_line.append(','.join(map(str, d['amplitude'])))	# amplitude column
        opt_line.append(','.join(map(str, d['phase'])))    # phase colume
        opt_line.append(arser.mean) # mean column
        opt_line.append(d['R2'])  # R_square column
        opt_line.append(d['R2adj']) # R2_adj column
        opt_line.append(d['coef_var'])	# coef_var column
        opt_line.append(d['pvalue'])	# pvalue column
        print >> fou, '\t'.join(map(str, opt_line))
        pvalues.append(d['pvalue'])
    fin.close()
    fou.close()
    
    # calculate q-value for all the probes
    fin = open(tempfile)
    fou = open(sys.argv[2],'w')
    head = fin.readline().split()
    #head.insert(12,'qvalue')
    head.insert(12,'fdr_BH')
    print >>fou, '\t'.join(head)
    #r.library('siggenes')
    #pi0 = r.pi0_est(pvalues)
    #if pi0['p0'] == 0:
    #    pi0['p0'] = 0.95
    #qvalues = r.qvalue_cal(pvalues,pi0['p0'])
    qvalues_BH = p_adjust(pvalues,'BH')
    for i,line in enumerate(fin):
        data = line.split()
        #data.insert(12,str(qvalues[i]))
        data.insert(12,str(qvalues_BH[i]))
        print >>fou, '\t'.join(data)
    os.remove(tempfile)
    print 'time used:', str(time.time()-run_start), 'seconds'
