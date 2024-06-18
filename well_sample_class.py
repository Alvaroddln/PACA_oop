import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats, integrate
from scipy.optimize import curve_fit
from copy import deepcopy
from functions import*

colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange', 'tab:cyan', 'tab:purple', 'tab:pink', 'tab:brown',
          'tab:olive', 'lightcoral', 'seagreen', 'mediumpurple', 'gold', 'sandybrown', 'deepskyblue', 'palegreen',
          'violet', 'darkturquoise', 'limegreen', 'lightsalmon']

class well:
    def __init__(self, well_name, growth_dataframe, sample_info_dict, blk_info_dict=False, trim=False):
        #name and wells of sample replicates
        self.name = well_name
        sample_info_dict_copy = deepcopy(sample_info_dict)
        for i, s in enumerate(sample_info_dict_copy):
            if self.name in sample_info_dict_copy[s]:
                self.sample = s
                self.color = colors[i]
        #compute time and bouts
        self.time = growth_dataframe.index.to_list()
        self.bouts = self.time[1]-self.time[0]
        #trimming option
        if trim != False:
            self.trim_time = trim
            #update time
            for i, value in enumerate(self.time):
                if value == self.trim_time:
                    self.time = self.time[:i+1]
            #update dataframe
            growth_dataframe = growth_dataframe[:self.time[-1]]
        else:
            pass   
        #raw data with blanks subtracted
        if blk_info_dict:
            blk_mean = growth_dataframe[blk_info_dict['blank']].mean(axis=1)
            self.data_raw = growth_dataframe[self.name].sub(blk_mean, axis=0).clip(lower=blk_mean)
        else:
            self.data_raw = growth_dataframe[self.name]
        self.raw_max_od = max(self.data_raw)
        self.raw_0 = list(self.data_raw)[0]
        #normalized data
        min_data = min(self.data_raw)
        self.data = self.data_raw/min_data #normalization by dividing by min value
        #basic parameters
        self.max_od = max(self.data)
        self.max_od_time_idx = self.data[self.data == self.max_od].index[0]
        self._0 = list(self.data)[0]
        #ln data
        self.data_ln = np.log(self.data)
        #rate
        self.rate_ln = pd.Series(np.gradient(self.data_ln, self.bouts), index=self.time)
        self.max_rate_ln = max(self.rate_ln)
        
    def compute_algorithm_points_richards(self):
        if self.R == True:
            x = np.linspace(self.time[0], self.time[-1], 1000)
            y = Richards(x, self.richards_A, self.richards_K, self.richards_r, self.richards_gamma, self.richards_nu)
            z = Richards_ln_diff(x, self.richards_A, self.richards_K, self.richards_r, self.richards_gamma, self.richards_nu)
    
            #2 * dt
            self.two_dt = 2 * self.dt
            #window calculation
            model_bout = self.time[-1]/1000        
            self.window_length = int(round_to_bout(self.two_dt, model_bout)/model_bout)
            
            windows_with_max = []
            for i in range(len(list(z))-self.window_length-1):
                window_temp = list(z)[i:i+self.window_length]
                if self.mu_max in window_temp:
                    windows_with_max.append((np.var(window_temp),
                                             list(range(i,i+self.window_length))))
            if len(windows_with_max) == 0:
                raise Exception('Not enough points for calculations')
            self.model_x_1, self.model_x_2 = min(windows_with_max)[1][0], min(windows_with_max)[1][-1] + 1
            self.model_alg_idx = slice(self.model_x_1, self.model_x_2) #output a slice object for list comprehension indices

            self.model_exponential = x[self.model_alg_idx]

            self.exponential = []
            for t in self.time:
                if x[self.model_x_1] < t < x[self.model_x_2]:
                    self.exponential.append(t)
            self.alg_idx = (self.exponential[0], self.exponential[-1])
    
    def compute_Richards(self):
        self.data_max = self.data[:self.max_od_time_idx]
        self.time_max = self.data_max.index
        self.R = False
        h = 1
        while self.R == False and h < 10:
            try:
                self.richards = curve_fit(Richards, self.time_max, self.data_max, full_output=True) #p0=[self._0, self.max_od, 1, 1, 1]
                self.richards_A, self.richards_K, self.richards_r, self.richards_gamma, self.richards_nu = self.richards[0]
                self.cost = np.sum(np.square(self.richards[2]['fvec']))
                if np.inf in self.richards[1]:
                    raise ValueError('Nope')
                #print('Fit found adding')
                self.R = True
                break
            except:
                self.data_max = self.data[:self.max_od_time_idx+self.bouts*h]
                self.time_max = self.data_max.index
                self.R = False
                #print('No fit was found for well %s adding' %self.name)
                if self.time_max[-1] == self.time[-1]:
                    h = 20
                h += 1

        if self.R == False:
            self.data_max = self.data[:self.max_od_time_idx]
            self.time_max = self.data_max.index
            k = 1
            while  self.R == False and k < 10:
    
                try:
                    self.richards = curve_fit(Richards, self.time_max, self.data_max, full_output=True) #p0=[self._0, self.max_od, 1, 1, 1]
                    self.richards_A, self.richards_K, self.richards_r, self.richards_gamma, self.richards_nu = self.richards[0]
                    self.cost = np.sum(np.square(self.richards[2]['fvec']))
                    if np.inf in self.richards[1]:
                        raise ValueError('Nope')
                    #print('Fit found subtracting')
                    self.R = True
                    break
                except:
                    self.data_max = self.data[:self.max_od_time_idx-self.bouts*k]
                    self.time_max = self.data_max.index
                    self.R = False
                    #print('No fit was found for well %s subtracting' %self.name)
                    k += 1
    
    def plot_richards(self):
        if self.R == True:
            x = np.linspace(self.time[0], self.time[-1], 1000)
            y = Richards(x, self.richards_A, self.richards_K, self.richards_r, self.richards_gamma, self.richards_nu)
            z = Richards_ln_diff(x, self.richards_A, self.richards_K, self.richards_r, self.richards_gamma, self.richards_nu)
            fig, (ax1, ax2) = plt.subplots(1,2,figsize=(12, 4))
            ax1.plot(x, y, '--', label='Fit', color='black', alpha=0.5)
            ax1.plot(self.time, self.data, 'c.', label='Data', c=self.color)
            ax1.legend()
            ax1.set_xlabel('Time (h)')
            ax1.set_ylabel('Normalized growth')
            
            #ax1.axvline(self.A_asymptote, color='cyan')
            #ax1.axvline(self.K_asymptote, color='green')
            ax2.plot(x, z, 'b--', label='Fit', color='black', alpha=0.5)
            ax2.plot(self.time, self.rate_ln, 'c.', label='Data', c=self.color)
            ax2.legend()
            ax2.set_xlabel('Time (h)')
            ax2.set_ylabel('Growth rate')

            #ax2.axvline(x[self.model_x_1], ls='--', color='black', alpha=0.5)
            #ax2.axvline(x[self.model_x_2], ls='--', color='black', alpha=0.5)
            
            fig.suptitle(self.name)
            plt.show()

    
    def kinetic_richards(self):
        if self.R == True:
            x = np.linspace(self.time[0], self.time[-1], 1000)
            y = Richards(x, self.richards_A, self.richards_K, self.richards_r, self.richards_gamma, self.richards_nu)
            z = Richards_ln_diff(x, self.richards_A, self.richards_K, self.richards_r, self.richards_gamma, self.richards_nu)

            self.fit_x, self.fit_y = x, y

            search_A, search_K = 1,1
            for i, j in zip(x,y):
                if j > self.richards_A + 0.1 and search_A ==1 :
                    self.A_asymptote = i
                    search_A = 0
                elif j > self.richards_K -0.1 and search_K == 1:
                    self.K_asymptote = i
                    search_K = 0
                else:
                    pass
            if search_A == 1:
                self.A_asymptote = 0
            if search_K == 1:
                self.K_asymptote = self.time[-1]
                
            self.model_width = self.K_asymptote - self.A_asymptote
            self.data_width = self.max_od_time_idx - self.A_asymptote
            
            # plt.axvline(self.A_asymptote, c='r')
            # plt.axvline(self.K_asymptote, c='b')
            # plt.axvline(self.max_od_time_idx, c='g')
            
            if self.data_width < self.model_width*(2/3):
                self.alert = 'O'
            elif self.data_width > self.model_width*(4/3):
                self.alert = 'U'
            else:
                self.alert = ''

            
            self.mu_max = max(z[np.isfinite(z)]) #find max ignoring nans and infs, sometimes it happens
            for i, j, k in zip(x,y,z):
                if k == self.mu_max:
                    self.mu_max_idx = i
                    self.ln_model_mu_max = np.log(j)
                    break
                else:
                    pass
            self.dt = np.log(2)/self.mu_max
            self.intercept = self.ln_model_mu_max - self.mu_max * self.mu_max_idx

            if self.richards_A > 0:
                self.lag = (np.log(self.richards_A) - self.intercept)/self.mu_max
            else:
                self.lag = (np.log(self._0) - self.intercept)/self.mu_max

            # plt.plot(x, np.log(y))
            # plt.axvline(self.lag)
            # xx = np.linspace(self.A_asymptote, self.K_asymptote)
            # plt.plot(xx, linear(xx, self.mu_max, self.intercept))
            # plt.plot(self.mu_max_idx, self.ln_model_mu_max, 'ro')
            # plt.axhline(np.log(self.richards_A))

            self.parameters = self.mu_max, self.max_od, self.lag, self.alert, self.cost

    def compute_activity(self, reporter_dataframe, use_model_data=True):
        self.reporter_data = reporter_dataframe[self.name]
        self.reporter_rate = pd.Series(np.gradient(self.reporter_data, self.bouts), index=self.time)
        if self.R == True:
            idcs = find_closest_data(self.fit_x, self.time)
            self.data_from_model = self.fit_y[idcs]
            self.dynact = self.reporter_rate/self.data_from_model
            
            if use_model_data == False:
                self.dynact = self.reporter_rate/self.data
                
            a, b = self.alg_idx
            self.dynact_exp = self.dynact[a:b]
            self.activity = np.mean(self.dynact_exp)
            
            #check linearity between growth and reporter
            self.RG_reg = stats.linregress(self.data[a:b], self.reporter_data[a:b])
            self.RG_rsq = self.RG_reg.rvalue ** 2       
            
            # plt.plot(self.time, self.dynact, 'b-')
            # plt.plot(self.exponential, self.dynact[a:b], 'ro')

class sample:
    def __init__(self, sample_name, growth_dataframe, sample_info_dict, blk_info_dict):
        self.name = sample_name
        self.time = growth_dataframe.index.to_list()
        self.wells = sample_info_dict[self.name]
        
        self.wells_dict = {}
        for w in self.wells:
            w_obj = well(w, growth_dataframe, sample_info_dict, blk_info_dict)
            self.wells_dict.update({w : w_obj})
            self.color = w_obj.color
            
        #raw well data
        self.data_raw = pd.DataFrame([self.wells_dict.get(w).data_raw for w in self.wells], index=self.wells).T
        self.data_raw_mean = self.data_raw.mean(axis=1)
        self.data_raw_std = self.data_raw.std(axis=1)
        
        self.raw_max_od = max(self.data_raw_mean)
        self.raw_0 = list(self.data_raw_mean)[0]
        
        #normalized well data
        self.data = pd.DataFrame([self.wells_dict.get(w).data for w in self.wells], index=self.wells).T
        self.data_mean = self.data.mean(axis=1)
        self.data_std = self.data.std(axis=1)        
        
        self.max_od = max(self.data_mean)
        self.max_od_time_idx = self.data_mean[self.data_mean == self.max_od].index[0]
        self._0 = list(self.data)[0]

        self.mean_fake_well = well(self.name, pd.DataFrame(self.data_mean, columns=[self.name]), {self.name:self.name})

    def compute_Richards(self):
        self.mean_fake_well.compute_Richards()
        if self.mean_fake_well.R == True:
            self.R = True
        else:
            self.R = False

    def kinetic_richards(self):
        if self.R == True:
            self.mean_fake_well.kinetic_richards()    

    def compute_algorithm_points_richards(self):
        if self.R == True:
            self.mean_fake_well.compute_algorithm_points_richards()

    def plot_richards(self):
        if self.R == True:
            self.mean_fake_well.plot_richards()

    def compute_activity(self, reporter_dataframe, use_model_data=True):
        if self.R == True:

            for w in self.wells_dict:
                w_obj = self.wells_dict[w]
                w_obj.compute_Richards()
                w_obj.kinetic_richards()
                w_obj.compute_algorithm_points_richards()
                w_obj.compute_activity(reporter_dataframe, use_model_data)
            
            self.reporter_data = pd.DataFrame([self.wells_dict.get(w).reporter_data for w in self.wells], index=self.wells).T
            self.reporter_data_mean = self.reporter_data.mean(axis=1)
            self.reporter_data_std = self.reporter_data.std(axis=1)
    
            self.mean_fake_well.compute_activity(pd.DataFrame(self.reporter_data_mean, columns=[self.name]), use_model_data)


    def plot(self, reporter_plot=False):
        if self.R == True:
            plt.rcParams.update({'font.size': 10})
            if reporter_plot:
                fig, ([ax1, ax2], [ax3, ax4]) = plt.subplots(2,2,figsize=(12, 8))
            else:
                fig, (ax1, ax2) = plt.subplots(1,2,figsize=(12, 4))
                
            #Replicates growth
            ax1.set_title('Replicates normalized growth')
            ax1.plot(self.time, self.data, label=self.data.columns)
            ax1.legend()
    
            #Mean + std and model fit
            ax2.set_title('Mean + std')
            ax2.plot(self.time, self.data_mean, '.', c=self.color, label='Data mean', ls='', ms=6)
            ax2.fill_between(self.time, self.data_mean+self.data_std, self.data_mean-self.data_std, alpha=0.2, color=self.color)
            #ax2.axvspan(self.mean_fake_well.alg_idx[0], self.mean_fake_well.alg_idx[1], color='black', alpha=0.2)
            a, b = self.mean_fake_well.alg_idx
            ax2.plot(self.mean_fake_well.exponential, self.data_mean[a:b], marker='o', c='black', alpha=0.5, fillstyle='none', ls='', ms=10, label='Exp. phase')
            ax2.plot(self.mean_fake_well.fit_x, self.mean_fake_well.fit_y, '--', c='black', alpha=0.5, label='Model fit')
            ax2.legend()
    
            if reporter_plot:
                #Growth vs. reporter
                ax3.set_title('Growth vs. reporter')
                ax3.plot(self.data_mean, self.reporter_data_mean, '.', color=self.color)
                ax3.plot(self.data_mean[a:b], self.reporter_data_mean[a:b], marker='o', c= 'black', alpha=0.5, fillstyle='none', ls='', ms=10, label='Exp. phase')
                ax3.set_title('Growth vs. reporter')
                ax3.text(0.8, 0.1,('rsq: %3.2f' %self.mean_fake_well.RG_rsq), transform=ax3.transAxes,
                             bbox=dict(facecolor='w', alpha=0.5), va='center', ha='center')
                xx = np.linspace(list(self.data_mean)[0], list(self.data_mean)[-1])
                ax3.plot(xx, linear(xx, self.mean_fake_well.RG_reg.slope, self.mean_fake_well.RG_reg.intercept), '--', c='black', alpha=0.5)
    
                #Dynact plot
                ax4.set_title('Dynamic activity')
                ax4.plot(self.time, self.mean_fake_well.dynact, c=self.color)
                ax4.plot(self.mean_fake_well.exponential, self.mean_fake_well.dynact[a:b], marker='o', c= 'black', alpha=0.5, fillstyle='none', ls='', ms=10)
                ax4.axhline(y=0, c='black', ls='--', alpha=0.2)

        
                                   