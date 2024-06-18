import numpy as np
import pandas as pd

#function to check that elements do not repeat within a list
def check_unique(query_list):
    if len(query_list) != len(set(query_list)):
        return True

# function to give the correct format to dataframes: wells in columns and timepoints for the indices. Also create a global time list
def correct_shape(dataframe, bouts):
    first = dataframe.columns[0]
    dataframe.index = list(dataframe[first])
    dataframe = dataframe.iloc[:, 1:]
    f = len(dataframe.columns) * bouts
    global time
    time = np.arange(0, f, bouts)
    dataframe.columns = time
    dataframe = dataframe.transpose()
    return dataframe

#function to get midpoint of a list
def len_list(list_data):
    l = len(list_data)
    if l%2 != 0:
        idx = int(l/2)
        return list_data[idx]
    else:
        idx1 = int(l/2)-1
        idx2 = int(l/2)
        return((list_data[idx1]+list_data[idx2])/2)

#round a number to nearest half (this assuming we measure each 0.5, could be optimized)
def round_half(number):
        return (round(number * 2) / 2)

#round number to given bout, optimized version of round_half
def round_to_bout(number, bout):
    return (round(number * (1/bout)) * bout)

#MAIN FUNCTION FOR KEREN ET AL. ALGORITHM
    #We apply the algorithm to individual samples instead of the whole plate
    #We minimize the variance in growthrate instead of dOD/dt
    #As there is the possibility that 2*doubling times exceeds the time of actual exponential phase, we can choose how many doubling times to use

#normalize values of a list to a value of that list
def normalization(list_input, normalization_idx, value=1):
    norm_target = list_input[normalization_idx]
    list_output = []
    for data in list_input:
        if data == norm_target:
            norm_value = value
        else:
            norm_value = (data * value)/norm_target
        list_output.append(norm_value)
    return list_output

#functions to generate plots from fits
def y_fit_points(x_data, fit_data):
    deg = len(fit_data)
    y_data = 0
    for coef, exp in zip(range(0, deg), reversed(range(0, deg))):
        y_data = y_data + fit_data[coef] * x_data ** exp
    return y_data

def fit_equation(fit_data):
    deg = len(fit_data)
    eq = '0'
    for coef, exp in zip(range(0, deg), reversed(range(0, deg))):
        eq = eq  + ' + %f * x ** %d' %(fit_data[coef], exp)
    return eq

def exp_func(x, slope, intercept):
    A = np.exp(intercept)
    m = slope
    return A * np.exp(m * x)

def fit_func(x, coef1, coef2, coef3):
    return coef1 * x ** 2 + coef2 * x ** 1 + coef3 * x ** 0

#function to get sample info from the output string info
def text_to_sample(string):
    samples = []    
    open_par = string.find('(') + 1
    close_par = string.find(')')    
    while True:
        out = string[open_par:close_par]
        if 'Samples' in out:
            break
        samples.append(out)        
        open_par = string.find('(', close_par) + 1
        close_par = string.find(')', open_par)
    return samples

#match sample info to wells
def sample_to_wells(string):
    wells = []
    ap = string.find("'") + 1
    ap_next = string.find("'", ap)
    while True:
        out = string[ap:ap_next]
        if out == '':
            break
        elif ' ' in out:
            ap = string.find("'", ap_next) + 1
            ap_next = string.find("'", ap)
        else:
            ap = string.find("'", ap_next) + 1
            ap_next = string.find("'", ap)
            wells.append(out)
    return wells

#function to convert algorithm output info to actual values when the algorithm
#computes only considering some cropped timespan of the whole data
def index_crop_to_index(initial_crop_val, w_idcs_from_algorithm, time):    
    n = list(time).index(initial_crop_val)
    i = w_idcs_from_algorithm[0]
    r = w_idcs_from_algorithm[-1]
    k = list(time)[n+i]
    p = list(time).index(k)
    s = list(time)[n+r]
    j = list(time).index(s)
    w_idcs_no_crop = np.arange(p, j+1, 1)
    return w_idcs_no_crop

#list indx to X,2 shape coordinates
def convert_idx_to_coord(idx):
    row = int(idx/2)
    col = idx-2*row
    return [row, col]

wells_in_plate = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12',
                  'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'B12',
                  'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12',
                  'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'D12',
                  'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12',
                  'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12',
                  'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12',
                  'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12']

#function to get wells that have no samples in the experiment
def spot_removed_wells(trimmed_df_dict):
    df_key = list(trimmed_df_dict.keys())[0]
    wells_in_df = list(trimmed_df_dict[df_key].keys())
    wells_not_in_df = []
    for well in wells_in_plate:
        if well not in wells_in_df:
            wells_not_in_df.append(well)
    return wells_not_in_df

colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange', 'tab:cyan', 'tab:purple', 'tab:pink', 'tab:brown',
          'tab:olive', 'lightcoral', 'seagreen', 'mediumpurple', 'gold', 'sandybrown', 'deepskyblue', 'palegreen',
          'violet', 'darkturquoise', 'limegreen', 'lightsalmon']

def compute_zscores(dataframe, zscore_threshold):
    #create zscore dataframe
    dataframe_transposed = dataframe.T
    zscore_df = pd.DataFrame(index=dataframe.index, columns=dataframe.columns)
    #iterate through timepoints
    for timepoint in dataframe_transposed.columns:
        data = dataframe_transposed[timepoint]
        timepoint_zscore = []
        #iterate through each point
        for i in range(len(data)):
            target = data[i]
            data_copy = list(data).copy()
            data_copy.remove(target)
            sample_mean = np.mean(data_copy)
            sample_std = np.std(data_copy)
            individual_zscore = (abs(target - sample_mean))/sample_std
            timepoint_zscore.append(individual_zscore)
        zscore_df.loc[timepoint] = timepoint_zscore
    #select outliers
    outliers = []
    zsum = zscore_df.sum()
    zmean, zstd = np.mean(zsum), np.std(zsum)
    for element, value in zip(zsum.index, zsum):
        if value > zmean + zscore_threshold * zstd:
            outliers.append(element)
    return outliers

#given a model data, returns de closest actual data indices
def find_closest_data(xmodel, xdata):
    idcs_list = []
    for i in list(xdata):
        abs_list = [abs(x-i) for x in list(xmodel)]
        target = min(abs_list)
        idcs_list.append(abs_list.index(target))
    return idcs_list

#line equation
def linear(x, m, n):
    f = m * x + n
    return f

#Modified Richards model, from Cerillo Labrador software
def Richards(t, A, K, r, gamma, nu):
    f = A + (K-A)/((1+np.exp(-r*(t-gamma)))**(1/nu))
    return f

#Richards equation derivative
def Richards_ln_diff(t, A, K, r, gamma, nu):
    factor1 = np.exp(-r*(-gamma+t))
    numerator = -A * r * factor1 * (factor1+1)**(-1+(2/nu)) + K * r * factor1 * (factor1+1)**(-1+(2/nu))
    denom = A * nu * (factor1+1)**(1/nu) * ((factor1+1)**(1/nu))**2 - A * nu * ((factor1+1)**(1/nu))**2 + K * nu * ((factor1+1)**(1/nu))**2
    f = numerator/denom
    return f