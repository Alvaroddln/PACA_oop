import os
import numpy as np
import pandas as pd
from functions import text_to_sample, sample_to_wells

###FUNCTION TO CREATE DATA DICTS###
def create_data_dicts(datasets_dictionary, sample_info_file_name):
    #Assigning new paths to each data file and assigning
    #the corresponding data structure for further analysis
    DF_def = {}
    dfs_raw_names = list(datasets_dictionary)

    for dataset in datasets_dictionary:
        dataset_path = '/content' + '/' + datasets_dictionary[dataset]
        DF_def[dataset] = pd.read_csv(dataset_path, index_col=0)

    samples_path = '/content' + '/' + sample_info_file_name
    f = open(samples_path, 'r')
    lines = f.readlines()
    f.close()
    samples_input_str = str(lines)

    #defining time and bouts
    time = np.array(DF_def[dataset].index)
    bouts = time[1]-time[0]

    print('INPUT SUMMARY:')
    print('--------------\n')
    print('File with %d datasets: ' %len(dfs_raw_names), dfs_raw_names)
    print('\nbouts = %.2f (h\u207B\u00B9)' %bouts)

    #Creating the sample-well strcuture
    samples = text_to_sample(samples_input_str)

    sample_names = []
    samples_def = {}
    blanks_def = {}
    wells_def = []

    for line in samples:
        wells = sample_to_wells(line)
        wells_temp =[]
        for i, well in enumerate(wells):
            if i == 0:
                sample_name = well
                if sample_name == 'blank':
                    pass
                else:
                    sample_names.append(sample_name)
            else:
                wells_temp.append(well)
                wells_def.append(well)
        if sample_name == 'blank':
            blanks_def.update({sample_name : wells_temp})
        else:
            samples_def.update({sample_name : wells_temp})
            
    #definition of wells without blanks
    wells_def_wo_blank = []
    for s in samples_def:
        for w in samples_def[s]:
            wells_def_wo_blank.append(w)
            
    print('\nSamples:\n')
    print(*samples_def.items(), sep='\n')
    print('\nBlanks:\n')
    print(*blanks_def.items(), sep='\n')

    return time, bouts, DF_def, sample_names, samples_def, blanks_def, wells_def, wells_def_wo_blank
