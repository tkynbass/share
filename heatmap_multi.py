
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from subprocess import getoutput
import itertools as itr
import random
import multiprocessing as mp
import sys

def Plot_hist2d (array, filename, color):
    
    fig, ax1 = plt.subplots ()
    
    plt.xlabel ('N-S Distance (um)', fontsize=20)
    plt.ylabel ('Distance between pos and x-axis', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=10)
    
    ax1.hist2d (array[0], array[1], alpha=0.7, cmap=color, bins=30)
    
    plt.savefig (filename, transparent=True, bbox_inches='tight')
    
    plt.close()

def subprocess (idx, send_rev, df, high_part, low_part, NS, nuc_pos):

    df.columns = ['chr', 'X', 'Y', 'Z']
    
    sub_high_pos = []
    sub_low_pos = []
    sub_chr_pos = [ [] for chr in range (3) ]
    
    pos_list = []
    
    LENGTH= 0.0966
    
    for index, row in df.iterrows ():
        
        NG = row['X':'Z'] - nuc_pos
        NG_norm = np.linalg.norm (NG)
        theta = np.arccos (np.dot (NS, NG) / ( np.linalg.norm (NS) * NG_norm ) )
        pos_list.append  ([NG_norm * np.cos(theta), NG_norm * np.sin (theta)])
    
    df = pd.concat ( [df, pd.DataFrame (pos_list, columns=['Proj_X', 'Proj_Y'])], axis=1)

    sub_high_pos.extend (df.iloc [high_part, 4:6].values * LENGTH)
    sub_low_pos.extend (df.iloc [low_part, 4:6].values * LENGTH)
    
    for chr_no in range (3):
        sub_chr_pos [chr_no].extend (df.loc [df.loc[:, 'chr']==chr_no, 'Proj_X':'Proj_Y'].values * LENGTH)

    send_rev.send ([sub_high_pos, sub_low_pos, sub_chr_pos])

#### main ####

def main ():
    
    argvs = sys.argv
    dir = argvs[1]
    status = int (argvs[2])
    
    low_gene = pd.read_csv ('low_gene_5k_center.txt', sep='\s+')
    high_gene = pd.read_csv ('high_gene_5k_center.txt', sep='\s+')

    high_part = high_gene.loc [:, 'Center'].values
    low_part = low_gene.loc [:, 'Center'].values

    file_list = getoutput (f'ls {dir}/sample_*.txt').split()
    df_list = [ pd.read_csv (file, header=None, usecols=[1,3,4,5], sep='\s+') for file in file_list ]
    
    stable_df = pd.read_csv (f'stable_status.txt', index_col=0, sep='\s+')
    nuc_pos = stable_df.iloc[status, 0:3]

    stable_spb = pd.read_csv (f'stable_spb.txt', index_col=0, sep='\s+')
    spb_pos = stable_spb.iloc[status, 0:3]

    NS = spb_pos - nuc_pos

    high_pos = []
    low_pos = []
    chr_pos = [ [] for chr in range (3) ]

    process_list = []
    pipe_list = []

    for idx in range (len (file_list)):
    
        get_rev, send_rev = mp.Pipe (False)
        p = mp.Process (target=subprocess, args=(idx, send_rev, df_list[idx], high_part, low_part, NS, nuc_pos) )
        process_list.append (p)
        pipe_list.append (get_rev)
        p.start ()
    
    for proc in process_list:
        proc.join ()
    
    for x in pipe_list:
        
        high_pos.extend (x.recv()[0])
        low_pos.extend (x.recv()[1])

        for chr_no in range (3):
            chr_pos [chr_no].extend (x.recv [2][chr_no])

    cmap_list = ['Reds', 'Blues', 'Greens']

    high_array = np.array (high_pos).T.reshape (2, len (high_pos))
    low_array = np.array (low_pos).T.reshape (2, len (low_pos))

    chr_arrays = [ np.array (pos).T.reshape (2, len (pos)) for pos in chr_pos ]

    Plot_hist2d (high_array, f'{dir}/hg_map.png', 'Reds')
    Plot_hist2d (low_array, f'{dir}/lg_map.png', 'Reds')

    for chr_no in range (3):
        Plot_hist2d (chr_arrays [chr_no], f'{dir}/chr{chr_no + 1}_map.png', cmap_list[chr_no])

    return 0

if __name__ == "__main__" :
    main()
