
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from subprocess import getoutput
import itertools as itr
import random
import multiprocessing as mp
import sys

def Plot_ppv (x, dir, high, low, random):
    
    fig, ax1 = plt.subplots ()
    
    plt.xlabel ('Distance (um)', fontsize=20)
    plt.ylabel ('Frequency', fontsize=20)
    plt.xticks(fontsize=15)
    #     plt.yticks(fontsize=10)
    
    ax1.set_yticklabels([]) #y軸目盛りの削除
    
    a=1.0
    
    ax1.plot (x, random, color='black')
    ax1.fill_between (x, random, 0, alpha=a, color='skyblue' )
    
    ax1.plot (x, low, color='black')
    ax1.fill_between (x, low, 0, alpha=a, color='darkblue' )
    
    ax1.plot (x, high, color='black')
    ax1.fill_between (x, high, 0, alpha=a, color='orange' )
    
    plt.savefig (f'{dir}/tppv.png', transparent=True, bbox_inches="tight")

    plt.close()

def Plot_hist (high, low, random):
    
    fig, ax1 = plt.subplots ()
    
    BINS = 50
    weights = np.ones_like (high) / float (len (high))
    
    h_array = ax1.hist (high, bins=BINS, weights=weights)
    l_array = ax1.hist (low, bins=BINS, weights=weights)
    r_array = ax1.hist (random, bins=BINS, weights=weights)
    
    plt.close()
    
    return h_array[1], h_array[0], l_array[0], r_array[0]

def ppv_func (x):

    return - np.log ( (x - 0.1713) / 1.0965) / 0.6865

def subprocess (idx, num_proc, send_rev, df_list, high_loop_list, low_loop_list, control_loop_list):

    high_tppv_list = []
    low_tppv_list = []
    control_tppv_list = []

    for df in df_list [idx::num_proc]:

        high_pd = pd.Series ([ np.linalg.norm (df.iloc [X, :] - df.iloc [Y, :] ) * 0.0966 for X, Y in high_loop_list ])
        low_pd = pd.Series ([ np.linalg.norm (df.iloc [X, :] - df.iloc [Y, :] ) * 0.0966 for X, Y in low_loop_list ])
        control_pd =  pd.Series ([ np.linalg.norm (df.iloc [X, :] - df.iloc [Y, :] ) * 0.0966 for X, Y in control_loop_list ])
        
        high_tppv_list.append (sum (ppv_func ( high_pd [ (high_pd > 0.1713) & (high_pd < 1.2678 )]) ))
        low_tppv_list.append (sum (ppv_func ( low_pd [ (low_pd > 0.1713) & (low_pd < 1.2678 )]) ))
        control_tppv_list.append (sum (ppv_func ( control_pd [ (control_pd > 0.1713) & (control_pd < 1.2678 )]) ))

    send_rev.send ([high_tppv_list, low_tppv_list, control_tppv_list])

#### main ####

def main ():
    
    argvs = sys.argv
    
    dir = argvs[1]
    
    num_proc = 6 #purosesusuu
    
    low_gene = pd.read_csv ('low_gene_5k_center.txt', sep='\s+')
    high_gene = pd.read_csv ('high_gene_5k_center.txt', sep='\s+')

    high_part = high_gene.loc [:, 'Center'].values
    low_part = low_gene.loc [:, 'Center'].values

    file_list = getoutput (f'ls {dir}/sample_*.txt').split()
    df_list = [ pd.read_csv (file, header=None, usecols=[3,4,5], sep='\s+') for file in file_list ]

    high_rand = random.sample (set (high_part), k=80)
    low_rand = random.sample (set (low_part), k=80)
    control_rand = random.sample ( range (2516), k=80)

    high_loop_list = [ (x, y) for x, y in itr.combinations (high_rand, 2)]
    low_loop_list = [ (x, y) for x, y in itr.combinations (low_rand, 2)]
    control_loop_list = [ (x, y) for x, y in itr.combinations (control_rand, 2) ]

    high_tppv = []
    low_tppv = []
    control_tppv = []

    process_list = []
    pipe_list = []

    for idx in range (num_proc):
    
        get_rev, send_rev = mp.Pipe (False)
        p = mp.Process (target=subprocess, args=(idx, num_proc, send_rev, df_list, high_loop_list, low_loop_list, control_loop_list) )
        process_list.append (p)
        pipe_list.append (get_rev)
        p.start ()
    
    for proc in process_list:
        proc.join ()
    
    result_list = [ np.ravel ([ x.recv()[type] for x in pipe_list ]) for type in range (3) ]
    tppv_array = np.array (result_list).reshape (3, len (result_list))

    bin_class,  h_array, l_array, r_array = Plot_hist (tppv_array[0], tppv_array[1], tppv_array[2])
    delta = bin_class[1] - bin_class[0]
    x = [ bin_class[num] + delta  for num in range (len (bin_class) -1) ]

    Plot_ppv (x, dir, h_array, l_array, r_array)

    return 0

if __name__ == "__main__" :
    main()
