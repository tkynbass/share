
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import sys

argvs = sys.argv
argc = len(argvs)

input_file = '%s' % argvs[1]
output_file = '%s' % argvs[2]

data=pd.read_table( input_file , delim_whitespace=True)

x=data.iloc[:,0]
high_y=data.iloc[:,1]
low_y=data.iloc[:,2]
ran_y=data.iloc[:,3]

#plt.figure(figsize=(5.5,4))

font=20
fig, ax =plt.subplots()

plt.xlabel("Physical proximity value",fontsize=font)
plt.ylabel("Frequency", fontsize=font)
plt.tight_layout()

a=1.0

plt.xlim (5880, 8250)
plt.ylim(0,70)

plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().yaxis.set_ticks_position('left')
plt.gca().xaxis.set_ticks_position('bottom')

plt.xticks(fontsize=13)
ax.set_yticklabels([]) #y軸目盛りの削除

plt.plot(x, ran_y, color='black')
plt.fill_between(x, ran_y, 0, alpha=a, color='skyblue')

plt.plot(x, low_y, color='black')
plt.fill_between(x, low_y, 0, alpha=a, color='darkblue')

plt.plot(x, high_y, color='black')
plt.fill_between(x,high_y, 0, alpha=a, color='orange')

plt.savefig (output_file, transparent=True)
#plt.show()

