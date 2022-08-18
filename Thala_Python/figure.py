


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

fat_plt1 = pd.read_csv("fout.txt", sep='\t')
mat_plt1 = pd.read_csv("Mout.txt", sep='\t')
x_fat1 = np.array(fat_plt1['POS'])
y_fat1 = np.array(fat_plt1['OddRatio_F'])
x_mat1 = np.array(mat_plt1['pos'])
y_mat1 = np.array(mat_plt1['OddRatio_M'])
target = [5246696,5248301]



def draw(x_mat, y_mat, x_fat, y_fat):
    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)
    plt.xlim(4.7e6,5.4e6, 0.0001)
    plt.ylim(-10,10)
    fig.supxlabel('Position')
    fig.supylabel('Odd Ratio')

    axs[0].axhline(0, linestyle=':' , color='black' )
    axs[0].plot(x_mat, y_mat, color='red', label='haplotype')
    axs[0].axvline(target[0], linestyle='--', label='HBB', color='yellow')
    axs[0].axvline(target[1], linestyle='--', color='yellow')
    axs[0].set_title("maternal inheritance")
    axs[0].text(4.55e6, 4.5, 'P', style='italic', bbox=dict(facecolor='red', alpha=0.9, boxstyle='circle'), size=40)
    axs[0].text(4.55e6, -5.5, 'N', style='italic', bbox=dict(facecolor='green', alpha=0.9, boxstyle='circle'), size=40)


    axs[1].axhline(0, linestyle=':', color='black')
    axs[1].plot(x_fat, y_fat, color='red')
    axs[1].axvline(target[0],linestyle='--', color='yellow')
    axs[1].axvline(target[1],linestyle='--', color='yellow')
    axs[1].set_title("paternal inheritance")

    fig.savefig(sys.argv[1])
    fig.legend()
   
    
    
  
    plt.show()

draw(x_mat1,y_mat1,x_fat1,y_fat1)