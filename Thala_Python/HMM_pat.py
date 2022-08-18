
import pandas as pd
import numpy as np
from scipy.stats import binom 
import math

#Lay nhung cap di hop bo, dong hop me tu file out.txt
read_out = pd.read_csv('out.txt', sep='\t')
for i in range (len(read_out)):
    if(read_out.at[i, 'F0']!=read_out.at[i, 'F1'] and 
    read_out.at[i, 'M0']==read_out.at[i, 'M1']):
        continue
    else:
        read_out.drop(i, inplace=True)
read_out.reset_index(inplace=True, drop=True)

gen_map = pd.read_csv('genetic_map_GRCh37_HBB_HBA_1M.txt', sep='\t', header=None)
#Neu chi co 1 loai NST 11 hoac 16
if (len(read_out.value_counts(['CHROM'])))==1:
    chr_type = read_out.at[0,'CHROM']
    for i in range(len(gen_map)):
        if(gen_map.at[i,0]!=chr_type):
            gen_map.drop(i, inplace=True)

gen_map.drop([0, 2], axis=1, inplace=True)
read_out['NAN'] = np.nan

fra_mean = float(open("ff.txt", 'r').read())

# Function to caculate Transmission probability
def getTransProb(genMap, snpSet, chrName): # DataFrame, DataFrame, String
    
    i0, i1 = 0, 0
    while(i1<len(snpSet) and snpSet.at[i1,'POS']<=genMap.at[0,1]):
        snpSet.at[i1, 'NAN']=0
        i1+=1
    while(i1<len(snpSet)):
        while(i0<len(genMap) and genMap.at[i0,1]<snpSet.at[i1,'POS']):
            i0+=1

        if(i0<len(genMap)):
            snpSet.at[i1,'NAN'] = ((genMap.at[i0, 3]-genMap.at[i0-1, 3])*
                (snpSet.at[i1,'POS']-genMap.at[i0-1,1]) / (genMap.at[i0,1]-genMap.at[i0-1,1]) 
                + genMap.at[i0-1,3])
            i1+=1
        else:
            while(i1<len(snpSet)):
                snpSet.at[i1,'NAN'] = genMap.at[len(genMap)-1,3]
                i1+=1
        Trans_Prob = []
    
        for i in range(len(snpSet)-1):
            p=((snpSet.at[i+1,'NAN']-snpSet.at[i,'NAN'])/100)
            
            if(p>0.5):
                p=0.5
            Trans_Prob.append(np.array([1-p, p, p, 1-p]).reshape(2,2))
    return Trans_Prob
                    

## Calculate emission probability 
uplimit=0.95
lowlimit=0.05
def getEmisProb(snpSet):
    #P0 P1 are emission probability
    snpSet['Hom'] = np.nan
    snpSet['Het'] = np.nan
    snpSet['P0'] = np.nan
    snpSet['P1'] = np.nan
    Emiss_Pro = []
    for i in range(len(snpSet)):
        depth = snpSet.at[i, 'Plasma_REF'] + snpSet.at[i, 'Plasma_ALT']
        alt = snpSet.at[i,'Plasma_ALT']
        ref = snpSet.at[i,'Plasma_REF']

        if snpSet.at[i, 'M0']==1:
            snpSet.at[i, 'Hom'] = binom.pmf(ref, depth, 0.02)
            snpSet.at[i, 'Het'] = binom.pmf(ref, depth, 0.5*fra_mean)
            if snpSet.at[i,'F0'] == 1:
                snpSet.at[i,'P0'] = snpSet.at[i, 'Hom']/(snpSet.at[i, 'Hom']+snpSet.at[i, 'Het'])
            else:
                snpSet.at[i,'P0'] = snpSet.at[i, 'Het']/(snpSet.at[i, 'Hom']+snpSet.at[i, 'Het'])
        
        elif snpSet.at[i, 'M0']==0:
            snpSet.at[i, 'Hom'] = binom.pmf(alt, depth, 0.02)
            snpSet.at[i, 'Het'] = binom.pmf(alt, depth, 0.5*fra_mean)
            if snpSet.at[i,'F0'] == 0:
                snpSet.at[i,'P0'] = snpSet.at[i, 'Hom']/(snpSet.at[i, 'Hom']+snpSet.at[i, 'Het'])
            else:
                snpSet.at[i,'P0'] = snpSet.at[i, 'Het']/(snpSet.at[i, 'Hom']+snpSet.at[i, 'Het'])
        
        snpSet.at[i,'P0'] = min(snpSet.at[i,'P0'], uplimit)
        snpSet.at[i,'P0'] = max(snpSet.at[i,'P0'], lowlimit)
        snpSet.at[i,'P1'] = 1 - snpSet.at[i,'P0']
        Emiss_Pro.append(np.array([snpSet.at[i,'P0'], snpSet.at[i,'P1']]))
    return Emiss_Pro

t = getTransProb(gen_map, read_out, 'chr11')
e = getEmisProb(read_out)



symbols =  read_out['CHROM']
observation = np.array(range(len(read_out)))

States=[0, 1]
startProbs=[0.5, 0.5]

def HMM(states, symbols, startProb, transProb, emissProb):
    nStates = len(states)
    nSymbols = len(symbols)
    S = [1/nStates]*nStates
    T = np.diag([0.5]*nStates) + 0.5/nStates
    E = np.zeros(shape=(nStates, nSymbols)) + 1/nSymbols
    if(startProb!=None):
        S = startProb
    if(transProb!=None):
        T = transProb
    if(emissProb!=None):
        E = emissProb
    return states, symbols, S, T, E



def viterbi(States,  startProb, transProb, emissProb, observation):
    nState = len(States)
    nObservation = len(observation)
    v = np.zeros([nState, nObservation])
    
    for i in States:
        v[i, 0] = math.log((startProb[i] * emissProb[observation[0]][i]))
        
   
    if nObservation >= 2:
        for k in range(1, nObservation):
            for state in States:
                maxi = -999999
                for preState in States:

                    temp = v[preState][k-1] + math.log(transProb[k-1][preState][state])
                    maxi = max(maxi, temp)
                    
                    
                v[state, k] = math.log(emissProb[k][state]) + maxi

    optimal_path = max(v[0, nObservation-1], v[1, nObservation-1]) 
    print(optimal_path)
    viterbiPath = np.zeros([nObservation])
    for state in States:
        if optimal_path == v[state, nObservation-1]:
            viterbiPath[nObservation-1] = state+1
            
            break
    
    if(nObservation>=2):
        for k in range((nObservation-1)-1, -1, -1): # k from n-1 -> 0
            for state in States:
                value = int(viterbiPath[k+1]-1)
               
                t0 = transProb[k][0,value]
                t1 = transProb[k][1,value]
                v0 = v[0,k]
                v1 = v[1,k]
                
                if v0 + math.log(t0) > v1 + math.log(t1):
                    viterbiPath[k] = 1
                    break
                else:
                    viterbiPath[k] = 2
                    break

    return viterbiPath
                
    

path = viterbi(States, startProbs, t, e, observation)
print(type(path))
read_out.drop(columns='NAN', inplace=True)
read_out['Fetal_F'] = path


#Odd Ratio
OR = np.zeros([len(read_out)])

iRow = int(path[0])-1
transP0 = t[0][iRow, 0]
transP1 = t[0][iRow, 1]

emissP0 = read_out.at[0, 'P0']
emissP1 = read_out.at[0, 'P1']

OR[0] = math.log10(transP0*emissP0/(transP1*emissP1))

for i in range(0,len(read_out)-1):
    iRow = int(path[i])-1
    transP0 = t[i][iRow, 0]
    transP1 = t[i][iRow, 1]
    emissP0 = read_out.at[i+1, 'P0']
    emissP1 = read_out.at[i+1, 'P1']
    OR[i+1] = math.log10(transP0*emissP0/(transP1*emissP1))


read_out['OddRatio_F'] = OR
read_out.to_csv('fout.txt', index=False, sep='\t')

