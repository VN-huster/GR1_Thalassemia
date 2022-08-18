import pandas as pd
import numpy as np
from scipy.stats import binom 
import math
import sys

mad = pd.read_csv("mad.txt", sep='\t')

fra = []
count = 0

for i in range(len(mad)):
    if mad.at[i, 'Genotype']=='0/1':
        fra.append(mad.at[i, 'Mat_ALT'] / ((mad.at[i, 'Mat_ALT'])+(mad.at[i, 'Mat_REF'])))
        


Q1,Q3 = np.percentile(fra , [25,75])

lower = Q1- (Q3-Q1)*1.5
upper = Q3 + (Q3-Q1)*1.5

for i in range(len(fra)-1, -1, -1):
    a = fra[i]
    if a>upper or a<lower:
        fra.pop(i)

e_ratio = np.array(fra).mean()

fra_mean = float(open("ff.txt", 'r').read())


uplimit=0.95
lowlimit=0.05

hap_file = sys.argv[3]
snpSet_mat = pd.read_csv(hap_file, sep='\t')    #Lay ca di hop me
for i in range (len(snpSet_mat)):
    if(snpSet_mat.at[i, 'M0']!=snpSet_mat.at[i, 'M1']):
        continue
    else:
        snpSet_mat.drop(i, inplace=True)
snpSet_mat.reset_index(inplace=True, drop=True)

gen_map = pd.read_csv('genetic_map_GRCh37_HBB_HBA_1M.txt', sep='\t', header=None)
#Neu chi co 1 loai NST 11 hoac 16
if (len(snpSet_mat.value_counts(['chr'])))==1:
    chr_type = str(snpSet_mat.at[0,'chr'])
    
    if(chr_type=='11'):
        chr_type='chr11'
    elif(chr_type=='16'):
        chr_type="chr16"

    for i in range(len(gen_map)):
        if(gen_map.at[i,0]!=chr_type):
            gen_map.drop(i, inplace=True)

gen_map.drop([0, 2], axis=1, inplace=True)
snpSet_mat['NAN'] = np.nan

def getEmisProb(snpSet):
    #P0 P1 are emission probability
    snpSet['Hom'] = np.nan
    snpSet['Het'] = np.nan
    snpSet['P0'] = np.nan
    snpSet['P1'] = np.nan
    Emiss_Pro = []
    for i in range(len(snpSet)):
        depth = snpSet.at[i, 'Mp_ref'] + snpSet.at[i, 'Mp_alt']
        alt = snpSet.at[i,'Mp_alt']
        ref = snpSet.at[i,'Mp_ref']

        if snpSet.at[i, 'F_origin']==0:
            snpSet.at[i, 'Hom'] = binom.pmf(alt, depth, (1-fra_mean)*e_ratio)
            snpSet.at[i, 'Het'] = binom.pmf(alt, depth, e_ratio)
            if snpSet.at[i,'M0'] == 0:
                snpSet.at[i,'P0'] = snpSet.at[i, 'Hom']/(snpSet.at[i, 'Hom']+snpSet.at[i, 'Het'])
            else:
                snpSet.at[i,'P0'] = snpSet.at[i, 'Het']/(snpSet.at[i, 'Hom']+snpSet.at[i, 'Het'])
        
        elif snpSet.at[i, 'F_origin']==1:
            snpSet.at[i, 'Hom'] = binom.pmf(alt, depth, (1-fra_mean)*e_ratio+fra_mean)
            snpSet.at[i, 'Het'] = binom.pmf(alt, depth, e_ratio)
            if snpSet.at[i,'M0'] == 1:
                snpSet.at[i,'P0'] = snpSet.at[i, 'Hom']/(snpSet.at[i, 'Hom']+snpSet.at[i, 'Het'])
            else:
                snpSet.at[i,'P0'] = snpSet.at[i, 'Het']/(snpSet.at[i, 'Hom']+snpSet.at[i, 'Het'])
        
        snpSet.at[i,'P0'] = min(snpSet.at[i,'P0'], uplimit)
        snpSet.at[i,'P0'] = max(snpSet.at[i,'P0'], lowlimit)
        snpSet.at[i,'P1'] = 1 - snpSet.at[i,'P0']
        Emiss_Pro.append(np.array([snpSet.at[i,'P0'], snpSet.at[i,'P1']]))
    return Emiss_Pro

def getTransProb(genMap, snpSet, chrName): # DataFrame, DataFrame, String
    if(len(snpSet)>2):
        i0, i1 = 0, 0
        while(i1<len(snpSet) and snpSet.at[i1,'pos']<=genMap.at[0,1]):
            snpSet.at[i1, 'NAN']=0
            i1+=1
        while(i1<len(snpSet)):
            while(i0<len(genMap) and genMap.at[i0,1]<snpSet.at[i1,'pos']):
                i0+=1

            if(i0<len(genMap)):
                snpSet.at[i1,'NAN'] = ((genMap.at[i0, 3]-genMap.at[i0-1, 3])*
                    (snpSet.at[i1,'pos']-genMap.at[i0-1,1]) / (genMap.at[i0,1]-genMap.at[i0-1,1]) 
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


States=[0, 1]
startProbs=[0.5, 0.5]
observation = np.array(range(len(snpSet_mat)))

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


t = getTransProb(gen_map, snpSet_mat, 'chr11')
e = getEmisProb(snpSet_mat)

path = viterbi(States, startProbs, t, e, observation)
snpSet_mat.drop(columns='NAN', inplace=True)
snpSet_mat['Fetal_M'] = path


#Odd Ratio
OR = np.zeros([len(snpSet_mat)])

iRow = int(path[0])-1
transP0 = t[0][iRow, 0]
transP1 = t[0][iRow, 1]

emissP0 = snpSet_mat.at[0, 'P0']
emissP1 = snpSet_mat.at[0, 'P1']


OR[0] = -math.log10(transP0*emissP0/(transP1*emissP1))

for i in range(0,len(snpSet_mat)-1):
    iRow = int(path[i])-1
    transP0 = t[i][iRow, 0]
    transP1 = t[i][iRow, 1]
    emissP0 = snpSet_mat.at[i+1, 'P0']
    emissP1 = snpSet_mat.at[i+1, 'P1']
    OR[i+1] = -math.log10(transP0*emissP0/(transP1*emissP1))

snpSet_mat['OddRatio_M'] = OR
snpSet_mat.to_csv('Mout.txt', index=False, sep='\t')
