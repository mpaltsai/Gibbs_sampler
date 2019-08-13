#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 13:39:27 2019

@author: ismini
"""

from random import randint
from collections import Counter, defaultdict
import pandas as pd 
import numpy as np 
import os


def Gibbs_assignment():

    def profile_matrix(list_motifs):
            
            '''function to generate the profile matrix of a list of motifs'''
            
            profile = []
            
            pos = range(len(list_motifs[0])) #range of positions in motifs
            
            n_seqs = range(len(list_motifs)) #range of motifs
            
            n_mot = len(list_motifs) #number of motifs
            
            for i in pos: #for each position count occurences of the four bases
                counter = Counter({'A':1, 'C':1, 'G':1, 'T':1}) # apply pseudocounts (initialise with 1 instead of 0)
                
                for n in n_seqs: #do it for each motif
                    counter[list_motifs[n][i]] += 1
                    
                freqs = {}
                for nuc in counter.keys(): #compute frequences for each base
                    freqs[nuc] = counter[nuc]/n_mot
                    
                profile.append(freqs)
                
            return profile
    
    def all_kmers(sequence, k):
        
        '''Generate all kmers from sequence'''
        
        pos_s = len(sequence)-k+1
        kmers = []
        for i in range(pos_s):
            kmer = sequence[i:i+k]
            kmers.append(kmer)
        return(kmers)
    
        
    def prob_kmer(kmer, profile_matrix):
        
        '''compute probability of a kmer using a profile matrix''' 
        
        ###initialize probability
        prob_kmer=1
        
        ###take probability as proposed from the profile matrix for each base at each position of the kmer
        l_kmer = len(kmer)
        
        for i in range(l_kmer):  
            prob = profile_matrix[i][kmer[i]]
            
            ###probability in total
            prob_kmer *= prob
            
        return prob_kmer
    
    def consensus(list_seq):
        
        '''find a consensus motif from a list of motifs'''
        
        consensus_list = []
        
        n_pos = len(list_seq[0])
        
        n_list = len(list_seq)
        
        for i in range(n_pos): #for each position count occurences of the four bases
            counter = Counter({'A':1, 'C':1, 'G':1, 'T':1})
            for n in range(n_list): #do it for each motif
                counter[list_seq[n][i]] += 1
            
            ###append the base with maximum frequency to consensus list    
            base = [x for x,y in counter.items() if y==max(counter.values())] 
    
            consensus_list.append(base)
            
            ###construct consensus motif from consensus_list
            ###np.random.randint(0,len(x)) was added so as to pick randomly a base 
            ###(some bases have the same frequency, resulting in more than one maximum)
            
            consensus_kmer = ''.join([x[np.random.randint(0,len(x))] for x in consensus_list])
        
        return consensus_kmer
    
    def dist_from_consensus(list_kmers, kmer_consensus):
        
        '''compute distance of a list of kmers from another kmer'''
        
        score = 0
       
        for kmer in list_kmers:  
            
            n_site = len(kmer)
        
            for i in range(n_site):
                if kmer[i] == kmer_consensus[i]:
                    continue
                else:
                    score += 1    
        return score
    
    
    
    def gibbs_motif_disc(sequences, k):
        
        '''Gibbs sampling algorithm'''
        
        ###generate randomly a motif from each sequence and store in motifs
        motifs=[]
        l = len(sequences)
        l_s = len(sequences[0])
        for i in range(l):
            start_kmer = randint(0, (l_s-k))
            kmer = sequences[i][start_kmer : start_kmer+k]
            motifs.append(kmer)
        #motifs
        
        ###find consensus motif from motifs
        consensus_old = consensus(motifs)
        
        ###compute the score from consensus_old
        score_old = dist_from_consensus(motifs,consensus_old)
        
        no_change = 0
        rounds = 0
        
        while no_change < 50 : #stop iterating if there is no change in scores after 
                               # 50 consecutive times
            no_change += 1
            rounds += 1
            
            ###pick a sequence randomly and exclude the corresponding motif from 
            ###motifs list
            xcl_int = randint(0, l-1)
            
            motifs_xcl = motifs[:xcl_int] + motifs[xcl_int+1:]
            
            
            ###profile matrix for motifs after the exclusion
            profile_xcl = profile_matrix(motifs_xcl)
            
            ###compute all kmers from the sequence we chose previously to exclude 
            ###the motif
            seq_xcl = sequences[xcl_int]
            all_k = all_kmers(seq_xcl,k)
            
            ###choose the one with the maximum probability
            kmer_lik = max([(prob_kmer(x, profile_xcl), x) for x in all_k])[1]
            
            ###add it to the list with the other motifs
            new_motifs = motifs_xcl[:] ###take a copy of motifs_xcl
            new_motifs.insert(xcl_int, kmer_lik)
            
            ###find consensus motif from new_motifs
            consensus_new = consensus(new_motifs)
            
            ###compute the score from consensus_new
            score_new = dist_from_consensus(new_motifs,consensus_new)
            
            ###compare old score with new score
            if score_new < score_old :
                if consensus_new != consensus_old:
                    motifs = new_motifs
                    consensus_old = consensus_new
                    score_old = score_new
                    no_change = 0
                else:
                    continue
               
        return consensus_new
              
    
    
    def sim_Gibbs(k1, k2, n, file):
        
        '''run Gibbs sampler mutiple times for different sizes of kmers'''
        
        '''k1 and k2 define respectively lower and upper range of different sizes 
           of kmers
           n determines repetitions
           file is the file where sequences are stored, each one in a new line, 
           without ANY symbols and all nucleotides in upper case
           The output of this function is a dictionary where every item has k size
           as a key value and different motifs of k size are stored as respective 
           values'''
        
        with open(file) as f:
            data= f.read()
     
        data = data.strip()    
        data=data.split('\n')
        
        dict_consensus = defaultdict(list)
        
        for k in range(k1,k2+1):
           
            for i in range(n):
                d = gibbs_motif_disc(data, k)
                dict_consensus[k].append(d) 
                
        return dict_consensus
    
    ###download sequences or not if file already exists in directory
    if not os.path.exists('motifs_in_sequence.fa'):
        print('downloading file')
        os.system("wget https://www.dropbox.com/s/w9cpq4bwsb90j1i/motifs_in_sequence.fa")
        
        
    ###run Gibbs sampler 1000 times for k=[3,4,5,6,7]
    
    d = sim_Gibbs(3, 7, 1000, 'motifs_in_sequence.fa')
    
        
    ###Position Weight Matrices for each kmer stored in a list
        
    pwm_list = []
    for key in d.keys():
        
        d_key = profile_matrix(d[key])
        
        pwm = pd.DataFrame(data = d_key, index=range(1, key+1))
        
        a = pwm.T
        
        pwm_list.append(a)
        
    ###write PWMs to csv file 
    
    for x in pwm_list:
        
        x.to_csv('pwm.csv',  mode = 'a')
        
    ###Calculate information content for each PWM
    
    def info_content_pwm(pwm_matrix):
        
        '''Calculation of the information content for a given Position Weight Matrix'''
        
        ###matrix as np.array
        pwm_array = np.array(pwm_matrix)
        
        ###dimensions of array
        nuc, k = pwm_array.shape
        
        ###background frequency of each base
        bf = (1/k)*np.sum(pwm_array , axis=1, keepdims=True)
        
        ###information content
        IC = np.sum( pwm_array * np.log2(pwm_array) - pwm_array * np.log2(bf) )
        
        return IC
    
    def ic_pwm_list(pwm_list):
        
        '''returns information content as dictionary entries for pwms in a list'''
        
        
        ic_dict = {}
        
        index_number = 0    
           
        for x in pwm_list:
            nuc, k = pwm_list[index_number].shape
            ic_dict[k] = info_content_pwm(x)
            index_number += 1
        return ic_dict
    
    ###information content for each pwm
    the_list = ic_pwm_list(pwm_list)
    
    
    ###find for which k we have the maximum information content
    the_k = [x for x,y in the_list.items() if y==max(the_list.values())]
    
    def pick_pwm(pwm_list):
        
        '''choose the PWM from a list of PWMs with maximum information content''' 
        '''uses the_k object found earlier'''
        
        for x in range(len(pwm_list)):
            nuc, k = pwm_list[x].shape
            if k == the_k[0]:
                return pwm_list[x]
    
    def max_ic_motif(pwm):
        
        ###matrix as np.array
        pwm_array = np.array(pwm)
        
        ###dimensions of array
        nuc, k = pwm_array.shape
        
        ###background frequency of each base
        bf = (1/k)*np.sum(pwm_array , axis=1, keepdims=True)   
    
        ###information content matrix for each base at each position
        ic_array = pwm_array * np.log2(pwm_array) - pwm_array*np.log2(bf)
        ic_array = ic_array.reshape(nuc,k)
        
    
        indices = np.where(ic_array==np.max(ic_array, axis=0))
        
        pos_nuc = sorted(list(zip(indices[1],indices[0])))
        
        motif = []
        base_dict = {0:'A', 1:'C', 2:'G', 3:'T'}
        
        for x in pos_nuc:
            base = base_dict[x[1]]
            motif.append(base)
        
        return ''.join(motif), ic_array
    
    
    
    
    ###PWM with maximum information content
    the_pwm = pick_pwm(pwm_list)
    
    ###find the motif with maximum information content
    the_motif = max_ic_motif(the_pwm)
    
    print('The motif with the highest information content is {} and it is a {}mer'.format(the_motif[0], the_k))
    return the_pwm, the_motif[0], the_k, the_motif[1]
    
result = Gibbs_assignment()

#np.savetxt('p_freq_m.csv', result[0], delimiter= ',')
#np.savetxt('array.csv', result[3], delimiter= ',')

#bf = (1/result[2][0])*np.sum(result[0] , axis=1)