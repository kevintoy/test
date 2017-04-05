# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 01:55:13 2017

@author: kevinhong
"""

from scipy.stats.stats import pearsonr
import math
import numpy as np
import random
import matplotlib.pyplot as plt
group_size=300
generation=100
beta=0.3#strength of genetic influence in gene_cul_type, can be thought of as the heritability of 
#educational attainment
beta_1=0.15 #coefficient in the genetic/cultural interaction in determining the phenotype
phi=1 #degree of genetic reproductive skew
phi_1=1 #degree of cultural reproductive skew
phi_1_shock=3
gamma_cul=0.2 #relative contribution of geno_cultural type in genetic fitness
gamma_cul_shock=0.8
gen_shock=20

beta_gen=1.0 #degree of "evolutionary override"
genetic_corr=0.9
gene_mut_freq=0.0
gene_mut_mag=0.0
cul_mut_freq=0.3
cul_mut_mag=0.2


repeat=5
def plot_tom():
    plt.plot(tom_stat(pop),pop_stat(pop,0),"ro")
def plot_fit():
    plt.plot(pop_stat(pop,1),pop_fitness(pop),"go")
def genetic_fitness(ind,pop):
    #the trait that is controlled by the locus but positively contributed to fitness#
    tom=np.random.normal(ind[0],3*(1-genetic_corr))     #3 is arbitrary here.
    if tom<0:
        tom_score1=0.00000001+random.uniform(-0.2,0.2)
    elif tom>1:
        tom_score1=0.9999999+random.uniform(-0.2,0.2)
    else:
        tom_score1=tom+random.uniform(-0.2,0.2)
    fitness=2+(1-gamma_cul)*tom_score1-gamma_cul*ind[1]+random.uniform(-0.5,0.5)
    if fitness<0:
        print "fitness smaller than zero"
        fitness=0.0
    
    return fitness**phi
    
def genetic_fitness_shock(ind,pop):
    tom=np.random.normal(ind[0],3*(1-genetic_corr))     #3 is arbitrary here.
    if tom<0:
        tom_score1=0.00000001+random.uniform(-0.2,0.2)
    elif tom>1:
        tom_score1=0.9999999+random.uniform(-0.2,0.2)
    else:
        tom_score1=tom+random.uniform(-0.2,0.2)
    fitness=2+(1-gamma_cul_shock)*tom_score1-gamma_cul_shock*ind[1]+random.uniform(-0.5,0.5)
    if fitness<0:
        print "fitness smaller than zero"
        fitness=0.0
    
    return fitness**phi


def tom_stat(pop):
    z=[]
    for i in pop:

        tom_score=np.random.normal(i[0],3*(1-genetic_corr))
        if tom_score<0:
            z.append(0.0000001+random.uniform(0,0.2))
        elif tom_score>1:
            z.append(0.99999999+random.uniform(-0.2,0))
        else:
            z.append(tom_score)
        
    return z
    
def pop_stat(pop,position):
    x=[]
    for i in pop:
        x.append(i[position])
    return x

def pop_fitness(pop):
    y=[]
    for i in pop:
        y.append(genetic_fitness(i,pop))
    return y
def pheno_vs_fitness():
    
    plt.plot(pop_stat(pop,1),pop_fitness(pop),"ro")
    plt.show()
    
    
def pop_mean(pop,position):
    total=0
    for i in pop:
        total=total+i[position]
    return total/len(pop)
        
    
def cultural_fitness(ind):
    return (ind[1])**phi_1

def cultural_fitness_shock(ind):
    return (ind[1])**phi_1_shock
    
def geno_cul_type(genotypic_value,pop,generation):
    ind_cul_fitness_list=[]    
    for ind in pop:
        
        ind_cul_fitness_list.append(cultural_fitness(ind))
    p_fit_weight_cul=[float(i)/sum(ind_cul_fitness_list) for i in ind_cul_fitness_list]

    cul_model_index=np.random.choice(range(len(pop)),1,replace=True,p=p_fit_weight_cul)
    
    return ((genotypic_value*beta)+(1-beta)*(pop[cul_model_index][1])+beta_1*genotypic_value*pop[cul_model_index][1])  *  (beta_gen**generation)
    
   #just taking the average of genotypic and model geno_cul_type value
def geno_cul_type_shock(genotypic_value,pop,generation):
    ind_cul_fitness_list=[]    
    for ind in pop:
        
        ind_cul_fitness_list.append(cultural_fitness_shock(ind))
    p_fit_weight_cul=[float(i)/sum(ind_cul_fitness_list) for i in ind_cul_fitness_list]

    cul_model_index=np.random.choice(range(len(pop)),1,replace=True,p=p_fit_weight_cul)
    
    return ((genotypic_value*beta)+(1-beta)*(pop[cul_model_index][1])+beta_1*genotypic_value*pop[cul_model_index][1])  *  (beta_gen**generation)

ave_geno=[]
ave_pheno=[]
for i in range(repeat):

    pop=[]
    for i in range(group_size):
        pop.append([round(random.random(),2),round(random.random(),4)])
        
        
    
    #now construct the f1 generation
    time_serie_geno=[]
    time_serie_pheno=[]


    gen=1
    for i in range(generation):
    
        ind_fitness_list=[]    
        for ind in pop:
            if gen<gen_shock:
                ind_fitness_list.append(genetic_fitness(ind,pop))
            else:
                ind_fitness_list.append(genetic_fitness_shock(ind,pop))
        
        
        p_fit_weight=[float(i)/sum(ind_fitness_list) for i in ind_fitness_list]
            
       
        
        
        
        parents_index=np.random.choice(range(len(pop)),group_size,replace=True,p=p_fit_weight) 
        #now we have the indices of individuals in the f1 generation that will reproduce
        naive_offspring_set=[]
                    
                    
        for i in parents_index:
            naive_offspring_set.append(pop[i])
        # the second element in each naive_offspring doesn't matter, we'll override it
        
        pop_f1=[]
        for i in naive_offspring_set:
            if random.random()<gene_mut_freq:
                if random.random()<0.5:
                    new_geno=i[0]*(1-gene_mut_mag)
                else:
                    new_geno=i[0]*(1+gene_mut_mag)
                
            else:
                new_geno=i[0]
            if gen<gen_shock:    
                if random.random()<cul_mut_freq:
                    if random.random()<0.5:                
                        new_cul=(geno_cul_type(i[0],pop,gen))*(1+cul_mut_mag)
                    else:
                        new_cul=(geno_cul_type(i[0],pop,gen))*(1-cul_mut_mag)
                    
                else:
                    new_cul=geno_cul_type(i[0],pop,gen)
            
            else:            
                if random.random()<cul_mut_freq:
                    if random.random()<0.5:                
                        new_cul=(geno_cul_type_shock(i[0],pop,gen))*(1+cul_mut_mag)
                    else:
                        new_cul=(geno_cul_type_shock(i[0],pop,gen))*(1-cul_mut_mag)
                    
                else:
                    new_cul=geno_cul_type_shock(i[0],pop,gen)
            
            
            
            pop_f1.append([round(new_geno,4),round(new_cul,4)])
               
       
                
        
        pop=pop_f1
        gen=gen+1
        print round(pop_mean(pop,0),4), round(pop_mean(pop,1),4)
        time_serie_geno.append(pop_mean(pop,0))
        time_serie_pheno.append(pop_mean(pop,1))
    plt.plot(range(1,generation+1),time_serie_geno,"r")
    plt.plot(range(1,generation+1),time_serie_pheno,"b")
    plt.ylim(0,5)
    ave_geno.append(pop_mean(pop,0))
    ave_pheno.append(pop_mean(pop,1))
print "average_geno=", sum(ave_geno)/repeat
print "average_pheno=", sum(ave_pheno)/repeat