import numpy as np
import pandas as pd
import re
import os
import pylab as pl
import sys
import pdb
def main():



    df = read_data()

    #print df

    #print df.index
    #print df.xs(0, level=0)

    plot_all(df)



def plot_all(df):


    
    grouped = df.groupby(['Function','Init. Condition'])
    
    
    
    Sr = []    
    x = []
    for key, item in grouped:
        
        
        
        
        func = item['Function'].unique()[0]
        cond = item['Init. Condition'].unique()[0]
	x.append(cond)
        #balancer = item['Balancer'].unique()[0]
        title = "{0}_{1}".format(func, cond)
        
        
        
        
        print item
        models = item["Model"].unique()
        n_models = len(models)

        print models
        #fig, axs = pl.subplots(max(2, n_models), 1, sharex=True, sharey=True)
        
        
        
        for i, model in enumerate(models):
            
            sliced = item.loc[item["Model"] == model]
        
            #error = 
            if(i==0):
                Srtmp=np.max(sliced["Mean"]) 
            else:
                Srtmp = Srtmp/np.max(sliced["Mean"])
            #axs[i].bar(sliced["Processor"], sliced["Mean"], yerr = sliced["Stdev"], label=model)
        
            #axs[i].legend(loc="best")
        Sr.append(Srtmp)

        #pl.title(title)
        #pl.savefig("{0}.png".format(title))
        #pl.show()
    fig,ax = pl.subplots(1)
    pl.grid()
    ax.bar(np.arange(len(Sr)),Sr)
    ax.set_xticks(np.arange(len(Sr)))
    x = ['$\omega$ =' +  str(i) for i in x]
    ax.set_xticklabels(x)
    ax.set_ylabel('Speed-up Ratio')
    pl.show()        
        
      
    

def plot_case(df):
    
    key = "Init. Condition"
    
    conditions = df[key]
    ranks = df.index.levels[0]

    for condition in conditions:
        
        i_df = df.loc[df[key] == condition]
        rank = i_df.index.levels[0]
        mean = i_df["Mean"]
        
        #print df.loc[df[key] == condition]
        

def read_data():

    names = get_fnames()

    frames = []
    ids = []
    for name in names:
        
        proc_id = int(re.search(r'\d+', name).group())
        
        df = pd.read_csv(name)
        ddf = df.assign(Processor=proc_id)
        frames.append(ddf)
        

    big = pd.concat(frames, ignore_index=True)
    return big


def get_fnames():

    all_files = os.listdir("./")
    good_files = []

    for f in all_files:
        if ("results_" in f):
            good_files.append(f)
    return good_files






main()