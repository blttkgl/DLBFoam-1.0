import numpy as np
import pandas as pd
import re
import os
import pylab as pl
import sys
import pdb
def main():

    cases = ['core_20','core_40','core_80','core_160','core_320']
    plot_bar(cases)
    

def plot_bar(cases):
    fig,ax = pl.subplots()
    width = 0.15  # the width of the bars
    i = 1
    for case in cases:
        df = read_data(case)
        Sr = get_Sr(df)
        x = np.arange(len(Sr))  # the label locations
        rects1 = ax.bar(x - width*i, Sr, width, label='N$_{CPU}$: ' + re.findall('\d+', case )[0])
        i+=1
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Speed-up Ratio [-]')
    ax.set_xticks([-0.5,2.5])
    ax.set_xticklabels(['Very Unbalanced','Very Balanced'])
    ax.legend(frameon=False)
    fig.savefig('Sr.png',bbox_inches='tight')
    pl.show()


def get_Sr(df):

    grouped = df.groupby(['Function','Init. Condition'])



    Sr = []    
    x = []
    for key, item in grouped:

        func = item['Function'].unique()[0]
        cond = item['Init. Condition'].unique()[0]
        x.append(cond)
        #balancer = item['Balancer'].unique()[0]
        title = "{0}_{1}".format(func, cond)




        print(item)
        models = item["Model"].unique()
        n_models = len(models)

        print(models)



        for i, model in enumerate(models):

            sliced = item.loc[item["Model"] == model]

            #error =
            if(i==0):
                Srtmp=np.max(sliced["Mean"])
            else:
                Srtmp = Srtmp/np.max(sliced["Mean"])

        Sr.append(Srtmp)

    return Sr
        

def read_data(case):

    names = get_fnames(case)
    frames = []
    ids = []
    for name in names:
        
        proc_id = int(re.search(r'\d+', name).group())
        df = pd.read_csv(case +'/' +  name)
        ddf = df.assign(Processor=proc_id)
        frames.append(ddf)
        

    big = pd.concat(frames, ignore_index=True)
    return big


def get_fnames(case):

    all_files = os.listdir(case+'/')
    good_files = []

    for f in all_files:
        if ("results_" in f):
            good_files.append(f)
    return good_files






main()
