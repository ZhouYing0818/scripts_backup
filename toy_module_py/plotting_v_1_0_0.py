# Date: 20231101
# Version: V1.0.0
# Author: ZHOU YING
# Description: This scripts contained plotting function in python
# Function List:
# 1. stack_bar_plot: plotting percentage/normal stack bar/hbar plot
#
# Modification explanation (format:{Date: stuff}):
#

# -----------------------------------------------Stack bar/hbar plot (percentage/normal)------------------------------------------------------#

def stack_bar_plot(data,num,xlabel,ylabel,savename,percentage=True,direction='v',barwidth=1,figsize='10,10'):
    '''Description: inputting data should form as pd.crosstable() results,like:
        data=pd.crosstab(adata.obs[['leiden','Lib']].leiden,adata.obs[['leiden','Lib']].Lib, margins=True)
        if need plotting percentage: data=data.div(data['All'],axis=0)
    * dataframe：columns are fraction, rows are clusters;
    * num：0：num are range of fraction;
    * xlabel, ylabel: setting x/y label;
    * savename: prefix of savename;
    * percentage: bool, default="True" plotting pecentage stack barplot, else normal form;
    * direction: {'v'/'h'}, default='v',v-vertical, h-horizontal;
    * barwidth: controlling bar width default=1;
    * figsize: a 'str' type inputs 2 numbers separate by comma, first number is figwidth parameter and second is height, default='10,10';
    '''
    import seaborn as sns
    import numpy as np
    import matplotlib.pyplot as plt
    # deciding percentage or normal barplot
    if percentage:
        # plotting percentage barplot
        # deciding horizotal or vertical barplot 
        if direction=='h':
            current_palette=sns.color_palette("coolwarm",num)
            params = {
                'figure.figsize': figsize
            }
            plt.rcParams.update(params)
            left = len(data) * [0]
            for i in range(0,num):
                plt.barh(data.index,data.iloc[:,i],left=left,
                         height=barwidth,
                         label=data.columns[i],
                         color=current_palette[i],
                         edgecolor=None,)
                left=left+data.iloc[:,i]
            plt.tick_params(axis='x',length=0)
            plt.xlabel(xlabel,fontdict={'size': 16})
            plt.ylabel(ylabel,fontdict={'size': 16})
            plt.yticks(size = 12)
            plt.xticks(size = 12)
            plt.ylim(bottom=-1)
            plt.xlim(left=0)
            plt.xticks(np.arange(0,1.2,0.2),[f'{i}%' for i in range(0,120,20)])
            plt.legend(frameon=False,fontsize=16,bbox_to_anchor=(1,1))
            plt.tight_layout()
            ax=plt.gca()
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            plt.savefig(savename+'.pdf')
        elif direction=='v':
            # plotting vertical barplot
            current_palette=sns.color_palette("coolwarm",num)
            params = {
                'figure.figsize': figsize
            }
            plt.rcParams.update(params)
            bottom = len(data) * [0]
            for i in range(0,num):
                plt.bar(data.index,data.iloc[:,i],bottom=bottom,
                        width=barwidth,
                        label=data.columns[i],
                        color=current_palette[i],
                        edgecolor=None,)
                bottom=bottom+data.iloc[:,i]
            plt.tick_params(axis='x',length=0)
            plt.ylabel(xlabel,fontdict={'size': 16})
            plt.xlabel(ylabel,fontdict={'size': 16})
            plt.yticks(size = 12)
            plt.xticks(size = 12)
            plt.ylim(bottom=0)
            plt.xlim(left=-1)
            plt.yticks(np.arange(0,1.2,0.2),[f'{i}%' for i in range(0,120,20)])
            plt.legend(frameon=False,fontsize=16,bbox_to_anchor=(1,1))
            plt.tight_layout()
            ax=plt.gca()
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            plt.savefig(savename+'.pdf')
        elif ~percentage:
        # plotting normal stacking barplot
            if direction=='h':
                current_palette=sns.color_palette("coolwarm",num)
                params = {
                    'figure.figsize': figsize
                }
                plt.rcParams.update(params)
                left = len(data) * [0]
                for i in range(0,num):
                    plt.barh(data.index,data.iloc[:,i],left=left,
                             height=barwidth,
                             label=data.columns[i],
                             color=current_palette[i],
                             edgecolor=None,)
                    left=left+data.iloc[:,i]
                plt.tick_params(axis='x',length=0)
                plt.xlabel(xlabel,fontdict={'size': 16})
                plt.ylabel(ylabel,fontdict={'size': 16})
                plt.yticks(size = 12)
                plt.xticks(size = 12)
                plt.ylim(bottom=-1)
                plt.xlim(left=0)
                #plt.xticks(np.arange(0,1.2,0.2),[f'{i}%' for i in range(0,120,20)])
                plt.legend(frameon=False,fontsize=16)
                plt.tight_layout()
                ax=plt.gca()
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                plt.savefig(savename+'.pdf')
            elif direction=='v':
                # plotting vertical barplot
                current_palette=sns.color_palette("coolwarm",num)
                params = {
                    'figure.figsize': figsize
                }
                plt.rcParams.update(params)
                bottom = len(data) * [0]
                for i in range(0,num):
                    plt.bar(data.index,data.iloc[:,i],bottom=bottom,
                            width=barwidth,
                            label=data.columns[i],
                            color=current_palette[i],
                            edgecolor=None,)
                    bottom=bottom+data.iloc[:,i]
                plt.tick_params(axis='x',length=0)
                plt.ylabel(xlabel,fontdict={'size': 16})
                plt.xlabel(ylabel,fontdict={'size': 16})
                plt.yticks(size = 12)
                plt.xticks(size = 12)
                plt.ylim(bottom=0)
                plt.xlim(left=-1)
                #plt.yticks(np.arange(0,1.2,0.2),[f'{i}%' for i in range(0,120,20)])
                plt.legend(frameon=False,fontsize=16)
                plt.tight_layout()
                ax=plt.gca()
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                plt.savefig(savename+'.pdf')

#---------------------------------------------------End of barplot-----------------------------------------------------------------------#
