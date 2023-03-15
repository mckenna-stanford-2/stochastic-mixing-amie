#==========================================================================
# Title: plot_precip_efficiency_fig9.py
# Author: McKenna W. Stanford
# Utility: Plots vertically integrated process rates.
# Used to make Fig. 9 of the manuscript.
#==========================================================================

#------------------------------------------
# Imports
#------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import glob
import xarray
import seaborn as sns
import matplotlib
from matplotlib.cm import get_cmap
from netCDF4 import Dataset
import wrf
import pickle
from matplotlib.lines import Line2D
import matplotlib.dates as mdates
#------------------------------------------
# Parameters
#------------------------------------------
dfmt = mdates.DateFormatter('%d-%H')
plt.rc('text',usetex=True)


path = '/glade/scratch/mckenna/AMIE/amie_post/proc_rates/'
files = sorted(glob.glob(path+'*_proc_rates.p'))
print(len(files))



# Optional to skip if already calcualted
icalc = 0
if icalc == 1.:
    
    sim_names = ['baseline',\
                 'stoch1_short',\
                 'stoch2_short',\
                 'stoch3_short',\
                 'stoch4_short',\
                 'stoch5_short',\
                 'stoch1_long',\
                 'stoch2_long',\
                 'stoch3_long',\
                 'stoch4_long',\
                 'stoch5_long',\
                 'theta_pert1',\
                 'theta_pert2',\
                 'theta_pert3',\
                 'theta_pert4',\
                 'theta_pert5',\
                 'no_mixing',\
                 '4x']
    
    pre_time_mean_dict = {}
    pre_domain_mean_dict = {}

    cond_time_mean_dict = {}
    cond_domain_mean_dict = {}

    evap_time_mean_dict = {}
    evap_domain_mean_dict = {}

    #pe_time_mean_dict = {}
    #pe_domain_mean_dict = {}
    
    for sim_ii in range(len(sim_names)):
        print('Simulation:',sim_names[sim_ii])
        file = path+sim_names[sim_ii]+'_proc_rates.p'
        pkl_file = open(file,'rb')
        proc_dict = pickle.load(pkl_file)
        pkl_file.close()

        col_int_cond = proc_dict['col_int_cond']
        col_int_evap = proc_dict['col_int_evap']
        rainnc = proc_dict['rainnc']
        snownc = proc_dict['snownc']
        time = proc_dict['time']

        pre = rainnc + snownc
        # Convert units
        pre = pre*4./3600.*997./1000. # convert to mass flux
        col_int_cond = col_int_cond*4./3600. # convert to per second
        col_int_evap = col_int_evap*4./3600. # convert to per second

        # Calculate precipitation efficient
        #pe = np.zeros(np.shape(pre))
        #dumid = np.where(col_int_cond > 0.)
        #pe[dumid] = pre[dumid]/col_int_cond[dumid]

        # Calculate time mean
        col_int_cond_time_mean = np.nanmean(col_int_cond,axis=(1,2))
        col_int_evap_time_mean = np.nanmean(col_int_evap,axis=(1,2))
        pre_time_mean = np.nanmean(pre,axis=(1,2))
        #pe_time_mean = np.nanmean(pe,axis=(1,2))

        # Calculate domain mean
        col_int_cond_domain_mean = np.nanmean(col_int_cond)
        col_int_evap_domain_mean = np.nanmean(col_int_evap)
        pre_domain_mean = np.nanmean(pre)
        #pe_domain_mean = np.nanmean(pe)
        
        dumstr = sim_names[sim_ii]

        cond_time_mean_dict[dumstr] = col_int_cond_time_mean
        evap_time_mean_dict[dumstr] = col_int_evap_time_mean
        pre_time_mean_dict[dumstr] = pre_time_mean
        #pe_time_mean_dict[dumstr] = pe_time_mean

        cond_domain_mean_dict[dumstr] = col_int_cond_domain_mean
        evap_domain_mean_dict[dumstr] = col_int_evap_domain_mean
        pre_domain_mean_dict[dumstr] = pre_domain_mean
        #pe_domain_mean_dict[dumstr] = pe_domain_mean


        iplot = False
        if iplot:
            fig = plt.figure(figsize=(12,12))
            Fontsize=14
            ax1 = fig.add_subplot(411)
            ax2 = fig.add_subplot(412)
            ax3 = fig.add_subplot(413)
            ax4 = fig.add_subplot(414)
            axlist = [ax1,ax2,ax3,ax4]
            for ax in axlist:
                ax.grid(ls='dotted',lw=1,c='dimgrey')
                ax.tick_params(labelsize=Fontsize)
                ax.set_xlabel('UTC Time',fontsize=Fontsize)

            ax1.plot(time,pre_time_mean,lw=2,c='k',ls='dashed')
            ax2.plot(time,col_int_cond_time_mean,lw=2,c='k',ls='dashed')
            ax3.plot(time,col_int_evap_time_mean,lw=2,c='k',ls='dashed')
            ax4.plot(time,pre_time_mean/col_int_cond_time_mean,lw=2,c='k',ls='dashed')           

            ax1.set_ylabel('PRE',fontsize=Fontsize)
            ax2.set_ylabel('PRE',fontsize=Fontsize)
            ax3.set_ylabel('EVAP',fontsize=Fontsize)
            ax4.set_ylabel('PE',fontsize=Fontsize)

            plt.subplots_adjust(hspace=0.3)
            plt.show()
            plt.close()


    out_dict = {'cond_time_mean':cond_time_mean_dict,\
                'evap_time_mean':evap_time_mean_dict,\
                'pre_time_mean':pre_time_mean_dict,\
                #'pe_time_mean':pe_time_mean_dict,\
                'cond_domain_mean':cond_domain_mean_dict,\
                'evap_domain_mean':evap_domain_mean_dict,\
                'pre_domain_mean':pre_domain_mean_dict,\
                #'pe_domain_mean':pe_domain_mean_dict,\
                'time':time,\
               }

    savdir = '/glade/scratch/mckenna/AMIE/amie_post/proc_rates/'
    f = open(savdir+'mean_proc_rates_all_sims.p','wb')
    pickle.dump(out_dict,f)
    f.close()

#==========================
# Read in file
#==========================
path = '/glade/scratch/mckenna/AMIE/amie_post/proc_rates/'

pkl_file = open(path+'mean_proc_rates_all_sims.p','rb')
proc_dict = pickle.load(pkl_file)
pkl_file.close()    
    
    
    
    
    
#==============================================
# Fig 9 of manuscript - time series of domain-mean
# vertically integrated process rates
#==============================================
fig = plt.figure(figsize=(8.3,16))
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)
axlist = [ax1,ax2,ax3,ax4]
time = proc_dict['time']
Fontsize=22
dum = 0
labs = ['(a)','(b)','(c)','(d)']
for ax in axlist:
    ax.grid(color='dimgrey',ls='dotted',lw=1)
    ax.set_xlabel('UTC Time [Day-Hour]',fontsize=Fontsize)
    ax.tick_params(labelsize=Fontsize)
    ax.xaxis.set_major_formatter(dfmt)
    ax.set_xticks(time[::24])
    ax.set_xlim(time[0],time[-1])
    ax.text(0.01,0.77,labs[dum],transform=ax.transAxes,fontsize=Fontsize*2.2)
    dum+=1

dumlw=3
ax1.set_ylabel('Domain Mean PRE\n[kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize)
ax2.set_ylabel('Domain Mean $\\langle$COND$\\rangle$\n[kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize)
ax3.set_ylabel('Domain Mean $\\langle$EVAP$\\rangle$\n[kg m$^{-2}$ s$^{-1}$]',fontsize=Fontsize)
ax4.set_ylabel('Domain Mean PE',fontsize=Fontsize)

# Pre
ax1.plot(proc_dict['time'],proc_dict['pre_time_mean']['baseline'],lw=dumlw,c='k',ls='solid')
ax1.plot(proc_dict['time'],proc_dict['pre_time_mean']['4x'],lw=dumlw,c='k',ls='dashed')
ax1.plot(proc_dict['time'],proc_dict['pre_time_mean']['no_mixing'],lw=dumlw,c='k',ls='dotted')

# Cond
ax2.plot(proc_dict['time'],proc_dict['cond_time_mean']['baseline'],lw=dumlw,c='k',ls='solid')
ax2.plot(proc_dict['time'],proc_dict['cond_time_mean']['4x'],lw=dumlw,c='k',ls='dashed')
ax2.plot(proc_dict['time'],proc_dict['cond_time_mean']['no_mixing'],lw=dumlw,c='k',ls='dotted')

# Evap
ax3.plot(proc_dict['time'],proc_dict['evap_time_mean']['baseline'],lw=dumlw,c='k',ls='solid')
ax3.plot(proc_dict['time'],proc_dict['evap_time_mean']['4x'],lw=dumlw,c='k',ls='dashed')
ax3.plot(proc_dict['time'],proc_dict['evap_time_mean']['no_mixing'],lw=dumlw,c='k',ls='dotted')

# PE
ax4.plot(proc_dict['time'],\
         proc_dict['pre_time_mean']['baseline']/proc_dict['cond_time_mean']['baseline'],\
         lw=dumlw,c='k',ls='solid')
ax4.plot(proc_dict['time'],\
         proc_dict['pre_time_mean']['4x']/proc_dict['cond_time_mean']['4x'],\
         lw=dumlw,c='k',ls='dashed')
ax4.plot(proc_dict['time'],\
         proc_dict['pre_time_mean']['no_mixing']/proc_dict['cond_time_mean']['no_mixing'],\
         lw=dumlw,c='k',ls='dotted')


custom_lines = [Line2D([0],[0],color = 'black',lw=3),\
                Line2D([0],[0],color='black',lw=3,ls='dotted'),\
                Line2D([0],[0],color='black',lw=3,ls='dashed')]

ax4.legend(custom_lines,['BASELINE','NO\_MIXING','4X'],fontsize=Fontsize,\
                loc='lower center',bbox_to_anchor=(0.5,-0.65),ncol=3)

ax1.set_ylim(0,0.00031)
ax2.set_ylim(0,0.00031)
ax3.set_ylim(0,0.00031)

plt.subplots_adjust(hspace=0.35)

outfile = 'fig_09.png'
plt.savefig('/glade/u/home/mckenna/figures/amie_paper/'+outfile,dpi=200,bbox_inches='tight')
#plt.show()
plt.close()
    

    
print(aaaa)   
    
    
    
#==========================
# Calculate time mean 
#==========================


mean_cond = {}
mean_evap = {}
mean_pre = {}
mean_pe = {}

rel_diff_cond = {}
rel_diff_evap = {}
rel_diff_pre = {}
rel_diff_pe = {}    


dum = 0
for key,val in pre.items():
    key2 = key
    key2 = key2[0:-4]
    tmp_conv_strat = conv_strat_dict[key2]
    conv_id = np.where(tmp_conv_strat == 1.)

    mean_cond[key] = np.nanmean(cond[key][conv_id])
    mean_pre[key] = np.nanmean(pre[key][conv_id])
    mean_evap[key] = np.nanmean(evap[key][conv_id])
    mean_pe[key] = mean_pre[key]/mean_cond[key]

    
    
    
for key,val in pre.items():
    if key != 'baseline_new':
        rel_diff_cond[key] = (mean_cond[key] - mean_cond['baseline_new'])/(mean_cond['baseline_new'])
        rel_diff_evap[key] = (mean_evap[key] - mean_evap['baseline_new'])/(mean_evap['baseline_new'])
        rel_diff_pre[key] = (mean_pre[key] - mean_pre['baseline_new'])/(mean_pre['baseline_new'])
        rel_diff_pe[key] = (mean_pe[key] - mean_pe['baseline_new'])/(mean_pe['baseline_new'])


sim_names = ['BASELINE',\
             'STOCH1\_SHORT',\
             'STOCH2\_SHORT',\
             'STOCH3\_SHORT',\
             'STOCH4\_SHORT',\
             'STOCH5\_SHORT',\
             'STOCH1\_LONG',\
             'STOCH2\_LONG',\
             'STOCH3\_LONG',\
             'STOCH4\_LONG',\
             'STOCH5\_LONG',\
             '$\\theta$-pert1',\
             '$\\theta$-pert2',\
             '$\\theta$-pert3',\
             '$\\theta$-pert4',\
             '$\\theta$-pert5',\
             'NO\_MIXING',\
             '4X',\
            ]

key_order = ['baseline','stoch1_short','stoch2_short','stoch3_short','stoch4_short',\
            'stoch5_short','stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long',\
            'theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5',\
            'no_mixing','4x']

sim_names = sim_names[0:6]
key_order = key_order[0:6]
    
dum=0
for key in key_order:
    if key == 'baseline':
        key = key+'_new'
        val = mean_pre[key]
        print(sim_names[dum],' & ',\
              np.round(val*1.e4,2),\
            ' & ',' - ',' & ',\
              np.round(mean_cond[key]*1.e4,2),
              ' & ',' - ',' & ',\
              np.round(mean_evap[key]*1.e4,2),\
              ' & ',' - ',' & ',\
              np.round(mean_pe[key],2),\
              ' & ',' - ',\
              ' \\\\',\
             )
        dum+=1
    else:
        key = key+'_new'
        val = mean_pre[key]
        print(sim_names[dum],' & ',np.round(val*1.e4,2),\
              ' & ',np.round(rel_diff_pre[key]*100.,2),' & ',\
              np.round(mean_cond[key]*1.e4,2),\
              ' & ',np.round(rel_diff_cond[key]*100.,2),' & ',\
              np.round(mean_evap[key]*1.e4,2),\
              ' & ',np.round(rel_diff_evap[key]*100.,2),' & ',\
              np.round(mean_pe[key],2),\
              ' & ',np.round(rel_diff_pe[key]*100.,2),\
              ' \\\\',\
             )
        dum+=1
        
print(aaaa)      
        


#==========================
# Calculate time time mean 
#and then do domin averaging
#==========================


#if True:
if False:
    
    mean_cond = {}
    mean_evap = {}
    mean_pre = {}
    mean_pe = {}

    rel_diff_cond = {}
    rel_diff_evap = {}
    rel_diff_pre = {}
    rel_diff_pe = {}    
    
    
    dum = 0
    for key,val in pre.items():
        
        mean_cond[key] = np.nanmean(cond[key])
        mean_pre[key] = np.nanmean(pre[key])
        mean_evap[key] = np.nanmean(evap[key])
        mean_pe[key] = mean_pre[key]/mean_cond[key]
        
    for key,val in pre.items():
        if key != 'baseline_new':
            rel_diff_cond[key] = (mean_cond[key] - mean_cond['baseline_new'])/(mean_cond['baseline_new'])
            rel_diff_evap[key] = (mean_evap[key] - mean_evap['baseline_new'])/(mean_evap['baseline_new'])
            rel_diff_pre[key] = (mean_pre[key] - mean_pre['baseline_new'])/(mean_pre['baseline_new'])
            rel_diff_pe[key] = (mean_pe[key] - mean_pe['baseline_new'])/(mean_pe['baseline_new'])
        
    #print(aaaaaaa)
        
    
    sim_names = ['BASELINE',\
                 'STOCH1\_SHORT',\
                 'STOCH2\_SHORT',\
                 'STOCH3\_SHORT',\
                 'STOCH4\_SHORT',\
                 'STOCH5\_SHORT',\
                 'STOCH1\_LONG',\
                 'STOCH2\_LONG',\
                 'STOCH3\_LONG',\
                 'STOCH4\_LONG',\
                 'STOCH5\_LONG',\
                 '$\\theta$-pert1',\
                 '$\\theta$-pert2',\
                 '$\\theta$-pert3',\
                 '$\\theta$-pert4',\
                 '$\\theta$-pert5',\
                 'NO\_MIXING',\
                 '4X',\
                ]
    
    key_order = ['baseline','stoch1_short','stoch2_short','stoch3_short','stoch4_short',\
                'stoch5_short','stoch1_long','stoch2_long','stoch3_long','stoch4_long','stoch5_long',\
                'theta_pert1','theta_pert2','theta_pert3','theta_pert4','theta_pert5',\
                'no_mixing','4x']
    dum=0
    for key in key_order:
        if key == 'baseline':
            key = key+'_new'
            val = mean_pre[key]
            print(sim_names[dum],' & ',\
                  np.round(val*1.e4,2),\
                ' & ',' - ',' & ',\
                  np.round(mean_cond[key]*1.e4,2),
                  ' & ',' - ',' & ',\
                  np.round(mean_evap[key]*1.e4,2),\
                  ' & ',' - ',' & ',\
                  np.round(mean_pe[key],2),\
                  ' & ',' - ',\
                  ' \\\\',\
                 )
            dum+=1
        else:
            key = key+'_new'
            val = mean_pre[key]
            print(sim_names[dum],' & ',np.round(val*1.e4,2),\
                  ' & ',np.round(rel_diff_pre[key]*100.,2),' & ',\
                  np.round(mean_cond[key]*1.e4,2),\
                  ' & ',np.round(rel_diff_cond[key]*100.,2),' & ',\
                  np.round(mean_evap[key]*1.e4,2),\
                  ' & ',np.round(rel_diff_evap[key]*100.,2),' & ',\
                  np.round(mean_pe[key],2),\
                  ' & ',np.round(rel_diff_pe[key]*100.,2),\
                  ' \\\\',\
                 )
            dum+=1
        #print('\\hline')        
    # print out lines that will be used for Table 5.3
    dum = 0
    for key,val in mean_pre.items():
        if key == 'baseline_new':
            print(sim_names[dum],' & ',\
                  np.round(val*1.e4,2),\
                ' & ',' - ',' & ',\
                  np.round(mean_cond[key]*1.e4,2),
                  ' & ',' - ',' & ',\
                  np.round(mean_evap[key]*1.e4,2),\
                  ' & ',' - ',' & ',\
                  np.round(mean_pe[key],2),\
                  ' & ',' - ',\
                  ' \\\\',\
                 )
        else:
            continue
            print(sim_names[dum],' & ',np.round(val*1.e4,2),\
                  ' & ',np.round(rel_diff_pre[key]*100.,2),' & ',\
                  np.round(mean_cond[key]*1.e4,2),\
                  ' & ',np.round(rel_diff_cond[key]*100.,2),' & ',\
                  np.round(mean_evap[key]*1.e4,2),\
                  ' & ',np.round(rel_diff_evap[key]*100.,2),' & ',\
                  np.round(mean_pe[key],2),\
                  ' & ',np.round(rel_diff_pe[key]*100.,2),\
                  ' \\\\',\
                 )
        #print('\\hline')
        
        dum+=1
        
        
    # print out lines that will be used for Table 5.3
    dum = 0
    for key,val in mean_pre.items():
        if key == 'baseline_new':
            print(sim_names[dum],' & ',np.round(val*1.e4,2),
                  ' \\\\',\
                 )
        else:
            print(sim_names[dum],' & ',np.round(val*1.e4,2),\
                  ' & ',np.round(rel_diff_pre[key]*100.,2),' & ',\
                  ' \\\\',\
                 )
        #print('\\hline')
        
        dum+=1        
        
        
        
    for key,val in rel_diff_cond.items():
        print(key,val*100.)


    
    
#==============================================
# Fig with only 4x, no mixing, and baseline
#==============================================
#if True:
if False:
    tmp_pre = {'baseline':pre['baseline_new'],'no_mixing':pre['no_mixing_new'],'4x':pre['4x_new']}
    tmp_cond = {'baseline':cond['baseline_new'],'no_mixing':cond['no_mixing_new'],'4x':cond['4x_new']}
    tmp_evap = {'baseline':evap['baseline_new'],'no_mixing':evap['no_mixing_new'],'4x':evap['4x_new']}


    fig = plt.figure(figsize=(6,16))
    ax1 = fig.add_subplot(411)
    ax2 = fig.add_subplot(412)
    ax3 = fig.add_subplot(413)
    ax4 = fig.add_subplot(414)
    #ax5 = fig.add_subplot(515)
    import matplotlib.dates as mdates

    dfmt = mdates.DateFormatter('%d-%H')

    colors = ['black','red','blue']
    lss = ['solid','dotted','dashed']
    dum = 0
    labs = ['(a)','(b)','(c)','(d)']
    axlist = [ax1,ax2,ax3,ax4]
    for ax in axlist:
        ax.grid(which='both',color='grey',ls='solid',lw=1)
        ax.set_xlabel('Time [UTC]',fontsize=20)
        ax.tick_params(labelsize=20)
        ax.xaxis.set_major_formatter(dfmt)
        ax.set_xticks([time[0],time[24],time[48],time[72],time[96]])
        ax.set_xlim(time[0],time[96])
        ax.text(0.02,0.8,labs[dum],transform=ax.transAxes,fontsize=40)
        dum+=1

    ax2.set_ylabel('Domain Mean\n$\\langle$COND$\\rangle$\n[kg m$^{-2}$ s$^{-1}$]',fontsize=20)
    ax3.set_ylabel('Domain Mean\n$\\langle$EVAP$\\rangle$\n[kg m$^{-2}$ s$^{-1}$]',fontsize=20)
    ax1.set_ylabel('Domain Mean \nPRE\n[kg m$^{-2}$ s$^{-1}$]',fontsize=20)
    #ax4.set_ylabel('Domain Mean \nRES\n[kg m$^{-2}$ s$^{-1}$]',fontsize=15)
    #ax4.set_ylabel('Domain Mean \nEVAP/COND',fontsize=15)
    ax4.set_ylabel('Domain Mean \nPE',fontsize=20)


    dum = 0
    for sim in tmp_cond.keys():
        print(sim)
        tmpcond = np.array(tmp_cond[sim])
        tmpevap = np.array(tmp_evap[sim])
        tmppre = np.array(tmp_pre[sim])
        pe = np.nanmean(tmppre,axis=(1,2))/np.nanmean(tmpcond,axis=(1,2))
        res = tmpcond - tmpevap - tmppre
        #evap_cond = tmpevap/tmp
        #print(aaa)

        #res = np.nanmean(tmpcond,axis=(1,2))-np.nanmean(tmpevap,axis=(1,2))-np.nanmean(tmppre,axis=(1,2))
        #res = np.nanmean(tmpcond,axis=(1,2))-np.nanmean(tmpevap,axis=(1,2))-np.nanmean(tmppre,axis=(1,2))



        ax1.plot(time,np.nanmean(tmppre,axis=(1,2)),color='k',ls=lss[dum],lw=3)
        ax2.plot(time,np.nanmean(tmpcond,axis=(1,2)),color='k',ls=lss[dum],lw=3)
        ax3.plot(time,np.nanmean(tmpevap,axis=(1,2)),color='k',ls=lss[dum],lw=3)
        #ax4.plot(time,np.nanmean(res,axis=(1,2)),color='k',ls=lss[dum],lw=3)
        #ax4.plot(time,res,color='k',ls=lss[dum],lw=3)
        #ax4.plot(time,np.nanmean(tmpevap,axis=(1,2))/np.nanmean(tmpcond,axis=(1,2)),color='k',ls=lss[dum],lw=3)
        ax4.plot(time,pe,color='k',ls=lss[dum],lw=3)
        dum = dum + 1

    custom_lines = [Line2D([0],[0],color = 'black',lw=4),\
                    Line2D([0],[0],color='black',lw=4,ls='dotted'),\
                    Line2D([0],[0],color='black',ls='dashed',lw=4)]

    ax4.legend(custom_lines,['BASELINE','NO\_MIXING','4X'],fontsize=20,\
                    loc='lower center',bbox_to_anchor=(0.5,-0.8),ncol=2)

    ax1.set_ylim(0,0.00031)
    ax2.set_ylim(0,0.00031)
    ax3.set_ylim(0,0.00031)
    #ax4.set_ylim(0,0.00031)
    
    plt.subplots_adjust(hspace=0.4)

    plt.show()

