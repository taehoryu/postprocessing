#!/usr/bin/env python
# -*- coding: utf-8 -*-

import yt
from yt.mods import *
import glob
import matplotlib
import os
#matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
np.set_printoptions(threshold=np.inf)

modelname=['radiation',\
          r'$     \kappa (\rho>10^{-14}{\mathrm{g}} {\mathrm{cm}^{-3}}) = \kappa$',\
          r'$     \kappa_(\rho<10^{-14}{\mathrm{g}} {\mathrm{cm}^{-3}}) = 0$',\
          r'$     F=1F_{\mathrm{irr}}$']
color=['black','black','black','black']
bbox_props = [dict(facecolor='none', edgecolor=color[0], pad=0.8),\
              dict(facecolor='none', edgecolor=color[1], pad=0.8),\
              dict(facecolor='none', edgecolor=color[2], pad=0.8),\
              dict(facecolor='none', edgecolor=color[3], pad=0.8)]
#box=1, no box=0
textbox=[0,0,0,0]
#Left top corner
textposition_LT=[0.15,0.83]
#Right top corner
textposition_RT=[0.6,0.83]
#Left bottom corner
textposition_LB=[0.15,0.3]
#Right bottom corner
textposition_RB=[0.6,0.3]

def index_exists(ls, i):
    return (0 <= i < len(ls)) or (-len(ls) <= i < 0)


numberoflines=5.0
Dimension=2
number_slice_x=100
max_time=10010.0
x_limit=[0.0,1.0e9]
P_range=[0.001,2.e2]
T_range=[2300,3000]
l=0
i=0
j=0
jj=0

every_nthfile=3
#print 'Listing all profile/dat files'
profilefilelist = sorted(glob.glob('plt*'), key=lambda name: int(name[3:]))
#raw_input('Press ENTER to continue...')
number_plt=0
time=0
radiation_field=1
for i in profilefilelist:
  number_plt=number_plt+1
  ds=load(i)
  time=np.maximum(time,float(ds.current_time))
  if number_plt==1:
    ds.index
    for j in sorted(ds.field_list):
      if j[1]=="rad":
        radiation_field=1
        print "rad"

print "radiation field", radiation_field
time_difference_x_den=time/numberoflines
time_difference_x_P=time/numberoflines
time_difference_x_T=time/numberoflines
time_difference_x_entropy=time/numberoflines
time_difference_x_rad=time/numberoflines
time_difference_T_P=time/numberoflines
time_difference_x_opacity=time/numberoflines
time_difference_x_Mach=time/numberoflines
time_difference_v_P=time/numberoflines

#print profilefilelist
print number_plt,i
average_T=0
average_P=0
average_P_in_x=0
average_T_in_x=0
average_Mach_in_x=0
T_01=0
T_10=0
T_1=0
KE_1=0
KE_10=0
HEAT_01=0
HEAT_1=0
HEAT_10=0
KE_1_list=0
HEAT_01_list=0
HEAT_1_list=0
HEAT_10_list=0
T_01_list=0
T_1_list=0
T_10_list=0
time=0
i=0
j=0
jj=0








if radiation_field==1:
  average_T_in_x=0
  average_P_in_x=0
  average_T=0
  average_P=0
  j=0
  jj=0
  for i in profilefilelist:
  #  if np.mod(j,every_nthfile)==0:
    ds=load(i)
    time=float(ds.current_time)
    filename_time=str(int(time))
    if time>=float(jj)*time_difference_x_rad and time<max_time:
      jj=jj+1
      print('Adding lines at t='+filename_time+'s for U-height plot')

      if Dimension==1:
        my_ray = ds.ortho_ray(0,(0,0))
        srt=np.argsort(my_ray['x'])
        x_coord=np.array(my_ray['x'][srt])/1.e5
      elif Dimension==2:
        my_ray = ds.ortho_ray(1,(0,0))
        srt=np.argsort(my_ray['y'])
        x_coord=np.array(my_ray['y'][srt])/1.e5
      radiation=np.array(my_ray['rad'][srt])
      plt.semilogy(x_coord,radiation,label='t='+filename_time+'s')
      plt.legend()
      plt.ylabel(r'$U [{\rm erg} {\rm cm}^{-3}]$')
      plt.xlabel(r'$y\/[\mathrm{km}]$')
    #plt.legend(loc=3)
      plt.yscale('log')
      text_num=0
      while (text_num<=3):
        if textbox[text_num]==0:
          plt.figtext(textposition_RT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
        elif textbox[text_num]==1:
          plt.figtext(textposition_RT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
        text_num=text_num+1
      plt.title('Radiation energy density $U$ vs height')
      plt.savefig("figure_x_radiation.png")
  
  
    j=j+1

  print("radiation plot made")
  plt.close()









if radiation_field==1:
  average_T_in_x=0
  average_P_in_x=0
  average_T=0
  average_P=0
  j=0
  jj=0
  for i in profilefilelist:
  #  if np.mod(j,every_nthfile)==0:
    ds=load(i)
    time=float(ds.current_time)
    filename_time=str(int(time))
    if time>=float(jj)*time_difference_x_rad and time<max_time:
      jj=jj+1
      print('Adding lines at t='+filename_time+'s for U-height plot')

      if Dimension==1:
        my_ray = ds.ortho_ray(0,(0,0))
        srt=np.argsort(my_ray['x'])
        x_coord=np.array(my_ray['x'][srt])/1.e5
      elif Dimension==2:
        my_ray = ds.ortho_ray(1,(0,0))
        srt=np.argsort(my_ray['y'])
        x_coord=np.array(my_ray['y'][srt])/1.e5
      radiation=np.array(my_ray['rad'][srt])
      plt.semilogy(x_coord,radiation,label='t='+filename_time+'s')
      plt.legend(loc=1)
      plt.ylabel(r'$U [{\rm erg} {\rm cm}^{-3}]$')
      plt.xlabel(r'$y\/[\mathrm{km}]$')
    #plt.legend(loc=3)
      plt.ylim(0.0,1)
      plt.xlim(0.5e5,0.7e5)
      plt.yscale('linear')
      text_num=0
      while (text_num<=3):
        if textbox[text_num]==0:
          plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
        elif textbox[text_num]==1:
          plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
        text_num=text_num+1
      plt.title('Radiation energy density $U$ vs height')
      plt.savefig("figure_x_radiation_zoom.png")
  
  
    j=j+1

  print("radiation plot made")
  plt.close()



average_T_in_x=0
average_P_in_x=0
average_T=0
average_P=0
old_temp=0
j=0
jj=0

for i in profilefilelist:
  #  if np.mod(j,every_nthfile)==0:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  
  
  
    #
    
    
  kappa=0.18*1e-9
  SB_constant=5.6704*1e-5
    
    
  if Dimension==1:
    my_ray = ds.ortho_ray(0,(0,0))
    srt=np.argsort(my_ray['x'])
    x_coord=np.array(my_ray['x'][srt])/1.e5
    average_Frlabx_in_x=np.array(my_ray['Frlabx'][srt])
    average_Frcomx_in_x=np.array(my_ray['Frcomx'][srt])/1.e6
  
  elif Dimension==2:
    my_ray = ds.ortho_ray(1,(0,0))
    srt=np.argsort(my_ray['y'])
    x_coord=np.array(my_ray['y'][srt])/1.e5

  if Dimension==2 :

    while (l<number_slice_x):
      position_x=x_limit[0]+float(l)*(x_limit[1]-x_limit[0])/float(number_slice_x)
        #     print "here",position_x,l,time
      my_ray = ds.ortho_ray(1,(0,position_x))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
      temp=np.array(my_ray['Temp'][srt])
      press=np.array(my_ray['pressure'][srt])/1.e6
      
      if l==0:
        max_T_in_x=temp
        min_T_in_x=temp
      else:
        max_T_in_x=np.maximum(max_T_in_x,temp)
        min_T_in_x=np.minimum(min_T_in_x,temp)

      average_T_in_x=(average_T_in_x*float(l)+temp)/(l+1)
      average_P_in_x=(average_P_in_x*float(l)+press)/(l+1)
      l=l+1
#  average_T_in_x=average_T_in_x/float(l-1)
#  average_P_in_x=average_P_in_x/float(l-1)

  l=0
  if j==0 :
    Frlabx0=np.array(my_ray['Frlabx'][srt])
    Frcomx0=np.array(my_ray['Frcomx'][srt])


    #plt.clf()
    #if j<=number_plt:
  if (time>=float(jj)*time_difference_T_P and time<max_time) or (j+1==number_plt):
    jj=jj+1
    plt.plot(x_coord,average_Frlabx_in_x,label='t='+filename_time+'s')
    if Dimension==2:
      plt.fill_between(average_P_in_x,max_T_in_x,min_T_in_x,alpha=0.5)
    plt.legend(loc=3)
    print('Adding lines at t='+filename_time+'s for T-P plot')
    plt.xlabel(r'height$[\mathrm{km}$]')
    plt.ylabel(r'lab flux$[{\mathrm{ergcm^{-2}s^{-1}}}]$')
    #plt.xlim(P_range[0],P_range[1])
    #plt.ylim(T_range[0],T_range[1])
    #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    plt.xscale('log')
    plt.yscale('log')
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
    plt.title('$x$ vs $F_{\mathrm{lab}}$')
    plt.savefig("figure_x_Frlabx.png")
  
  j=j+1
print("Frlabx plot average made")
plt.close()


average_T_in_x=0
average_P_in_x=0
average_T=0
average_P=0
old_temp=0
j=0
jj=0

for i in profilefilelist:
  #  if np.mod(j,every_nthfile)==0:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  
  
  
    #
    
    
  kappa=0.18*1e-9
  SB_constant=5.6704*1e-5
    
    
  if Dimension==1:
    my_ray = ds.ortho_ray(0,(0,0))
    srt=np.argsort(my_ray['x'])
    x_coord=np.array(my_ray['x'][srt])/1.e5
    average_Frlabx_in_x=np.array(my_ray['Frlabx'][srt])
    average_Frcomx_in_x=np.array(my_ray['Frcomx'][srt])/1.e6
  
  elif Dimension==2:
    my_ray = ds.ortho_ray(1,(0,0))
    srt=np.argsort(my_ray['y'])
    x_coord=np.array(my_ray['y'][srt])/1.e5

  if Dimension==2 :

    while (l<number_slice_x):
      position_x=x_limit[0]+float(l)*(x_limit[1]-x_limit[0])/float(number_slice_x)
        #     print "here",position_x,l,time
      my_ray = ds.ortho_ray(1,(0,position_x))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
      temp=np.array(my_ray['Temp'][srt])
      press=np.array(my_ray['pressure'][srt])/1.e6
      
      if l==0:
        max_T_in_x=temp
        min_T_in_x=temp
      else:
        max_T_in_x=np.maximum(max_T_in_x,temp)
        min_T_in_x=np.minimum(min_T_in_x,temp)

      average_T_in_x=(average_T_in_x*float(l)+temp)/(l+1)
      average_P_in_x=(average_P_in_x*float(l)+press)/(l+1)
      l=l+1
#  average_T_in_x=average_T_in_x/float(l-1)
#  average_P_in_x=average_P_in_x/float(l-1)

  l=0
  if j==0 :
    Frlabx0=np.array(my_ray['Frlabx'][srt])
    Frcomx0=np.array(my_ray['Frcomx'][srt])


    #plt.clf()
    #if j<=number_plt:
  if (time>=float(jj)*time_difference_T_P and time<max_time) or (j+1==number_plt):
    jj=jj+1
    plt.plot(x_coord,average_Frlabx_in_x,label='t='+filename_time+'s')
    if Dimension==2:
      plt.fill_between(average_P_in_x,max_T_in_x,min_T_in_x,alpha=0.5)
    plt.legend(loc=3)
    print('Adding lines at t='+filename_time+'s for T-P plot')
    plt.xlabel(r'height$[\mathrm{km}$]')
    plt.ylabel(r'lab flux$[{\mathrm{ergcm^{-2}s^{-1}}}]$')
    plt.xlim(0.5e5,0.8e5)
    #plt.ylim(T_range[0],T_range[1])
    #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    #plt.xscale('log')
    plt.yscale('log')
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
    plt.title('$x$ vs $F_{\mathrm{lab}}$')
    plt.savefig("figure_x_Frlabx_zoom.png")
  
  j=j+1
print("Frlabx plot average made")
plt.close()


