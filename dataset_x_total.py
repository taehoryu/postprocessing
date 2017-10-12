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
          r'$     \kappa_{R} = \kappa_{0}$',\
          r'$     \kappa_{P} = \kappa_{0}$',\
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
max_time=100000.0
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





for i in profilefilelist:
  # if np.mod(j,every_nthfile)==0:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  if (time>=float(jj)*time_difference_x_den and time<max_time) or j+1==number_plt :
    jj=jj+1
    print('Adding lines at t='+filename_time+'s for density-height plot')
    location1=ds.find_max('Temp')
    location2=ds.find_min('Temp')
  
    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5
  
    ad=ds.all_data()
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    
    
    
    
    
    
    
    
    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])
  
    my_ray = ds.ortho_ray(1,(0,0))
    plt.semilogy(x_coord,dens,label='t='+filename_time+'s')
    plt.legend()
    plt.ylabel(r'$\rho\/[\mathrm{g}/\mathrm{cm}^{3}]$')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
  #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    plt.xscale('log')
  #  plt.yscale('linear')
    plt.title('Density vs height')
    plt.savefig("figure_x_den_full.png")
  j=j+1

print("density plot made")
plt.close()


for i in profilefilelist:
  # if np.mod(j,every_nthfile)==0:
  ds=load(i)
  time=float(ds.current_time)
  filename_time=str(int(time))
  if (time>=float(jj)*time_difference_x_den and time<max_time) or j+1==number_plt :
    jj=jj+1
    print('Adding lines at t='+filename_time+'s for density-height plot')
    location1=ds.find_max('Temp')
    location2=ds.find_min('Temp')
  
    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5
  
    ad=ds.all_data()
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    
    
    
    
    
    
    
    
    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])
  
    my_ray = ds.ortho_ray(1,(0,0))
    plt.semilogy(x_coord,dens,label='t='+filename_time+'s')
    plt.legend()
    plt.ylabel(r'$\rho\/[\mathrm{g}/\mathrm{cm}^{3}]$')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
  #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    plt.xscale('log')
    plt.xlim(0.0,1e4)
    plt.ylim(1e-8,1e-2)
  #  plt.yscale('linear')
    plt.title('Density vs height')
    plt.savefig("figure_x_den.png")
  j=j+1

print("density plot made")
plt.close()



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
  if time>=float(jj)*time_difference_x_P and time<max_time:
    jj=jj+1
    print('Adding lines at t='+filename_time+'s for P-height plot')
    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5
  
  
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])
    press=np.array(my_ray['pressure'][srt])/1.e6


    if j==0 :
      temp0=np.array(my_ray['Temp'][srt])
      ent0=np.array(my_ray['entropy'][srt])
      pressure0=np.array(my_ray['pressure'][srt])
      density0=np.array(my_ray['density'][srt])
      location1=np.amin(pressure0)/1e10
      location2=np.amax(pressure0)

#plt.clf()
    plt.plot(x_coord,press,label='t='+filename_time+'s')
    plt.legend()

    plt.ylabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
  #plt.ylim(location1,location2)
  #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
  #  plt.xscale('linear')
    plt.title('P vs height')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.0,1e4)
    plt.ylim(1e-4,1e3)
    plt.savefig("figure_x_P.png")
  
  
  j=j+1


print("x P plot made")
plt.close()




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
  if time>=float(jj)*time_difference_x_P and time<max_time:
    jj=jj+1
    print('Adding lines at t='+filename_time+'s for P-height plot')
    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5
  
  
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])
    ent=np.array(my_ray['entropy'][srt])
    dens=np.array(my_ray['density'][srt])
    press=np.array(my_ray['pressure'][srt])/1.e6


    if j==0 :
      temp0=np.array(my_ray['Temp'][srt])
      ent0=np.array(my_ray['entropy'][srt])
      pressure0=np.array(my_ray['pressure'][srt])
      density0=np.array(my_ray['density'][srt])
      location1=np.amin(pressure0)/1e10
      location2=np.amax(pressure0)

#plt.clf()
    plt.plot(x_coord,press,label='t='+filename_time+'s')
    plt.legend()

    plt.ylabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
  #plt.ylim(location1,location2)
  #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
  #  plt.xscale('linear')
    plt.title('P vs height')
    plt.xscale('log')		
    plt.yscale('log')
    #plt.xlim(0.0,1e4)
    #plt.ylim(1e-4,1e3)
    plt.savefig("figure_x_P_full.png")
  
  
  j=j+1


print("x P plot made")
plt.close()







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
  if time>=float(jj)*time_difference_x_T and time<max_time:
    jj=jj+1
    print('Adding lines at t='+filename_time+'s for T-height plot')
    location1=ds.find_max('Temp')
    location2=ds.find_min('Temp')
  
    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5
  
  
    ad=ds.all_data()
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])

    plt.semilogy(x_coord,temp,label='t='+filename_time+'s')
    plt.legend()
    plt.ylabel(r'$T$ [K]')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
    plt.title('T vs height')
    plt.xlim(0.0,1e4)
    plt.ylim(2300,3e3)

  #plt.legend(loc=3)
  #plt.yscale('linear')
    plt.savefig("figure_x_T.png")
  
  
  j=j+1

print("T plot made")
plt.close()


average_T_in_x=0
average_vy_in_x=0
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
    average_vy_in_x=np.array(my_ray['y_velocity'][srt])/1.e5
    average_P_in_x=np.array(my_ray['pressure'][srt])/1.e6

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
      press=np.array(my_ray['pressure'][srt])/1.e6
      y_vel=np.array(my_ray['y_velocity'][srt])/1.e5
      
      if l==0:
        max_vy_in_x=y_vel
        min_vy_in_x=y_vel
      else:
        max_vy_in_x=np.maximum(max_vy_in_x,y_vel)
        min_vy_in_x=np.minimum(min_vy_in_x,y_vel)

      average_vy_in_x=(average_vy_in_x*float(l)+y_vel)/(l+1)
      average_P_in_x=(average_P_in_x*float(l)+press)/(l+1)
      l=l+1
#  average_T_in_x=average_T_in_x/float(l-1)
#  average_P_in_x=average_P_in_x/float(l-1)

  l=0
  if j==0 :
    pressure0=np.array(my_ray['pressure'][srt])
    y_vel0=np.array(my_ray['y_velocity'][srt])/1.e5

    #plt.clf()
    #if j<=number_plt:
  if (time>=float(jj)*time_difference_v_P and time<max_time) or (j+1==number_plt and time<=max_time):
    jj=jj+1
    plt.plot(average_P_in_x,average_vy_in_x,label='t='+filename_time+'s')
    if Dimension==2:
      plt.fill_between(average_P_in_x,max_vy_in_x,min_vy_in_x,alpha=0.5)
    plt.legend(loc=1)
    print('Adding lines at t='+filename_time+'s for T-P plot')
    plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.ylabel(r'$v\/[\mathrm{km s^{-1}}]$')
    plt.xlim(P_range[0],P_range[1])
    plt.ylim(-10,10)
    #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    plt.xscale('log')
    #plt.yscale('log')
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
    plt.title(r'$v_{y}$ vs $P$')
    plt.savefig("figure_vy_P.png")
  
  j=j+1
print("v_y-P plot average made")
plt.close()




if Dimension==2:



  average_T_in_x=0
  average_vy_in_x=0
  average_vx_in_x=0
  average_P_in_x=0
  average_T=0
  average_P=0
  old_temp=0
  j=0
  jj=0
  old_time=-time_difference_v_P
  for i in profilefilelist:
  #  if np.mod(j,every_nthfile)==0:
    ds=load(i)
    time=float(ds.current_time)
    filename_time=str(int(time))
  
  
  
    #
    
    

    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
      average_vx_in_x=np.array(my_ray['x_velocity'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5

    if Dimension==2 :

      while (l<number_slice_x):
        position_x=x_limit[0]+float(l)*(x_limit[1]-x_limit[0])/float(number_slice_x)
        my_ray = ds.ortho_ray(1,(0,position_x))
        srt=np.argsort(my_ray['y'])
        x_coord=np.array(my_ray['y'][srt])/1.e5
        press=np.array(my_ray['pressure'][srt])/1.e6
        x_vel=np.array(my_ray['x_velocity'][srt])/1.e5
			
        if l==0:
          max_vx_in_x=x_vel
          min_vx_in_x=x_vel
        else:
          max_vx_in_x=np.maximum(max_vx_in_x,x_vel)
          min_vx_in_x=np.minimum(min_vx_in_x,x_vel)

        average_vx_in_x=(average_vx_in_x*float(l)+x_vel)/(l+1)
        average_P_in_x=(average_P_in_x*float(l)+press)/(l+1)
        l=l+1

    l=0
    if j==0 :
      pressure0=np.array(my_ray['pressure'][srt])
      x_vel0=np.array(my_ray['x_velocity'][srt])/1.e5

    if (time>=float(jj)*time_difference_v_P and time<max_time) or (j+1==number_plt and time<=max_time):
      jj=jj+1
      plt.plot(average_P_in_x,average_vx_in_x,label='t='+filename_time+'s')
      plt.fill_between(average_P_in_x,max_vx_in_x,min_vx_in_x,alpha=0.5)
      plt.legend(loc=1)
      print('Adding lines at t='+filename_time+'s for T-P plot')
      plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
      plt.ylabel(r'$v\/[\mathrm{km s^{-1}}]$')
      plt.xlim(P_range[0],P_range[1])
      plt.ylim(-10,10)
      plt.xscale('log')
      #plt.yscale('log')
      text_num=0
      while (text_num<=3):
        if textbox[text_num]==0:
          plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
        elif textbox[text_num]==1:
          plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
        text_num=text_num+1
      plt.title(r'$v_{x}$ vs $P$')
      plt.savefig("figure_vx_P.png")
  
    j=j+1
  print("v_x-P plot average made")
  plt.close()








average_T_in_x=0
average_P_in_x=0
average_density_in_x=0
average_T=0
average_P=0
old_temp=0
j=0
jj=0
jjj=0
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
    x_coord=np.array(my_ray['x'][srt])
    average_T_in_x=np.array(my_ray['Temp'][srt])
    average_density_in_x=np.array(my_ray['density'][srt])
    average_column_mass_in_x=np.zeros(len(average_density_in_x))
    jjj=len(average_density_in_x)
    print jjj
    while(jjj>0):
			if jjj==len(average_density_in_x):
				average_column_mass_in_x[jjj-1]=average_density_in_x[jjj-1]*x_coord[jjj-1]
			else:
				average_column_mass_in_x[jjj-1]=average_column_mass_in_x[jjj]+average_density_in_x[jjj-1]*abs(x_coord[jjj-1]-x_coord[jjj])
			jjj=jjj-1
    jjj=0
    average_P_in_x=np.array(my_ray['pressure'][srt])/1.e6
  
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
    temp0=np.array(my_ray['Temp'][srt])
    pressure0=np.array(my_ray['pressure'][srt])


    #plt.clf()
    #if j<=number_plt:
  if (time>=float(jj)*time_difference_T_P and time<max_time) or (j+1==number_plt):
    jj=jj+1
    plt.plot(average_column_mass_in_x,average_T_in_x,label='t='+filename_time+'s')
    if Dimension==2:
      plt.fill_between(average_P_in_x,max_T_in_x,min_T_in_x,alpha=0.5)
    plt.legend(loc=1)
    print('Adding lines at t='+filename_time+'s for T-P plot')
    plt.xlabel(r'column mass$\/[\mathrm{g}/\mathrm{cm}^2$]')
    plt.ylabel(r'$T\/[\mathrm{K}]$')
    plt.xlim(1e-2,1e5)
    plt.ylim(T_range[0],T_range[1])
    #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    plt.xscale('log')
    #plt.yscale('log')
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
    plt.title('$T$ vs $P$')
    plt.savefig("figure_T_columnmass.png")
  
  j=j+1
print("T-column mass plot average made")
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
    average_T_in_x=np.array(my_ray['Temp'][srt])
    average_P_in_x=np.array(my_ray['pressure'][srt])/1.e6
  
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
    temp0=np.array(my_ray['Temp'][srt])
    pressure0=np.array(my_ray['pressure'][srt])


    #plt.clf()
    #if j<=number_plt:
  if (time>=float(jj)*time_difference_T_P and time<max_time) or (j+1==number_plt):
    jj=jj+1
    plt.plot(average_P_in_x,average_T_in_x,label='t='+filename_time+'s')
    if Dimension==2:
      plt.fill_between(average_P_in_x,max_T_in_x,min_T_in_x,alpha=0.5)
    plt.legend(loc=1)
    print('Adding lines at t='+filename_time+'s for T-P plot')
    plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.ylabel(r'$T\/[\mathrm{K}]$')
    plt.xlim(P_range[0],P_range[1])
    plt.ylim(T_range[0],T_range[1])
    #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    plt.xscale('log')
    #plt.yscale('log')
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
    plt.title('$T$ vs $P$')
    plt.savefig("figure_T_P.png")
  
  j=j+1
print("T-P plot average made")
plt.close()








average_T_in_x=0
average_P_in_x=0
average_KE_in_x=0
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
    average_P_in_x=np.array(my_ray['pressure'][srt])/1.e6
    average_KE_in_x=np.array(my_ray['kineng'][srt])
  
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
      press=np.array(my_ray['pressure'][srt])/1.e6
      kinetic_E=np.array(my_ray['kineng'][srt])
      if l==0:
        max_KE_in_x=kinetic_E
        min_KE_in_x=kinetic_E
      else:
        max_KE_in_x=np.maximum(max_KE_in_x,kinetic_E)
        min_KE_in_x=np.minimum(min_KE_in_x,kinetic_E)
      
      average_KE_in_x=(average_KE_in_x*float(l)+kinetic_E)/(l+1)
      average_P_in_x=(average_P_in_x*float(l)+press)/(l+1)
      l=l+1
  #  average_T_in_x=average_T_in_x/float(l-1)
#  average_P_in_x=average_P_in_x/float(l-1)

  l=0
  if j==0 :
    pressure0=np.array(my_ray['pressure'][srt])
    average_KE_in_x0=average_KE_in_x


#plt.clf()
#if j<=number_plt:
  if (time>=float(jj)*time_difference_T_P and time<max_time and j>0) or (j+1==number_plt):
    jj=jj+1
    plt.plot(average_P_in_x,average_KE_in_x,label='t='+filename_time+'s')
    if Dimension==2:
      plt.fill_between(average_P_in_x,max_KE_in_x,min_KE_in_x,alpha=0.5)
    plt.legend(loc=4)
    print('Adding lines at t='+filename_time+'s for T-P plot')
    plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.ylabel(r'specific $KE$ $[{\rm g/cm s^{2}}$]')
    plt.xlim(P_range[0],P_range[1])
    plt.ylim(1e-8,1e5)
    #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    plt.xscale('log')
    plt.yscale('log')
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LB[0],textposition_LB[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LB[0],textposition_LB[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
    plt.title('$KE$ vs $P$')
    plt.savefig("figure_KE_P.png")
  
  j=j+1
print("K-P plot average made")
plt.close()





average_T_in_x=0
average_P_in_x=0
average_KE_in_x=0
average_HEAT_in_x=0
average_KE_in_x=0
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
    average_P_in_x=np.array(my_ray['pressure'][srt])/1.e6
    average_KE_in_x=np.array(my_ray['kineng'][srt])

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
      press=np.array(my_ray['pressure'][srt])/1.e6
      kinetic_E=np.array(my_ray['kineng'][srt])
      if l==0:
        max_KE_in_x=0.0
        min_KE_in_x=1e100
      else:
        max_KE_in_x=np.maximum(max_KE_in_x,abs(average_KE_in_x-average_KE_in_x0)/average_KE_in_x0)
        min_KE_in_x=np.minimum(min_KE_in_x,abs(average_KE_in_x-average_KE_in_x0)/average_KE_in_x0)
      
      average_KE_in_x=(average_KE_in_x*float(l)+kinetic_E)/(l+1)
      average_P_in_x=(average_P_in_x*float(l)+press)/(l+1)

      l=l+1

  l=0
#Notice this is j=1. At t=0, the kinentic energy is 0, so you need a flie other than plt00000 at time close enough to 0
  if j==1 :
    pressure0=np.array(my_ray['pressure'][srt])/1.e6
    average_KE_in_x0=average_KE_in_x

  diff_KE=abs(average_KE_in_x-average_KE_in_x0)/average_KE_in_x0
# diff_pressure=abs(average_P-average_P_in_x0)

#plt.clf()
#if j<=number_plt:
  if (time>=float(jj)*time_difference_T_P and time<max_time and j>0) or (j+1==number_plt):
    jj=jj+1
    plt.plot(average_P_in_x,diff_KE,label='t='+filename_time+'s')
    
    if Dimension==2:
      plt.fill_between(average_P_in_x,max_KE_in_x,min_KE_in_x,alpha=0.5)
    plt.legend(loc=1)
    print('Adding lines at t='+filename_time+'s for Delta E_th/E_th-P plot')
    plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.ylabel(r'$|\Delta KE/KE(t=0)|$')
    plt.xlim(P_range[0],P_range[1])
    plt.ylim(1e-5,1e8)
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
    plt.title(r'$|\Delta KE/KE(t=0)|$ vs $P$')

    plt.savefig("figure_KE_P_diff.png")
  
  j=j+1
print("Delta KE -P plot average made")
plt.close()








average_T_in_x=0
average_P_in_x=0
average_KE_in_x=0
average_HEAT_in_x=0
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
    average_P_in_x=np.array(my_ray['pressure'][srt])/1.e6
    average_HEAT_in_x=np.array(my_ray['eint_E'][srt])

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
      press=np.array(my_ray['pressure'][srt])/1.e6
      thermal_E=np.array(my_ray['eint_E'][srt])
      if l==0:
        max_HEAT_in_x=thermal_E
        min_HEAT_in_x=thermal_E
      else:
        max_HEAT_in_x=np.maximum(max_HEAT_in_x,thermal_E)
        min_HEAT_in_x=np.minimum(min_HEAT_in_x,thermal_E)
      
      average_HEAT_in_x=(average_HEAT_in_x*float(l)+thermal_E)/(l+1)
      average_P_in_x=(average_P_in_x*float(l)+press)/(l+1)
      l=l+1

  l=0
  if j==0 :
    pressure0=np.array(my_ray['pressure'][srt])
    average_HEAT_in_x0=average_HEAT_in_x


#plt.clf()
#if j<=number_plt:
  if (time>=float(jj)*time_difference_T_P and time<max_time) or (j+1==number_plt):
    jj=jj+1
    plt.plot(average_P_in_x,average_HEAT_in_x,label='t='+filename_time+'s')
    
    if Dimension==2:
      plt.fill_between(average_P_in_x,max_HEAT_in_x,min_HEAT_in_x,alpha=0.5)
    plt.legend(loc=1)
    print('Adding lines at t='+filename_time+'s for E_th-P plot')
    plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.ylabel(r'Thermal energy $E_{\rm th}$ $[{\rm erg}$]')
    plt.xlim(P_range[0],P_range[1])
    plt.ylim(1e11,2.0e11)
    #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    plt.xscale('log')
    #plt.yscale('log')
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
    plt.title(r'Thermal energy $E_{\rm{ th}}$ vs $P$')
    plt.savefig("figure_thermalE_P.png")
  
  j=j+1
print("Thermal E -P plot average made")
plt.close()






average_T_in_x=0
average_P_in_x=0
average_KE_in_x=0
average_HEAT_in_x=0

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
    average_P_in_x=np.array(my_ray['pressure'][srt])/1.e6
    press=np.array(my_ray['pressure'][srt])/1.e6
    average_HEAT_in_x=np.array(my_ray['eint_E'][srt])

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
      press=np.array(my_ray['pressure'][srt])/1.e6
      thermal_E=np.array(my_ray['eint_E'][srt])
      if l==0:
        max_HEAT_in_x=0.0
        min_HEAT_in_x=1e100
      else:
        max_HEAT_in_x=np.maximum(max_HEAT_in_x,abs(average_HEAT_in_x-average_HEAT_in_x0)/average_HEAT_in_x0)
        min_HEAT_in_x=np.minimum(min_HEAT_in_x,abs(average_HEAT_in_x-average_HEAT_in_x0)/average_HEAT_in_x0)
      
      average_HEAT_in_x=(average_HEAT_in_x*float(l)+thermal_E)/(l+1)
      average_P_in_x=(average_P_in_x*float(l)+press)/(l+1)

      l=l+1

  l=0
  if j==0 :
    temp0=np.array(my_ray['Temp'][srt])
    ent0=np.array(my_ray['entropy'][srt])
    pressure0=np.array(my_ray['pressure'][srt])/1.e6
    density0=np.array(my_ray['density'][srt])
    average_HEAT_in_x0=average_HEAT_in_x

  diff_HEAT=abs(average_HEAT_in_x-average_HEAT_in_x0)/average_HEAT_in_x0
# diff_pressure=abs(average_P-average_P_in_x0)
  if(j>0):
    HEAT_01_list= [ii for ii,v in enumerate(diff_HEAT) if v>=0.001]
    HEAT_1_list= [ii for ii,v in enumerate(diff_HEAT) if v>=0.01]
    HEAT_10_list= [ii for ii,v in enumerate(diff_HEAT) if v>=0.1]
    for ll in HEAT_01_list:
      HEAT_01=np.maximum(press[ll],HEAT_01)
    for ll in HEAT_1_list:
      HEAT_1=np.maximum(press[ll],HEAT_1)
   
    for ll in HEAT_10_list :
      HEAT_10=np.maximum(press[ll],HEAT_10)
    
    print HEAT_01,HEAT_1,HEAT_10
  
#if j<=number_plt:
  if (time>=float(jj)*time_difference_T_P and time<max_time and j>0) or (j+1==number_plt):
    jj=jj+1
    plt.plot(average_P_in_x,diff_HEAT,label='t='+filename_time+'s')
    
    if Dimension==2:
      plt.fill_between(average_P_in_x,max_HEAT_in_x,min_HEAT_in_x,alpha=0.5)
    plt.legend(loc=1)
    print('Adding lines at t='+filename_time+'s for Delta E_th/E_th-P plot')
    plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.ylabel(r'$|\Delta E_{\rm th}/E_{\rm th}(t=0)|$]')
    plt.xlim(P_range[0],P_range[1])
    plt.ylim(1e-8,0.1)
    #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    plt.xscale('log')
    plt.yscale('log')
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LB[0],textposition_LB[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LB[0],textposition_LB[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
    plt.title(r'$|\Delta E_{\rm{ th}}/E_{\rm{ th}}(t=0)|$ vs $P$')

    plt.savefig("figure_thermalE_P_diff.png")
  
  j=j+1
print("Delta Thermal E -P plot average made")
plt.close()





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
  
  
  
    #
    
    
  kappa=0.18*1e-9
  SB_constant=5.6704*1e-5
    
    
  if Dimension==1:
    my_ray = ds.ortho_ray(0,(0,0))
    srt=np.argsort(my_ray['x'])
    x_coord=np.array(my_ray['x'][srt])/1.e5
    average_T_in_x=np.array(my_ray['Temp'][srt])
    average_P_in_x=np.array(my_ray['pressure'][srt])/1.e6

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

  l=0

  average_T=(average_T*float(j)+average_T_in_x)/float(j+1)
  average_P=(average_P*float(j)+average_P_in_x)/float(j+1)

  
  if j==0 :
    temp0=np.array(my_ray['Temp'][srt])
    pressure0=np.array(my_ray['pressure'][srt])/1.e6
    average_T_in_x0= average_T_in_x
    average_P_in_x0= average_P_in_x
  
    #plt.clf()
    #if j<=number_plt:
  if (time>=float(jj)*time_difference_T_P and time<max_time ) or (j+1==number_plt):
    jj=jj+1
    plt.plot(average_P,average_T,label='t='+filename_time+'s')
    
    if Dimension==2:
      plt.fill_between(average_P_in_x,max_T_in_x,min_T_in_x,alpha=0.5)
    plt.legend(loc=1)
    print('Adding lines at t='+filename_time+'s for <T>-<P> plot')
    plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.ylabel(r'$T\/[\mathrm{K}]$')
    plt.xlim(P_range[0],P_range[1])
    plt.ylim(T_range[0],T_range[1])
    #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    plt.xscale('log')
    #plt.yscale('log')
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
    plt.title('<$T$> vs <$P$> (average up to $t$)')
    plt.savefig("figure_T_P_average.png")
  
  j=j+1
print("T-P plot average made")
plt.close()


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
  
  
  
  #
  
  
  kappa=0.18*1e-9
  SB_constant=5.6704*1e-5
  
  
  if Dimension==1:
    my_ray = ds.ortho_ray(0,(0,0))
    srt=np.argsort(my_ray['x'])
    x_coord=np.array(my_ray['x'][srt])/1.e5
    average_T_in_x=np.array(my_ray['Temp'][srt])
    average_P_in_x=np.array(my_ray['pressure'][srt])/1.e6

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

  l=0

  average_T=(average_T*float(j)+average_T_in_x)/float(j+1)
  average_P=(average_P*float(j)+average_P_in_x)/float(j+1)
  
  
  if j==0 :
    temp0=np.array(my_ray['Temp'][srt])
    pressure0=np.array(my_ray['pressure'][srt])/1.e6
    average_T_in_x0=average_T_in_x
    average_P_in_x0=average_P_in_x

  diff_temp=abs(average_T-average_T_in_x0)/average_T_in_x0
  diff_pressure=abs(average_P-average_P_in_x0)

  #plt.clf()
  #if j<=number_plt:
  if (time>=float(jj)*time_difference_T_P and time<max_time ) or (j+1==number_plt):
    jj=jj+1
    plt.plot(average_P,diff_temp,label='t='+filename_time+'s')

    plt.legend(loc=1)
    print('Adding lines at t='+filename_time+'s for |<Delta_T>/<T>|-P plot')
    plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.ylabel(r'$|\frac{<\Delta T>}{T_{0}}|$')
    plt.xlim(P_range[0],P_range[1])
    #plt.ylim(1400,1500)
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
    #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    plt.xscale('log')
    plt.yscale('log')
    plt.title(r'|$\frac{<\Delta T>}{<T>}$| vs <$P$> (average up to $t$)')
    plt.savefig("figure_T_P_diff_average.png")
  
  j=j+1
print("T-P difference average plot made")
plt.close()






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
  if time>=float(jj)*time_difference_T_P and time<max_time:
    jj=jj+1
    print('Adding lines at t='+filename_time+'s for delta_T - P plot')
    
    
    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5
    
    
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
      average_T_in_x=np.array(my_ray['Temp'][srt])
      average_P_in_x=np.array(my_ray['pressure'][srt])/1.e6

    elif Dimension==2:
    #    my_ray = ds.ortho_ray(1,(0,0))
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
        average_T_in_x=(average_T_in_x*float(l)+temp)/(l+1)
        average_P_in_x=(average_P_in_x*float(l)+press)/(l+1)

        if l==0:
          max_T_in_x=0.0
          min_T_in_x=1e10
        else:
          max_T_in_x=np.maximum(max_T_in_x,abs(average_T_in_x-average_T_in_x0)/average_T_in_x0)
          min_T_in_x=np.minimum(min_T_in_x,abs(average_T_in_x-average_T_in_x0)/average_T_in_x0)

        l=l+1
          #   average_T_in_x=average_T_in_x/float(l)
          #average_P_in_x=average_P_in_x/float(l)

    l=0

    if j==0 :
      temp0=np.array(my_ray['Temp'][srt])
      pressure0=np.array(my_ray['pressure'][srt])
      average_T_in_x0=average_T_in_x
  
    diff_temp=abs(average_T_in_x-average_T_in_x0)/average_T_in_x0
    
    if(j>0):
      T_01_list= [ii for ii,v in enumerate(diff_temp) if v>=0.001]
      T_1_list= [ii for ii,v in enumerate(diff_temp) if v>=0.01]
      T_10_list= [ii for ii,v in enumerate(diff_temp) if v>=0.1]
      for ll in T_01_list:
        T_01=np.maximum(press[ll],T_01)
      
      for ll in T_1_list:
        T_1=np.maximum(press[ll],T_1)
   
      for ll in T_10_list :
        T_10=np.maximum(press[ll],T_10)
    
      print T_01,T_1,T_10

    #   diff_temp_max=abs(max_T_in_x-average_T_in_x0)/average_T_in_x
    # diff_temp_min=abs(min_T_in_x-average_T_in_x0)/average_T_in_x
    #diff_temp_error=np.maximum(diff_temp_max,diff_temp_min)-diff_temp
    print('Adding lines at t='+filename_time+'s for |Delta_T/T|-P plot')

    #plt.clf()
    plt.plot(average_P_in_x,diff_temp,label='t='+filename_time+'s')
    plt.legend(loc=1)
    
    if Dimension==2:
      plt.fill_between(average_P_in_x,max_T_in_x,min_T_in_x,alpha=0.5)

    plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
    plt.ylabel(r'$|\frac{\Delta T}{T_{0}}| $')
    plt.xlim(P_range[0],P_range[1])
    plt.ylim(1e-7,1)
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
    plt.title(r'$|\frac{\Delta T}{T}|$ vs $P$')
    #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    plt.xscale('log')
    plt.yscale('log')

    plt.savefig("figure_T_P_diff.png")
  j=j+1
print("T-P_diff plot made")
plt.close()




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
  if time>=float(jj)*time_difference_x_T and time<max_time:
    jj=jj+1
    print('Adding lines at t='+filename_time+'s for T-height plot')
    location1=ds.find_max('Temp')
    location2=ds.find_min('Temp')
  
    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5
  
  
    ad=ds.all_data()
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])

    plt.semilogy(x_coord,temp,label='t='+filename_time+'s')
    plt.legend()
    plt.ylabel(r'$T$ [K]')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
    plt.title('T vs height')
    plt.xscale('log')
  #plt.legend(loc=3)
  #plt.yscale('linear')
    plt.savefig("figure_x_T_full.png")
  
  
  j=j+1

print("T plot made")
plt.close()


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
  if time>=float(jj)*time_difference_x_opacity and time<max_time:
    jj=jj+1
    print('Adding lines at t='+filename_time+'s for opacity-height plot')
    location1=ds.find_max('Temp')
    location2=ds.find_min('Temp')
    
    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5
    
    
    ad=ds.all_data()
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    opacity=np.array(my_ray['pressure'][srt])*0.18/1.e9
    plt.semilogy(x_coord,opacity,label='t='+filename_time+'s')
    plt.legend()
    plt.ylabel(r'$\kappa$ [$\mathrm{cm}^{2}\mathrm{g}^{-1}$]')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
    plt.title('Opacity $\kappa$ vs height')
    plt.xscale('log')
    #plt.legend(loc=3)
    #plt.yscale('linear')
    plt.savefig("figure_x_opacity.png")
  
  
  j=j+1

print("opacity plot made")
plt.close()


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
  if time>=float(jj)*time_difference_x_Mach and time<max_time:
    jj=jj+1
    print('Adding lines at t='+filename_time+'s for Mach number-height plot')
    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5
    
    
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
      average_P_in_x=np.array(my_ray['pressure'][srt])/1.e6
      average_Mach_in_x=np.array(my_ray['MachNumber'][srt])

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
        press=np.array(my_ray['pressure'][srt])/1.e6
        Mach_number=np.array(my_ray['MachNumber'][srt])
        if l==0:
          max_Mach_in_x=Mach_number
          min_Mach_in_x=Mach_number
        else:
          max_Mach_in_x=np.maximum(max_Mach_in_x,Mach_number)
          min_Mach_in_x=np.minimum(min_Mach_in_x,Mach_number)

        average_Mach_in_x=(average_Mach_in_x*float(l)+Mach_number)/(l+1)
        average_P_in_x=(average_P_in_x*float(l)+press)/(l+1)
        l=l+1
#  average_T_in_x=average_T_in_x/float(l-1)
#  average_P_in_x=average_P_in_x/float(l-1)

    l=0


    
    if j==0 :
      pressure0=np.array(my_ray['pressure'][srt])

    #plt.clf()
    plt.plot(press,average_Mach_in_x,label='t='+filename_time+'s')
    plt.legend(loc=1)
    
    if Dimension==2:
      plt.fill_between(average_P_in_x,max_Mach_in_x,min_Mach_in_x,alpha=0.5)

    plt.ylabel(r'Mach number')
    plt.xlabel(r'$P\/[\mathrm{bar}]$')
    #plt.ylim(location1,location2)
    #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
    #  plt.xscale('linear')
    plt.xlim(P_range[0],P_range[1])
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LB[0],textposition_LB[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LB[0],textposition_LB[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
    plt.ylim(1e-7,100)

    plt.axhline(y=1.0, color='black', linestyle=':')
    plt.yscale('log')
    plt.xscale('log')
    plt.title(r'Mach number $\mathcal{M}$ vs P')
    plt.savefig("figure_P_Mach.png")


  j=j+1


print("Mach plot made")
plt.close()


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
  if time>=float(jj)*time_difference_x_opacity and time<max_time:
    jj=jj+1
    print('Adding lines at t='+filename_time+'s for opacity-T plot')
    location1=ds.find_max('Temp')
    location2=ds.find_min('Temp')
    
    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5
    
    
    ad=ds.all_data()
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    temp=np.array(my_ray['Temp'][srt])
    opacity=np.array(my_ray['pressure'][srt])*0.18/1.e9
    plt.semilogy(temp,opacity,label='t='+filename_time+'s')
    plt.legend()
    plt.ylabel(r'$\kappa$ [$\mathrm{cm}^{2}\mathrm{g}^{-1}$]')
    plt.xlabel(r'$T\/[\mathrm{K}]$')
    #plt.legend(loc=3)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(400,10000)
    plt.ylim(-100,1000)
    plt.title('Opacity $\kappa $vs T')
    plt.savefig("figure_T_opacity.png")
  
  
  j=j+1

print("opacity-T plot made")
plt.close()


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
  if time>=float(jj)*time_difference_x_entropy and time<max_time:
    jj=jj+1
    print('Adding lines at t='+filename_time+'s for entropy-height plot')
    location1=ds.find_max('Temp')
    location2=ds.find_min('Temp')
    
    kappa=0.18*1e-9
    SB_constant=5.6704*1e-5
    
    
    ad=ds.all_data()
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
    elif Dimension==2:
      my_ray = ds.ortho_ray(1,(0,0))
      srt=np.argsort(my_ray['y'])
      x_coord=np.array(my_ray['y'][srt])/1.e5
    ent=np.array(my_ray['entropy'][srt])
    
    plt.semilogy(x_coord,ent,label='t='+filename_time+'s')
    plt.legend()
    plt.ylabel(r'$\mathrm{Entropy}$')
    plt.xlabel(r'$y\/[\mathrm{km}]$')
    plt.xscale('log')
    text_num=0
    while (text_num<=3):
      if textbox[text_num]==0:
        plt.figtext(textposition_LB[0],textposition_LB[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
      elif textbox[text_num]==1:
        plt.figtext(textposition_LB[0],textposition_LB[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
      text_num=text_num+1
    #plt.legend(loc=3)
    #plt.yscale('linear')
    plt.title('Entropy vs height')
    plt.savefig("figure_x_entropy.png")
  
  
  j=j+1

print("entropy plot made")
plt.close()



if radiation_field==1:
  average_T_in_x=0
  average_P_in_x=0
  average_KE_in_x=0
  average_U_in_x=0
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
  
  
    if Dimension==1:
      my_ray = ds.ortho_ray(0,(0,0))
      srt=np.argsort(my_ray['x'])
      x_coord=np.array(my_ray['x'][srt])/1.e5
      average_P_in_x=np.array(my_ray['pressure'][srt])/1.e6
      average_U_in_x=np.array(my_ray['rad'][srt])

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
        press=np.array(my_ray['pressure'][srt])/1.e6
        radiation=np.array(my_ray['rad'][srt])
        if l==0:
          max_U_in_x=radiation
          min_U_in_x=radiation
        else:
          max_U_in_x=np.maximum(max_U_in_x,radiation)
          min_U_in_x=np.minimum(min_U_in_x,radiation)
      
        average_U_in_x=(average_U_in_x*float(l)+radiation)/(l+1)
        average_P_in_x=(average_P_in_x*float(l)+press)/(l+1)
        l=l+1

    l=0
    if j==0 :
      pressure0=np.array(my_ray['pressure'][srt])


#plt.clf()
#if j<=number_plt:
    if ( time>=float(jj)*time_difference_T_P and time<max_time) or (j+1==number_plt):
      jj=jj+1
      plt.plot(average_P_in_x,average_U_in_x,label='t='+filename_time+'s')
    
      if Dimension==2:
        plt.fill_between(average_P_in_x,max_U_in_x,min_U_in_x,alpha=0.5)
      plt.legend(loc=4)
      print('Adding lines at t='+filename_time+'s for T-P plot')
      plt.xlabel(r'$P\/[\mathrm{bar}=10^{6}\mathrm{dyne}/\mathrm{cm}^2$]')
      plt.ylabel(r'$U$ $[{\rm erg}/{\rm cm}^{3}]$')
      plt.xlim(P_range[0],P_range[1])
     # plt.ylim(1e-3,1e-1)
    #plt.xaxis.set_major_formatter(FormatStrFormatter('10^{%T}'))
      plt.xscale('log')
      plt.yscale('symlog')
      text_num=0
      while (text_num<=3):
        if textbox[text_num]==0:
          plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
        elif textbox[text_num]==1:
          plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
        text_num=text_num+1
      plt.title('$U$ vs $P$')
      plt.savefig("figure_U_P.png")
  
    j=j+1
  print("U-P plot average made")
  plt.close()





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
      plt.yscale('symlog')
      plt.xscale('log')
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

print "P at HEAT_0.1%",HEAT_01
print "P at HEAT_1%",HEAT_1
print "P at HEAT_10%",HEAT_10
print "P at T_01%",T_01
print "P at T_1%",T_1
print "P at T_10%",T_10
