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
import colorsys
np.set_printoptions(threshold=np.inf)


def get_color_red(color):
  for hue in range(color):
    hue=1.0*hue/color
    col=[int(x) for x in colorsys.hsv_to_rgb(hue,1.0,230)]
    yield "#{0:02x}{1:02x}{2:02x}".format(*col)

def get_color_blue(color):
  for hue in range(color):
    hue=1.0*hue/color
    col=[int(x) for x in colorsys.hsv_to_rgb(1.0,hue,230)]
    yield "#{0:02x}{1:02x}{2:02x}".format(*col)



color_red=get_color_red(2)
color_blue=get_color_blue(2)

print (color_red)
modelname=['fixed $\Delta t=10^{-6}$',\
          r'$\sigma_{\mathrm{vol}}/y=5\sigma_{0}$',\
          r'$A_{\mathrm{vol}}/y=A_{0}$',\
          r' $N_{\mathrm{vol}}=100$   ']
color=['black','black','black','black']
bbox_props = [dict(facecolor='none', edgecolor=color[0], pad=0.8),\
              dict(facecolor='none', edgecolor=color[1], pad=0.8),\
              dict(facecolor='none', edgecolor=color[2], pad=0.8),\
              dict(facecolor='none', edgecolor=color[3], pad=0.8)]
#box=1, no box=0
fluxlabel=['SF=102','SF=103','SF=104','SF=1002','SF=1003(F)']
linetype=['-','--','-.',':']
linethickness=[0.5,1,1,1.5]
linec=['k','k','k','k','k','k','k','k']#,'orange','r','r']
fontsize1=9
textbox=[0,0,0,0]
#Left top corner
textposition_LT=[0.15,0.83]
#Right top corner
textposition_RT=[0.45,0.83]
#Left bottom corner
textposition_LB=[0.15,0.3]
#Right bottom corner
textposition_RB=[0.45,0.3]

def index_exists(ls, i):
    return (0 <= i < len(ls)) or (-len(ls) <= i < 0)


numberoflines=5.0
Dimension=1
number_slice_x=100
x_range=[0.0,1.0e9]
P_range=[0.001,2.e2]
T_range=[2300,3000]
rho_range=[1e-10,0.1]
l=0
i=0
j=0
jj=0
Dimension = int(raw_input("Dimension :"))
numberoflines = int(raw_input("numberoflines :"))
every_nthfile=3
#print 'Listing all profile/dat files'
profilefilelist = sorted(glob.glob('plt*'), key=lambda name: int(name[3:]))
#raw_input('Press ENTER to continue...')
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

number_plt=0
time=0
radiation_field=1







for i in profilefilelist:
    #  if np.mod(j,every_nthfile)==0:
  ds=load(i)
  number_plt=number_plt+1
  max_time=np.maximum(time,float(ds.current_time))
  if number_plt==1:
    ds.index
    for j in sorted(ds.field_list):
      if j[1]=="rad":
        radiation_field=1
        print "rad"
    print i
    print "radiation field", radiation_field
    my_ray_x = ds.ortho_ray(0,(0,0))
    srt_x=np.argsort(my_ray_x['x'])
    my_ray_y = ds.ortho_ray(1,(0,0))
    srt_y=np.argsort(my_ray_y['y'])
    number_cell_height=float(len(np.array(my_ray_x['x'][srt_x])))
    number_cell_width=float(len(np.array(my_ray_y['y'][srt_y])))
    if Dimension==1:
      print ("1Dimension")
      y_range=[np.min(np.array(my_ray_x['x'][srt_x])),np.amax(np.array(my_ray_x['x'][srt_x]))]
      T_range=[np.min(np.array(my_ray_x['Temp'][srt_x])),np.amax(np.array(my_ray_x['Temp'][srt_x]))]
      rho_range=[np.min(np.array(my_ray_x['density'][srt_x])),np.amax(np.array(my_ray_x['density'][srt_x]))]
      P_range=[np.min(np.array(my_ray_x['pressure'][srt_x])),np.amax(np.array(my_ray_x['pressure'][srt_x]))]
    elif Dimension==2:
      print ("2Dimension")
      x_range=[np.min(np.array(my_ray_x['x'][srt_x])),np.amax(np.array(my_ray_x['x'][srt_x]))]
      
      y_range=[np.min(np.array(my_ray_y['y'][srt_y])),np.amax(np.array(my_ray_y['y'][srt_y]))]
      T_range=[np.min(np.array(my_ray_y['Temp'][srt_y])),np.amax(np.array(my_ray_y['Temp'][srt_y]))]
      rho_range=[np.min(np.array(my_ray_y['density'][srt_y])),np.amax(np.array(my_ray_y['density'][srt_y]))]
      P_range=[np.min(np.array(my_ray_y['pressure'][srt_y])),np.amax(np.array(my_ray_y['pressure'][srt_y]))]
    delta_x=(y_range[1]-y_range[0])/number_cell_height
#P_range=[x for x in P_range] #[bar]
#x_range=[x for x in x_range] #[km]
#y_range=[x for x in y_range] #[km]
print "preliminary step"
print "x range [km] :",x_range[0]/1e5,x_range[1]/1e5
print "y range [km] :",y_range[0]/1e5,y_range[1]/1e5
print "cell number  :",number_cell_height
print "T range [T]  :",T_range
print "K range [bar]:",P_range[0]/1e6,P_range[1]/1e6
print "rho range [bar]:",rho_range[0],rho_range[1]
print "Are you sure about the dimension?:",Dimension
time_difference_x_den=max_time/numberoflines
time_difference_x_P=max_time/numberoflines
time_difference_x_T=max_time/numberoflines
time_difference_x_entropy=max_time/numberoflines
time_difference_x_rad=max_time/numberoflines
time_difference_T_P=max_time/numberoflines
time_difference_x_opacity=max_time/numberoflines
time_difference_x_Mach=max_time/numberoflines
time_difference_v_P=max_time/numberoflines


plot_type = np.zeros(20)
print "    X-axis -   Y-axis"
print "---------------------"
print "1 : height - rho plot"
print "2 : height - P plot"
print "3 : height - T plot"
print "4 : height - v_x plot"
print "5 : height - v_y plot"
print "6 : T      - P plot"
print "7 : v_x    - P plot"
print "8 : v_y    - P plot"
print "9 : particle position - t  plot"
print "10: Favre average turbulent kinetic energy  - P plot "
plot_type = [int(n) for n in raw_input('Enter numbers: ').split()]
print 'plot_type',plot_type
  #for i in range(int(num)):
  #n = raw_input("num :")
#  plot_type.append(int(n))
#print 'ARRAY: ',plot_type

#print profilefilelist
for plot_type_index in plot_type :
  if plot_type_index==1:
    print  "1 : height - rho plot"
    average_P_in_x=0
    average_rho_in_x=0
    j=0
    jj=0

    for i in profilefilelist:
      l=0

      ds=load(i)
      time=float(ds.current_time)
      filename_time=str(time)
      if time>=float(jj)*time_difference_x_den and time<=max_time:
        jj=jj+1
        print('Adding lines at t='+filename_time+'s for heigt-rho plot')
        my_ray_x = ds.ortho_ray(0,(0,0))
        srt_x=np.argsort(my_ray_x['x'])
        my_ray_y = ds.ortho_ray(1,(0,0))
        srt_y=np.argsort(my_ray_y['y'])

        if Dimension==1:
          x_coord=np.array(my_ray_x['x'][srt_x])
          rho=np.array(my_ray_x['density'][srt_x])
          average_rho_in_x=rho
        elif Dimension==2:
          x_coord=np.array(my_ray_y['y'][srt_y])
          while (l<number_slice_x):
            position_x=x_range[0]+float(l)*(x_range[1]-x_range[0])/float(number_slice_x)
            my_ray_y = ds.ortho_ray(1,(0,position_x))
            srt_y=np.argsort(my_ray_y['y'])
            rho=np.array(my_ray_y['density'][srt_y])
    
            if l==0:
              max_rho_in_x=rho
              min_rho_in_x=rho
            else:
              max_rho_in_x=np.maximum(max_rho_in_x,rho)
              min_rho_in_x=np.minimum(min_rho_in_x,rho)

            average_rho_in_x=(average_rho_in_x*float(l)+rho)/float(l+1)
            l=l+1

        if j==0 :
          rho0=rho

        plt.plot(x_coord/1e5,average_rho_in_x,label='t='+filename_time[:6]+'s')#
              #,label=fluxlabel[j],linestyle=linetype[j], linewidth=linethickness[j], color=linec[j])
        plt.legend(loc=1)#,prop={'size':fontsize1})
        plt.ylabel(r'$\rho\/[\mathrm{g}/\mathrm{cm}^{3}]$')
        plt.xlabel(r'$y\/[\mathrm{km}]$')
        plt.xlim(y_range[0]*0.9/1e5,y_range[1]*1.1/1e5)
        plt.ylim(rho_range[0]*0.9,rho_range[1]*1.1)
        plt.yscale('log')
      #plt.xscale('log')
        text_num=0
        while (text_num<=3):
          if textbox[text_num]==0:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
          elif textbox[text_num]==1:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
          text_num=text_num+1
        plt.title(r'density $\rho$ vs distance $h$')
        plt.savefig("figure_h_density.png")
      j=j+1

    print("density-height plot made")
    plt.close()


  if plot_type_index==2:
    print  "1 : height - P plot"
    average_press_in_x=0
    average_rho_in_x=0
      
    jj=0
    j=0
    for i in profilefilelist:

      l=0

      ds=load(i)
      time=float(ds.current_time)
      filename_time=str(time)
      if time>=float(jj)*time_difference_x_P and time<=max_time:
        my_ray_x = ds.ortho_ray(0,(0,0))
        srt_x=np.argsort(my_ray_x['x'])
        my_ray_y = ds.ortho_ray(1,(0,0))
        srt_y=np.argsort(my_ray_y['y'])


        jj=jj+1
        print('Adding lines at t='+filename_time+'s for heigt-rho plot')
        if Dimension==1:
          x_coord=np.array(my_ray_x['x'][srt_x])
          press=np.array(my_ray_x['pressure'][srt_x])
          average_press_in_x=press
        elif Dimension==2:
          x_coord=np.array(my_ray_y['y'][srt_y])
          while (l<number_slice_x):
            position_x=x_range[0]+float(l)*(x_range[1]-x_range[0])/float(number_slice_x)
            my_ray_y = ds.ortho_ray(1,(0,position_x))
            srt_y=np.argsort(my_ray_y['y'])
            my_ray_y = ds.ortho_ray(1,(0,position_x))
            srt_y=np.argsort(my_ray_y['y'])

            press=np.array(my_ray_y['pressure'][srt_y])
    
            if l==0:
              max_press_in_x=press
              min_press_in_x=press
            else:
              max_press_in_x=np.maximum(max_press_in_x,press)
              min_press_in_x=np.minimum(min_press_in_x,press)

            average_press_in_x=(average_press_in_x*float(l)+press)/float(l+1)
            l=l+1

        if j==0 :
          press0=press

        plt.plot(x_coord/1e5,average_press_in_x/1e6,label='t='+filename_time[:6]+'s')#
              #,label=fluxlabel[j],linestyle=linetype[j], linewidth=linethickness[j], color=linec[j])
        plt.legend(loc=1)#,prop={'size':fontsize1})
        plt.ylabel(r'$P\/[\mathrm{bar}]$')
        plt.xlabel(r'$y\/[\mathrm{km}]$')
        plt.xlim(y_range[0]*0.9/1e5,y_range[1]*1.1/1e5)
        plt.ylim(P_range[0]*0.9/1e6,P_range[1]*1.1/1e6)
        plt.yscale('log')
      #plt.xscale('log')
        text_num=0
        while (text_num<=3):
          if textbox[text_num]==0:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
          elif textbox[text_num]==1:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
          text_num=text_num+1
        plt.title(r'Pressure $P$ vs distance $h$')
        plt.savefig("figure_h_pressure.png")
      j=j+1

    print("pressure-height plot made")
    plt.close()




  if plot_type_index==3:
    print  "1 : height - T plot"
    average_temperature_in_x=0
    jj=0
    j=0
    for i in profilefilelist:
      ds=load(i)
      time=float(ds.current_time)
      filename_time=str(time)
      l=0

      if time>=float(jj)*time_difference_x_T and time<=max_time:
        my_ray_x = ds.ortho_ray(0,(0,0))
        srt_x=np.argsort(my_ray_x['x'])
        my_ray_y = ds.ortho_ray(1,(0,0))
        srt_y=np.argsort(my_ray_y['y'])


        jj=jj+1
        print('Adding lines at t='+filename_time+'s for heigt-rho plot')
        if Dimension==1:
          x_coord=np.array(my_ray_x['x'][srt_x])
          temperature=np.array(my_ray_x['Temp'][srt_x])
          average_temperature_in_x=temperature
        elif Dimension==2:
          x_coord=np.array(my_ray_y['y'][srt_y])
          while (l<number_slice_x):
            position_x=x_range[0]+float(l)*(x_range[1]-x_range[0])/float(number_slice_x)
            my_ray_y = ds.ortho_ray(1,(0,position_x))
            srt_y=np.argsort(my_ray_y['y'])
            temperature=np.array(my_ray_y['Temp'][srt_y])
    
            if l==0:
              max_temperature_in_x=temperature
              min_temperature_in_x=temperature
            else:
              max_temperature_in_x=np.maximum(max_temperature_in_x,temperature)
              min_temperature_in_x=np.minimum(min_temperature_in_x,temperature)

            average_temperature_in_x=(average_temperature_in_x*float(l)+temperature)/float(l+1)
            l=l+1

        if j==0 :
          temperature0=temperature

        plt.plot(x_coord/1e5,average_temperature_in_x,label='t='+filename_time[:6]+'s')#
              #,label=fluxlabel[j],linestyle=linetype[j], linewidth=linethickness[j], color=linec[j])
        plt.legend(loc=1)#,prop={'size':fontsize1})
        plt.ylabel(r'$T\/[\mathrm{K}]$')
        plt.xlabel(r'$y\/[\mathrm{km}]$')
        plt.xlim(y_range[0]*0.9/1e5,y_range[1]*1.1/1e5)
        plt.ylim(T_range[0]*0.9,T_range[1]*1.1)
        plt.yscale('linear')
      #plt.xscale('log')
        text_num=0
        while (text_num<=3):
          if textbox[text_num]==0:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
          elif textbox[text_num]==1:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
          text_num=text_num+1
        plt.title(r'temperature $T$ vs distance $h$')
        plt.savefig("figure_h_T.png")
      j=j+1

    print("T-height plot made")
    plt.close()



  if plot_type_index==4:
    print  "1 : height - v_x plot"
    average_vx_in_x=0
    jj=0
    j=0
    for i in profilefilelist:
      ds=load(i)
      time=float(ds.current_time)
      filename_time=str(time)
      l=0

      if time>=float(jj)*time_difference_x_T and time<=max_time:
        my_ray_x = ds.ortho_ray(0,(0,0))
        srt_x=np.argsort(my_ray_x['x'])
        my_ray_y = ds.ortho_ray(1,(0,0))
        srt_y=np.argsort(my_ray_y['y'])


        jj=jj+1
        print('Adding lines at t='+filename_time+'s for heigt-rho plot')
        if Dimension==1:
          x_coord=np.array(my_ray_x['x'][srt_x])
          vx=np.array(my_ray_x['x_velocity'][srt_x])
          average_vx_in_x=vx
        elif Dimension==2:
          x_coord=np.array(my_ray_y['y'][srt_y])
         
          while (l<number_slice_x):
            position_x=x_range[0]+float(l)*(x_range[1]-x_range[0])/float(number_slice_x)
            my_ray_y = ds.ortho_ray(1,(0,position_x))
            srt_y=np.argsort(my_ray_y['y'])
            vx=np.array(my_ray_y['x_velocity'][srt_y])
    
            if l==0:
              max_vx_in_x=vx
              min_vx_in_x=vx
            else:
              max_vx_in_x=np.maximum(max_vx_in_x,vx)
              min_vx_in_x=np.minimum(min_vx_in_x,vx)

            average_vx_in_x=(average_vx_in_x*float(l)+vx)/float(l+1)
            l=l+1

        if j==0 :
          vx0=vx

        plt.plot(x_coord/1e5,average_vx_in_x/1e5,label='t='+filename_time[:6]+'s')#
              #,label=fluxlabel[j],linestyle=linetype[j], linewidth=linethickness[j], color=linec[j])
        plt.legend(loc=1)#,prop={'size':fontsize1})
        plt.ylabel(r'$v\/[\mathrm{km/s}]$')
        plt.xlabel(r'$y\/[\mathrm{km}]$')
        plt.xlim(y_range[0]*0.9/1e5,y_range[1]*1.1/1e5)
        #  plt.ylim(T_range[0]*0.9,T_range[1]*1.1)
        plt.yscale('linear')
      #plt.xscale('log')
        text_num=0
        while (text_num<=3):
          if textbox[text_num]==0:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
          elif textbox[text_num]==1:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
          text_num=text_num+1
        plt.title(r'velocity $v_{x}$ vs distance $h$')
        plt.savefig("figure_h_vx.png")
      j=j+1

    print("vx-height plot made")
    plt.close()



  if plot_type_index==5 and Dimension==2:
    print  "1 : height - v_y plot"
    jj=0
    j=0
    average_vy_in_x=0
    for i in profilefilelist:
      l=0

      ds=load(i)
      time=float(ds.current_time)
      filename_time=str(time)
  
      if time>=float(jj)*time_difference_x_T and time<=max_time:
        my_ray_x = ds.ortho_ray(0,(0,0))
        srt_x=np.argsort(my_ray_x['x'])
        my_ray_y = ds.ortho_ray(1,(0,0))
        srt_y=np.argsort(my_ray_y['y'])


        jj=jj+1
        print('Adding lines at t='+filename_time+'s for heigt-rho plot')
        if Dimension==2:

          x_coord=np.array(my_ray_y['y'][srt_y])
          while (l<number_slice_x):
            position_x=x_range[0]+float(l)*(x_range[1]-x_range[0])/float(number_slice_x)
            my_ray_y = ds.ortho_ray(1,(0,position_x))
            srt_y=np.argsort(my_ray_y['y'])
            vy=np.array(my_ray_y['y_velocity'][srt_y])
    
            if l==0:
              max_vy_in_x=vy
              min_vy_in_x=vy
            else:
              max_vy_in_x=np.maximum(max_vy_in_x,vy)
              min_vy_in_x=np.minimum(min_vy_in_x,vy)

            average_vy_in_x=(average_vy_in_x*float(l)+vy)/float(l+1)
            l=l+1

        if j==0 :
          vy0=vy

        plt.plot(x_coord/1e5,average_vy_in_x/1e5,label='t='+filename_time[:6]+'s')#
              #,label=fluxlabel[j],linestyle=linetype[j], linewidth=linethickness[j], color=linec[j])
        plt.legend(loc=1)#,prop={'size':fontsize1})
        plt.ylabel(r'$v\/[\mathrm{km/s}]$')
        plt.xlabel(r'$y\/[\mathrm{km}]$')
        plt.xlim(y_range[0]*0.9/1e5,y_range[1]*1.1/1e5)
        #  plt.ylim(T_range[0]*0.9,T_range[1]*1.1)
        plt.yscale('linear')
      #plt.xscale('log')
        text_num=0
        while (text_num<=3):
          if textbox[text_num]==0:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
          elif textbox[text_num]==1:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
          text_num=text_num+1
        plt.title(r'velocity $v_{y}$ vs distance $h$')
        plt.savefig("figure_h_vy.png")
      j=j+1

    print("vy-height plot made")
    plt.close()





  if plot_type_index==6:
    print  "1 : height - T_P plot"
    average_temperature_in_x=0
    average_press_in_x=0
    jj=0
    j=0

    for i in profilefilelist:
      ds=load(i)
      time=float(ds.current_time)
      filename_time=str(time)
      l=0
 
      if time>=float(jj)*time_difference_T_P and time<=max_time:
        jj=jj+1
        my_ray_x = ds.ortho_ray(0,(0,0))
        srt_x=np.argsort(my_ray_x['x'])
        my_ray_y = ds.ortho_ray(1,(0,0))
        srt_y=np.argsort(my_ray_y['y'])

        print('Adding lines at t='+filename_time+'s for heigt-rho plot')
        if Dimension==1:
          x_coord=np.array(my_ray_x['x'][srt_x])
          temperature=np.array(my_ray_x['Temp'][srt_x])
          press=np.array(my_ray_x['pressure'][srt_x])
          average_temperature_in_x=temperature
          average_press_in_x=press
        if Dimension==2:
          x_coord=np.array(my_ray_y['y'][srt_y])
          

          while (l<number_slice_x):
            position_x=x_range[0]+float(l)*(x_range[1]-x_range[0])/float(number_slice_x)
            my_ray_y = ds.ortho_ray(1,(0,position_x))
            srt_y=np.argsort(my_ray_y['y'])
            temperature=np.array(my_ray_y['Temp'][srt_y])
            press=np.array(my_ray_y['pressure'][srt_y])
            if l==0:
              max_temperature_in_x=temperature
              min_temperature_in_x=temperature
              max_press_in_x=press
              min_press_in_x=press
            else:
              max_temperature_in_x=np.maximum(max_temperature_in_x,temperature)
              min_temperature_in_x=np.minimum(min_temperature_in_x,temperature)
              max_press_in_x=np.maximum(max_press_in_x,press)
              min_press_in_x=np.minimum(min_press_in_x,press)

            average_temperature_in_x=(average_temperature_in_x*float(l)+temperature)/float(l+1)
            average_press_in_x=(average_press_in_x*float(l)+press)/float(l+1)
            l=l+1

        if j==0 :
          temperature0=temperature
          press0=press

        plt.plot(average_press_in_x/1e6,average_temperature_in_x,label='t='+filename_time[:6]+'s')#
              #,label=fluxlabel[j],linestyle=linetype[j], linewidth=linethickness[j], color=linec[j])
        plt.fill_between(average_press_in_x/1e6,max_temperature_in_x,min_temperature_in_x,alpha=0.5)


        plt.legend(loc=1)#,prop={'size':fontsize1})
        plt.ylabel(r'$T\/[\mathrm{K}]$')
        plt.xlabel(r'$P\/[\mathrm{bar}]$')
        plt.xlim(0.001,P_range[1]*1.1/1e6)
        plt.ylim(2800,T_range[1]*1.1)
        plt.yscale('linear')
        plt.xscale('log')
        text_num=0
        while (text_num<=3):
          if textbox[text_num]==0:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
          elif textbox[text_num]==1:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
          text_num=text_num+1
        plt.title(r'Temperature $T$ vs Pressure $P$')
        plt.savefig("figure_T_P.png")
      j=j+1

    print("T-P plot made")
    plt.close()

  if plot_type_index==7:
    print  "1 : height - vx_P plot"
    average_vx_in_x=0
    average_press_in_x=0
    jj=0
    j=0

    for i in profilefilelist:
      ds=load(i)
      time=float(ds.current_time)
      filename_time=str(time)
      l=0
 
      if time>=float(jj)*time_difference_T_P and time<=max_time:
        jj=jj+1
        my_ray_x = ds.ortho_ray(0,(0,0))
        srt_x=np.argsort(my_ray_x['x'])
        my_ray_y = ds.ortho_ray(1,(0,0))
        srt_y=np.argsort(my_ray_y['y'])

        print('Adding lines at t='+filename_time+'s for heigt-rho plot')
        if Dimension==1:
          x_coord=np.array(my_ray_x['x'][srt_x])
          vx=np.array(my_ray_x['x_velocity'][srt_x])
          press=np.array(my_ray_x['pressure'][srt_x])
          average_vx_in_x=vx
          average_press_in_x=press
        if Dimension==2:
          x_coord=np.array(my_ray_y['y'][srt_y])
          

          while (l<number_slice_x):
            position_x=x_range[0]+float(l)*(x_range[1]-x_range[0])/float(number_slice_x)
            my_ray_y = ds.ortho_ray(1,(0,position_x))
            srt_y=np.argsort(my_ray_y['y'])
            vx=np.array(my_ray_y['x_velocity'][srt_y])
            press=np.array(my_ray_y['pressure'][srt_y])
    
            if l==0:
              max_vx_in_x=vx
              min_vx_in_x=vx
              max_press_in_x=press
              min_press_in_x=press

            else:
              max_vx_in_x=np.maximum(max_vx_in_x,vx)
              min_vx_in_x=np.minimum(min_vx_in_x,vx)
              max_press_in_x=np.maximum(max_press_in_x,press)
              min_press_in_x=np.minimum(min_press_in_x,press)

            average_vx_in_x=(average_vx_in_x*float(l)+vx)/float(l+1)
            average_press_in_x=(average_press_in_x*float(l)+press)/float(l+1)

            l=l+1

        if j==0 :
          vx0=vx
          press0=press

        plt.plot(average_press_in_x/1e6,average_vx_in_x/1e5,label='t='+filename_time[:6]+'s')#
              #,label=fluxlabel[j],linestyle=linetype[j], linewidth=linethickness[j], color=linec[j])
        plt.fill_between(average_press_in_x/1e6,max_vx_in_x,min_vx_in_x,alpha=0.5)

        plt.legend(loc=1)#,prop={'size':fontsize1})
        plt.ylabel(r'$v\/[\mathrm{km/s}]$')
        plt.xlabel(r'$P\/[\mathrm{bar}]$')
        plt.xlim(0.001,P_range[1]*1.1/1e6)
        # plt.ylim(2800,T_range[1]*1.1)
        plt.yscale('linear')
        plt.xscale('log')
        text_num=0
        while (text_num<=3):
          if textbox[text_num]==0:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
          elif textbox[text_num]==1:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
          text_num=text_num+1
        plt.title(r'v_{\mathrm{x}} vs Pressure $P$')
        plt.savefig("figure_vx_P.png")
      j=j+1

    print("v_x-P plot made")
    plt.close()


  if plot_type_index==8 and Dimension==2:
    print  "1 : height - vy_P plot"
    average_vy_in_x=0
    average_press_in_x=0
    jj=0
    j=0

    for i in profilefilelist:
      ds=load(i)
      time=float(ds.current_time)
      filename_time=str(time)
      l=0
 
      if time>=float(jj)*time_difference_T_P and time<=max_time:
        jj=jj+1
        my_ray_x = ds.ortho_ray(0,(0,0))
        srt_x=np.argsort(my_ray_x['x'])
        my_ray_y = ds.ortho_ray(1,(0,0))
        srt_y=np.argsort(my_ray_y['y'])

        print('Adding lines at t='+filename_time+'s for heigt-rho plot')
        if Dimension==1:
          x_coord=np.array(my_ray_x['x'][srt_x])
          vy=np.array(my_ray_x['y_velocity'][srt_x])
          press=np.array(my_ray_x['pressure'][srt_x])
          average_vy_in_x=vy
          average_press_in_x=press
        if Dimension==2:
          x_coord=np.array(my_ray_y['y'][srt_y])
          

          while (l<number_slice_x):
            position_x=x_range[0]+float(l)*(x_range[1]-x_range[0])/float(number_slice_x)
            my_ray_y = ds.ortho_ray(1,(0,position_x))
            srt_y=np.argsort(my_ray_y['y'])
            vy=np.array(my_ray_y['y_velocity'][srt_y])
            press=np.array(my_ray_y['pressure'][srt_y])
    
            if l==0:
              max_vy_in_x=vy
              min_vy_in_x=vy
              max_press_in_x=press
              min_press_in_x=press

            else:
              max_vy_in_x=np.maximum(max_vy_in_x,vy)
              min_vy_in_x=np.minimum(min_vy_in_x,vy)
              max_press_in_x=np.maximum(max_press_in_x,press)
              min_press_in_x=np.minimum(min_press_in_x,press)

            average_vy_in_x=(average_vy_in_x*float(l)+vy)/float(l+1)
            average_press_in_x=(average_press_in_x*float(l)+press)/float(l+1)

            l=l+1

        if j==0 :
          vy0=vy
          press0=press

        plt.plot(average_press_in_x/1e6,average_vy_in_x/1e5,label='t='+filename_time[:6]+'s')#
              #,label=fluxlabel[j],linestyle=linetype[j], linewidth=linethickness[j], color=linec[j])
        plt.fill_between(average_press_in_x/1e6,max_vy_in_x,min_vy_in_x,alpha=0.5)


        plt.legend(loc=1)#,prop={'size':fontsize1})
        plt.ylabel(r'$v\/[\mathrm{km/s}]$')
        plt.xlabel(r'$P\/[\mathrm{bar}]$')
        plt.xlim(0.001,P_range[1]*1.1/1e6)
        # plt.ylim(2800,T_range[1]*1.1)
        plt.yscale('linear')
        plt.xscale('log')
        text_num=0
        while (text_num<=3):
          if textbox[text_num]==0:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
          elif textbox[text_num]==1:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
          text_num=text_num+1
        plt.title(r'v_{\mathrm{y}} vs Pressure $P$')
        plt.savefig("figure_vy_P.png")
      j=j+1

    print("v_x-P plot made")
    plt.close()


  if plot_type_index==9 and Dimension ==2:
    print  "1 : height - particle position plot"
    
    
    particlefilelist = sorted(glob.glob('*P_Timestamp_average.txt'), key=lambda name: int(name[:1]))
    
    for i in particlefilelist:
    
      particle_input=open(i,"r")
      particle_data=particle_input.readlines()
    
      number_particle=[]
      t=[]
      ini_particle_rho=[]
      average_particle_rho=[]
      max_particle_rho=[]
      min_particle_rho=[]
      ini_particle_T=[]
      average_particle_T=[]
      max_particle_T=[]
      min_particle_T=[]

      ini_particle_P=[]
      average_particle_P=[]
      max_particle_P=[]
      min_particle_P=[]

      ini_particle_x=[]
      average_particle_x=[]
      max_particle_x=[]
      min_particle_x=[]

      ini_particle_y=[]
      average_particle_y=[]
      max_particle_y=[]
      min_particle_y=[]

      ini_particle_vx=[]
      average_particle_vx=[]
      max_particle_vx=[]
      min_particle_vx=[]

      ini_particle_vy=[]
      average_particle_vy=[]
      max_particle_vy=[]
      min_particle_vy=[]

    
      for line in particle_data:
        p=line.split()
        number_particle.append(float(p[0]))
        t.append(float(p[1]))
        ini_particle_rho.append(float(p[2]))
        average_particle_rho.append(float(p[3]))
        max_particle_rho.append(float(p[4]))
        min_particle_rho.append(float(p[5]))
        ini_particle_T.append(float(p[6]))
        average_particle_T.append(float(p[7]))
        max_particle_T.append(float(p[8]))
        min_particle_T.append(float(p[9]))
        ini_particle_P.append(float(p[10]))
        average_particle_P.append(float(p[11]))
        max_particle_P.append(float(p[12]))
        min_particle_P.append(float(p[13]))
        ini_particle_x.append(float(p[14]))
        average_particle_x.append(float(p[15]))
        max_particle_x.append(float(p[16]))
        min_particle_x.append(float(p[17]))
        ini_particle_y.append(float(p[18]))
        average_particle_y.append(float(p[19]))
        max_particle_y.append(float(p[20]))
        min_particle_y.append(float(p[21]))
        ini_particle_vx.append(float(p[22]))
        average_particle_vx.append(float(p[23]))
        max_particle_vx.append(float(p[24]))
        min_particle_vx.append(float(p[25]))
        ini_particle_vy.append(float(p[26]))
        average_particle_vy.append(float(p[27]))
        max_particle_vy.append(float(p[28]))
        min_particle_vy.append(float(p[29]))
    
      number_particle=np.array(number_particle)
      t=np.array(t)
      ini_particle_rho=np.array(ini_particle_rho)
      average_particle_rho=np.array(average_particle_rho)
      max_particle_rho=np.array(max_particle_rho)
      min_particle_rho=np.array(min_particle_rho)
      ini_particle_T=np.array(ini_particle_T)
      average_particle_T=np.array(average_particle_T)
      max_particle_T=np.array(max_particle_T)
      min_particle_T=np.array(min_particle_T)
      ini_particle_P=np.array(ini_particle_P)
      average_particle_P=np.array(average_particle_P)
      max_particle_P=np.array(max_particle_P)
      min_particle_P=np.array(min_particle_P)
      ini_particle_x=np.array(ini_particle_x)
      average_particle_x=np.array(average_particle_x)
      max_particle_x=np.array(max_particle_x)
      min_particle_x=np.array(min_particle_x)
      ini_particle_y=np.array(ini_particle_y)
      average_particle_y=np.array(average_particle_y)
      max_particle_y=np.array(max_particle_y)
      min_particle_y=np.array(min_particle_y)
      ini_particle_vx=np.array(ini_particle_vx)
      average_particle_vx=np.array(average_particle_vx)
      max_particle_vx=np.array(max_particle_vx)
      min_particle_vx=np.array(min_particle_vx)
      ini_particle_vy=np.array(ini_particle_vy)
      average_particle_vy=np.array(average_particle_vy)
      max_particle_vy=np.array(max_particle_vy)
      min_particle_vy=np.array(min_particle_vy)

      average_press_in_x=0
      jj=0
      j=0

      filename_time=str(ini_particle_P)
      plt.plot(t,average_particle_y-ini_particle_y,label=r'$P_{0}=$'+filename_time[1:7]+'bar')#
              #,label=fluxlabel[j],linestyle=linetype[j], linewidth=linethickness[j], color=linec[j])
      plt.fill_between(t,max_particle_y-ini_particle_y,min_particle_y-ini_particle_y,alpha=0.5)
      plt.legend(loc=1)#,prop={'size':fontsize1})
      plt.ylabel(r'$\Delta h (=h-h_{0})$[km]')
      plt.xlabel(r'$t\/[\mathrm{s}]$')
      text_num=0
      while (text_num<=3):
        if textbox[text_num]==0:
          plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
        elif textbox[text_num]==1:
          plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
        text_num=text_num+1
      plt.title(r'particle height vs $t$')
      plt.savefig("figure_particleheight_t.png")
  
    print("particle height plot made")
    plt.close()


  if plot_type_index==10 and Dimension==2:
    print  "1 : turbulent heat flux - P plot"
    average_rhovy_in_x=0.0
    average_rhovx_in_x=0.0
    average_press_in_x=0.0
    average_KE_in_x=0.0
    average_rho2_in_x=0.0
    jj=0
    j=0

    for i in profilefilelist:
      ds=load(i)
      time=float(ds.current_time)
      filename_time=str(time)
      l=0
 
      if time>=float(jj)*time_difference_T_P and time<=max_time:
        jj=jj+1
        my_ray_x = ds.ortho_ray(0,(0,0))
        srt_x=np.argsort(my_ray_x['x'])
        my_ray_y = ds.ortho_ray(1,(0,0))
        srt_y=np.argsort(my_ray_y['y'])

        print('Adding lines at t='+filename_time+'s for heigt-rho plot')
        if Dimension==1:
          x_coord=np.array(my_ray_x['x'][srt_x])
          vy=np.array(my_ray_x['y_velocity'][srt_x])
          press=np.array(my_ray_x['pressure'][srt_x])
          average_vy_in_x=vy
          average_press_in_x=press
        if Dimension==2:
          x_coord=np.array(my_ray_y['y'][srt_y])
          

          while (l<number_slice_x):
            position_x=x_range[0]+float(l)*(x_range[1]-x_range[0])/float(number_slice_x)
            my_ray_y = ds.ortho_ray(1,(0,position_x))
            srt_y=np.argsort(my_ray_y['y'])
            vy=np.array(my_ray_y['y_velocity'][srt_y])
            vx=np.array(my_ray_y['x_velocity'][srt_y])
            press=np.array(my_ray_y['pressure'][srt_y])
            rho=np.array(my_ray_y['density'][srt_y])
    
            rhovx=np.multiply(rho,vx)
            rhovy=np.multiply(rho,vy)
            average_rhovy_in_x=(average_rhovy_in_x*float(l)+rhovy)/float(l+1)
            average_rhovx_in_x=(average_rhovx_in_x*float(l)+rhovx)/float(l+1)
            average_press_in_x=(average_press_in_x*float(l)+press)/float(l+1)
            average_rho_in_x=(average_rho2_in_x*float(l)+rho)/float(l+1)

            l=l+1
          l=0
          while (l<number_slice_x):
            position_x=x_range[0]+float(l)*(x_range[1]-x_range[0])/float(number_slice_x)
            my_ray_y = ds.ortho_ray(1,(0,position_x))
            srt_y=np.argsort(my_ray_y['y'])
            vy=np.array(my_ray_y['y_velocity'][srt_y])
            vx=np.array(my_ray_y['x_velocity'][srt_y])
            rho=np.array(my_ray_y['density'][srt_y])
            v2=np.multiply(vx,vx)+np.multiply(vy,vy)
            rho2=np.multiply(rho,rho)
            rho2v2=np.multiply(rho2,v2)
            
            KE=0.5*(rho2v2-2.0*np.multiply(rho,(np.multiply(vx,average_rhovx_in_x)+np.multiply(vy,average_rhovy_in_x)))+(np.multiply(average_rhovx_in_x,average_rhovx_in_x)+np.multiply(average_rhovy_in_x,average_rhovy_in_x)))/np.multiply(average_rho_in_x,average_rho_in_x)
            print KE
            if l==0:
              max_KE_in_x=KE
              min_KE_in_x=KE

            else:
              max_KE_in_x=np.maximum(max_KE_in_x,KE)
              min_KE_in_x=np.minimum(min_KE_in_x,KE)


            average_KE_in_x=(average_KE_in_x*float(l)+KE)/float(l+1)
            l=l+1
          l=0

        if j==0 :
          vy0=vy
          press0=press

        plt.plot(average_press_in_x/1e6,average_KE_in_x,label='t='+filename_time[:6]+'s')#
              #,label=fluxlabel[j],linestyle=linetype[j], linewidth=linethickness[j], color=linec[j])
        plt.fill_between(average_press_in_x/1e6,max_KE_in_x,min_KE_in_x,alpha=0.5)


        plt.legend(loc=1)#,prop={'size':fontsize1})
        plt.ylabel(r'$\tilde{k}\/[\mathrm{cm}^{2}\mathrm{s}^{-2}]$')
        plt.xlabel(r'$P\/[\mathrm{bar}]$')
        plt.xlim(0.001,P_range[1]*1.1/1e6)
        # plt.ylim(2800,T_range[1]*1.1)
        plt.yscale('symlog')
        plt.xscale('linear')
        text_num=0
        while (text_num<=3):
          if textbox[text_num]==0:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13)
          elif textbox[text_num]==1:
            plt.figtext(textposition_LT[0],textposition_LT[1]-(0.05*(text_num-1)+0.05),modelname[text_num],color=color[text_num], fontsize=13,bbox=bbox_props[text_num])
          text_num=text_num+1
        plt.title(r'average turbulent kinetic energy $\tilde{k}$ vs Pressure $P$')
        plt.savefig("figure_turbulentkineticE_P.png")
      j=j+1

    print("turbulentkineticE-P plot made")
    plt.close()



