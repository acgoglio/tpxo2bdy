#
#
# by AC Goglio November 2019
#
#source activate mappyenv # activate my virtual environment
#
# imports
import os # File and directory handling
import time # for waiting n seconds
#
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as NC
from scipy.optimize import curve_fit
from scipy import stats
import collections
import pandas as pd
import csv
import math
import datetime
from operator import itemgetter 
#
from numpy import ma
from scipy.stats import ks_2samp # Kolmogorov-smirnof test
from statsmodels.distributions.empirical_distribution import ECDF # empirical distribution functions
from mpl_toolkits.basemap import Basemap # For plotting maps
from matplotlib.colors import LogNorm # For log scale in pcolor 
#
#################################################################
# The user should modify the following lines to set his run
#################################################################
# General run paraeters:
# General run parameters:
#---------------------
# work dir path (WARNING: file in the directory will be removed..)
workdir_path = '/work/oda/ag15419/tmp/tidal_bdy'+'/'
#bathy_file='/work/oda/ag15419/PHYSW24_DATA/TIDES/DATA0/bathy_meter.nc' # tidal bathimetry
#
#---------------------
# field(s) to be plotted
#-----------------------
mod_depth_var='deptht'
#
# INPUTS
# Path and name of inputs datasets
#
# MODEL DATASETS
# TPXO8
tpxo_path='/users_home/oda/ag15419/tpxo2bdy/'
tpxo_fileprename='mfs_bdytide_'
tpxo_label='TPXO9' 
print ('TPXO9 file templates: ',tpxo_path,tpxo_fileprename,'%tidal_comp%','_grid_','%grid%','.nc')
# FES2014
fes_path='/work/oda/ag15419/PHYSW24_DATA/TIDES/FES2014/'
fes_fileprename='mfs_bdytide_'
fes_label='FES2014'
print ('Fes2014 file templates: ',fes_path,fes_fileprename,'%tidal_comp%','_grid_','%grid%','.nc')
#
#-----------------------
# how to compute field(s) 
#-----------------------
plot_type='abs'
# Options are:
# abs  --> absolute values of Re and Im parts
# sgn  --> values (original sign) of Re and Im parts

########################################################
# DO NOT CHANGE THE CODE BELOW THIS LINES
########################################################

# Move to the work dir and clean it if it exists otherwise create it
if os.path.exists(workdir_path) : 
   os.chdir(workdir_path) # mv to work dir
   print('I am moving in the directory: ',os.getcwd()) # pwd
   os.listdir() # ls
   # We do not want to remove input files if externally produced!!!!
   #print('WARNING: I am going to clean this directory in a while..')
   #time.sleep(10) # sleep for 10 seconds before removig everything in the work dir!
   #file2berm=os.path.join(workdir_path,'*.*') 
   #os.rm(file2berm) # clean dir
   #print ('Cleaning:',file2berm)
   #print('I am in the clean directory: ',os.getcwd()) # pwd
else:
  os.mkdir(workdir_path)
  os.chdir(workdir_path)
  print('I am in a new work directory: ',os.getcwd()) 
  

##################################
#
# Read the bathymetry
#
#bathy_field = NC.Dataset(bathy_file,'r')
#vals_bathy=bathy_field.variables['Bathymetry'][:]
    
# ===============================
# Loop on U/V grids:
# ===============================
for grid in ('U','V','T'):
  if grid=='U':
     field_prename='u'
     field_units='m/s'
     fielddiff_units='m/s'
     #field_vallim=0.04
     #fielddiff_vallim=1 
  elif grid=='V':
     field_prename='v'
     field_units='m/s'
     fielddiff_units='m/s'
     #field_vallim=0.04
     #fielddiff_vallim=1 
  elif grid=='T':
     field_prename='z'
     field_units='m'
     fielddiff_units='m'
     #field_vallim=0.8 
     #fielddiff_vallim=1 
  print ('Working on grid: ',grid)
  # loop on Re Im parts
  for ReIm_idx in (1,2,3,4,5):
   if ReIm_idx==1:
      part='Re'
      field_name=field_prename+str(ReIm_idx)
   elif ReIm_idx==2:
      part='Im'
      field_name=field_prename+str(ReIm_idx)
   elif ReIm_idx==3:
      part='Amp_'+field_prename
      field_name1=field_prename+'1'
      field_name2=field_prename+'2'
      field_name=part
   elif ReIm_idx==4:
      part='Re_comp'
      field_name1=field_prename+'1'
      field_name2=field_prename+'2'
      field_name=field_prename+field_prename+'1'+'_comp'
   elif ReIm_idx==5:
      part='Im_comp'
      field_name1=field_prename+'1'
      field_name2=field_prename+'2'
      field_name=field_prename+field_prename+'2'+'_comp'
   print ('Working on ',part,' part')
   print ('Working on field: ',field_name)
   # Loop on Atlantic Box sides
   for side in ('North','West','South'): 
    print ('Working on Atlantic Box side ',side)
    if side=='North':
      idx_side_min=0
      idx_side_max=320 
      var2plot_xax='lons'
      coo_side_min=-18.200
      coo_side_max=-4.000
      x_label='Longitude'
    elif side=='South':
      idx_side_min=1320 
      idx_side_max=1520 
      var2plot_xax='lons'
      coo_side_min=-18.100
      coo_side_max=-8.000
      x_label='Longitude'
    elif side=='West':
      idx_side_min=2620 
      idx_side_max=2950 
      var2plot_xax='lats'
      coo_side_min=30.000
      coo_side_max=46.000
      x_label='Latitude'
    print ('xranges ',idx_side_min,idx_side_max)
    # Loop on tidal comp
    # Plots Per grid and per part
    plotname=workdir_path+'tpxoVsfes_'+field_name+'_'+side+'.jpg'
    print ('Plot path/name: ',plotname)
    plt.figure(figsize=(18,12))
    plt.rc('font', size=11)
    plt.subplot(2,2,1)
    # Plot 1
    plt.title ('Comparison TPXO9/FES2014 along the '+side+' side of the Atlantic bdy --- '+part+' --- '+field_name+' field')
    #
    for comp in ('M2','S2','K1','O1','N2','K2','P1','Q1'):
        print ('Open tidal comp: ',comp)
        # Build the path/name of the nc file and open it 
        tpxo2open=tpxo_path+tpxo_fileprename+comp+'_grid_'+grid+'.nc'
        fes2open=fes_path+fes_fileprename+comp+'_grid_'+grid+'.nc'
        print ('Input files = ',tpxo2open,fes2open)
        modtpxo=NC.Dataset(tpxo2open,'r')
        modfes=NC.Dataset(fes2open,'r')
        #
        print('Open fields..')
        #xb=modtpxo.variables['xb'+grid][:]
        if ReIm_idx==1 or ReIm_idx==2: 
           fieldtpxo=modtpxo.variables[field_name][:]
           lons = np.squeeze(modtpxo.variables['nav_lon'][:]) 
           lats = np.squeeze(modtpxo.variables['nav_lat'][:]) 
           fieldfes=modfes.variables[field_name][:]
        elif ReIm_idx==3:
           fieldtpxo1=modtpxo.variables[field_name1][:]
           fieldtpxo2=modtpxo.variables[field_name2][:]
           fieldtpxo=np.sqrt(np.power(fieldtpxo1,2)+np.power(fieldtpxo2,2))
           lons = np.squeeze(modtpxo.variables['nav_lon'][:]) 
           lats = np.squeeze(modtpxo.variables['nav_lat'][:]) 
           fieldfes1=modfes.variables[field_name1][:]
           fieldfes2=modfes.variables[field_name2][:]
           fieldfes=np.sqrt(np.power(fieldfes1,2)+np.power(fieldfes2,2))
        elif ReIm_idx==4 or ReIm_idx==5:
           fieldtpxo1=modtpxo.variables[field_name1][:]
           fieldtpxo2=modtpxo.variables[field_name2][:]
           amptpxo=np.sqrt(np.power(fieldtpxo1,2)+np.power(fieldtpxo2,2))
           phatpxo=np.arctan2(fieldtpxo2,fieldtpxo1)
           np.where(phatpxo<0,phatpxo+360,phatpxo)
           np.where(phatpxo>=360,phatpxo-360,phatpxo)
           fieldtpxo1=amptpxo*np.cos(phatpxo)
           fieldtpxo2=amptpxo*np.sin(phatpxo)
           lons = np.squeeze(modtpxo.variables['nav_lon'][:])
           lats = np.squeeze(modtpxo.variables['nav_lat'][:])
           fieldfes1=modfes.variables[field_name1][:]
           fieldfes2=modfes.variables[field_name2][:]
           if ReIm_idx==4:
              fieldtpxo=fieldtpxo1
              fieldfes=fieldfes1
           elif ReIm_idx==5:
              fieldtpxo=fieldtpxo2
              fieldfes=fieldfes2

        x_vals=globals()[var2plot_xax][idx_side_min:idx_side_max]
        if plot_type=='abs':
           y_vals=np.abs(np.squeeze(fieldtpxo))[idx_side_min:idx_side_max]
           y_vals2=np.abs(np.squeeze(fieldfes))[idx_side_min:idx_side_max]
        elif plot_type=='sgn':
           y_vals=np.squeeze(fieldtpxo)[idx_side_min:idx_side_max]
           y_vals2=np.squeeze(fieldfes)[idx_side_min:idx_side_max]

        plt.plot(x_vals,y_vals,label = 'TPXO9 '+comp)
        plt.plot(x_vals,y_vals2,'--',color=plt.gca().lines[-1].get_color(), label = 'FES14 ' +comp)
        print('Done..', comp)
 
    

    print ('Grid and legend..')
    plt.grid ()
    plt.xlim(coo_side_min,coo_side_max)
    #plt.ylim(0,field_vallim)
    #plt.axhline(y=0 color='black')
    plt.ylabel (field_name+' '+field_units)
    plt.xlabel (x_label)
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    # Plot 2
    plt.subplot(2,2,3)
    plt.title ('Diff TPXO9-FES2014 along the '+side+' side of the Atlantic bdy --- '+part+' part --- '+field_name+' field')
    #
    for comp in ('M2','S2','K1','O1','N2','K2','P1','Q1'):
        print ('Open tidal comp: ',comp)
        # Build the path/name of the nc file and open it 
        tpxo2open=tpxo_path+tpxo_fileprename+comp+'_grid_'+grid+'.nc'
        fes2open=fes_path+fes_fileprename+comp+'_grid_'+grid+'.nc'
        print ('Input files = ',tpxo2open,fes2open)
        modtpxo=NC.Dataset(tpxo2open,'r')
        modfes=NC.Dataset(fes2open,'r')
        #
        print('Open fields and compute diffs..')
        if ReIm_idx==1 or ReIm_idx==2:
           fieldtpxo=modtpxo.variables[field_name][:]
           fieldfes=modfes.variables[field_name][:]
        elif ReIm_idx==3:
           fieldtpxo1=modtpxo.variables[field_name1][:]
           fieldtpxo2=modtpxo.variables[field_name2][:]
           fieldtpxo=np.sqrt(np.power(fieldtpxo1,2)+np.power(fieldtpxo2,2))
           fieldfes1=modfes.variables[field_name1][:]
           fieldfes2=modfes.variables[field_name2][:]
           fieldfes=np.sqrt(np.power(fieldfes1,2)+np.power(fieldfes2,2))
        elif ReIm_idx==4 or ReIm_idx==5:
           fieldtpxo1=modtpxo.variables[field_name1][:]
           fieldtpxo2=modtpxo.variables[field_name2][:]
           amptpxo=np.sqrt(np.power(fieldtpxo1,2)+np.power(fieldtpxo2,2))
           phatpxo=np.arctan2(fieldtpxo2,fieldtpxo1)
           np.where(phatpxo<0,phatpxo+360,phatpxo)
           np.where(phatpxo>=360,phatpxo-360,phatpxo)
           fieldtpxo1=amptpxo*np.cos(phatpxo)
           fieldtpxo2=amptpxo*np.sin(phatpxo)
           fieldfes1=modfes.variables[field_name1][:]
           fieldfes2=modfes.variables[field_name2][:]
           if ReIm_idx==4:
              fieldtpxo=fieldtpxo1
              fieldfes=fieldfes1
           elif ReIm_idx==5:
              fieldtpxo=fieldtpxo2
              fieldfes=fieldfes2

        x_vals=globals()[var2plot_xax][idx_side_min:idx_side_max]
        if plot_type=='abs':
           y_vals=(np.abs(np.squeeze(fieldtpxo))-np.abs(np.squeeze(fieldfes)))[idx_side_min:idx_side_max]
        elif plot_type=='sgn':
           y_vals=(np.squeeze(fieldtpxo)-np.squeeze(fieldfes))[idx_side_min:idx_side_max]

        plt.plot(x_vals,y_vals,label = 'TPXO9-FES2014 '+comp)
        print('Done..', comp)
    print ('Grid and legend..')
    plt.grid ()
    plt.xlim(coo_side_min,coo_side_max)
    #plt.axhline(y=0 color='black')
    plt.ylabel (field_name+' diff ['+fielddiff_units+']')
    plt.xlabel (x_label)
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    # Save and close 
    plt.savefig(plotname)
    plt.clf()

