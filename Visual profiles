
#Duarte Oom#
# Code to visualize the plots temporal profiles..dont include the criteria, MIn,moving window and new tests (T-test, Kde and MAnn Whitney)

try:
  import sys
except ImportError:
  print 'Sys is not installed.'
  exit(-1)

try:
  import time
except ImportError:
  print 'time is not installed.'
  exit(-1)
  
try:
  import GTr1ReadFunc_vs_3_1_wings
except ImportError:
  print 'bacafunc is not installed.'
  exit(-1)
  
try:
  import os
except ImportError:
  print 'os is not installed'
  exit(-1)  
 
try:
  import GTr1AlgFunc_vs_3_3_wings_vs2
except ImportError:
  print 'bacafunc is not installed.'
  exit(-1) 
  
try:
  import dba
except ImportError:
  print 'dba is not installed.'
  exit(-1)   
    
 
try:
  import scipy.io
except ImportError:
  print 'scipy.io is not installed'
  exit(-1)
  
  
try:
  import pandas as pnds
except ImportError:
  print 'scipy.io is not installed'
  exit(-1)  
   
 
 
try:
  import datetime as DT
except ImportError:
  print 'datetime is not installed.'
  exit(-1)
  
  
  
try:
  from matplotlib.dates import date2num
except ImportError:
  print 'Matplotlib.dates.date2num is not installed.'
  exit(-1)

try:
  from matplotlib.dates import num2date

except ImportError:
  print 'Matplotlib.dates.num2date is not installed.'
  exit(-1) 
 
from scipy import stats  #for t-test

from scipy.interpolate import interp1d

File1= sys.argv[1]  #config file
   
ficheiro_leitura=File1
  
  
#i=int(i)
#j=int(j)

#Create the pixel profile file
#WorkDir = os.path.dirname(File2)
#FirefileName = os.path.basename(File2)[0:-4] + '_' + str(i) + '_' + str(j) + '_profile.txt'
#f2 = open(WorkDir + '/' + FirefileName, 'w')




class Log(object):
  def write(self, msg):
    print "LOG: %s" % msg
    
    
    
def normalize_data(x):
  import numpy 
  x2=x*0.0
  mean_x=numpy.mean(x)
  std_x=numpy.std(x)
  for z in range(0,len(x)):
    x2[z]=(x[z]-mean_x)/std_x
  
  return x2


def GetStrDoY(DoY):

  if len(str(int(DoY))) == 1:
    strDoY = '00' + str(int(DoY))
  elif len(str(int(DoY))) == 2:
    strDoY = '0' + str(int(DoY))
  else:
    strDoY = str(int(DoY))

  return strDoY

import scipy.stats as stat

def bootstrap(sample, samplesize, nsamples, statfunc):
    """
    Arguments:
       sample - input sample of values
       nsamples - number of samples to generate
       samplesize - sample size of each generated sample
       statfunc- statistical function to apply to each generated sample.
 example - bootstrap(sample, samplesize = None, nsamples = 1000, statfunc = mean):
    Performs resampling from sample with replacement, gathers
    statistic in a list computed by statfunc on the each generated sample.
    """
    if samplesize is None:                                                                   
        samplesize=len(sample)
    #print "input sample = ",  sample
    n = len(sample)
    X = []
    #print "range",range(nsamples)
    for i in range(nsamples):
        
        #print "i = ",  i, 
        resample = [sample[j] for j in stat.randint.rvs(0, n-1, size=samplesize)] 
        x = statfunc(resample)
        X.append(x)
    return X




def calCPD_GT(plot_nb,data_nir,numleituras, diajul, cellref, dist,GregDay, YearJulianDay,ano,year,indices_core_year,Greg_day_core_year,i,j):
	import numpy
	import matplotlib.pyplot as plt
	from scipy.stats import t
	import rpy2.robjects as robjects
	from rpy2.robjects.packages import importr
	import GTr1AlgFunc_vs_3_3_wings_vs2
	import dba
	
	import GTr1ReadFunc_vs_3_1_wings
	#reload(GTr1AlgFunc_vs_3_1)
	import GTr1Par_vs_3_1
	changepoint = importr("changepoint")
	Vector = robjects.FloatVector # object to convert python vector into R vector
	mmc=changepoint.multiple_mean_norm
	#mmc=changepoint.multiple_var_norm
	
	factor=2.0
		
	
	ficheiros, ano, mes, dia, numlinhas, numcolunas, numleituras, caminho_Mask,GregDay,YearJulianDay,year,indices_core_year,Greg_day_year_core = GTr1ReadFunc_vs_3_1_wings.loadDataspecs(str(ficheiro_leitura))
	GregDay=numpy.reshape(GregDay, len(GregDay),order='F')
	validos, val_dias, val_nir, numval,validos_year_core,val_dias_year_core,val_nir_year_core,numval_year_core=GTr1AlgFunc_vs_3_3_wings_vs2.validar(data_nir[i-1,j-1,:], numleituras, diajul,GregDay,ano,year,indices_core_year,Greg_day_year_core)
	
	
	if numval_year_core>GTr1Par_vs_3_1.numcasos:
	  print ('To process') 
	else:
	  print ('Not enough valid data') 
	  sys.exit("Not enough valid data")	    	  
	
	#convert to julian date
	val_dias_jul = numpy.zeros(val_dias.shape,numpy.int32)
	for jj in range(0,val_dias.shape[0]):	  
	  val_dias_jul[jj]=int(DT.datetime.strftime( num2date(val_dias[jj]), "%j" ))	
	
	
	  

	

	  	  


	###################CALCULATE THE MOVING WINDOW TO THE MINIMUM SERIES#######################
		
	####MINIMUM############Filter CO(min-max-max-min) ############################
	
	novoy=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(val_nir,2,'min')
	novoy1=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy,2,'max')
	novoy2=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy1,2,'max')
	novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy2,2,'min')
	
	novoy3_1=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3,3,'min')
	novoy3_1=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3_1,3,'max')
	novoy3_1=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3_1,3,'max')
	novoy3_1=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3_1,3,'min')
	
	
	novoy3_2=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3_1,4,'min')
	novoy3_2=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3_2,4,'max')
	novoy3_2=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3_2,4,'max')
	novoy3_2=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3_2,4,'min')
	
	
	novoy3_3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3_2,5,'min')
	novoy3_3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3_3,5,'max')
	novoy3_3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3_3,5,'max')
	novoy3_3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3_3,5,'min')	
	
	
	#novoy=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(val_nir,2,'min')
	#novoy1=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy,2,'max')
	#novoy2=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy1,2,'max')
	#novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy2,2,'min')
	
	#novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3,3,'min')
	#novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3,3,'max')
	#novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3,3,'max')
	#novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3,3,'min')
	
	
	#novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3,4,'min')
	#novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3,4,'max')
	#novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3,4,'max')
	#novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3,4,'min')
	
	
	#novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3,5,'min')
	#novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3,5,'max')
	#novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3,5,'max')
	#novoy3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left_morphological(novoy3,5,'min')	
	
	
	
	
	
	
	#novoy,total=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(val_nir,3,0,0)#Moving window for time series (minimum).
	#novoy1,total=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy,3,2,0)#Moving window for maxima time series (maximum)--see pseudo code g)
	#novoy2,total=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy1,3,2,0)#Moving window for maxima time series (maximum)--see pseudo code g)
	#novoy3,total=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy2,3,0,0)#Moving window for maxima time series (minimum)--see pseudo code g)
	
	#novoy3,total3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy3,5,0,0)#Moving window for time series (minimum)
	#novoy3,total3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy3,5,4,0)#Moving window for maxima time series (maximum)--see pseudo code g)
	#novoy3,total3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy3,5,4,0)#Moving window for maxima time series (maximum)--see pseudo code g)
	#novoy3,total3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy3,5,0,0)#Moving window for maxima time series (minimum)--see pseudo code g)	
	
	
	#novoy3,total3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy3,7,0,0)#Moving window for time series (minimum)
	#novoy3,total3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy3,7,6,0)#Moving window for maxima time series (maximum)--see pseudo code g)
	#novoy3,total3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy3,7,6,0)#Moving window for maxima time series (maximum)--see pseudo code g)
	#novoy3,total3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy3,7,0,0)#Moving window for maxima time series (minimum)--see pseudo code g)	
	
	
	#create the nan's vector
	#first task - convert values above 25 as nan - create a new vector
	
	novoy4=numpy.array(novoy3_3)
	#novoy4=numpy.array(novoy3)
	
	

	index_invalid=numpy.where(novoy4>30) #in the core year only the valids
	if (index_invalid[0].shape[0]>0): #there are values above 25
	  
	  novoy4[index_invalid]=numpy.NAN #replace with nan
	else: 
	  #novoy4=novoy3 #STAY THE SAME
	  novoy4=novoy3_3

	year_number=int(year)#check the year
	year_previous_number=str(year_number-1) #the previous year for wings
	year_after_number=str(year_number+1) #the following year for wings
	
	
	Day_one_wings=(year_previous_number+'12'+'01')
	Day_last_wings=(year_after_number+'01'+'31')

	greg_date_all_start= date2num(DT.datetime.strptime(Day_one_wings, "%Y%m%d"))
	greg_date_all_stop= date2num(DT.datetime.strptime(Day_last_wings, "%Y%m%d"))

	
	#date array
	GregDay_all=numpy.arange(greg_date_all_start,greg_date_all_stop+1,1) #see numpy.arange


	#convert to julian date
	GregDay_all_jul = numpy.zeros(GregDay_all.shape,numpy.int32)
	for jj in range(0,GregDay_all.shape[0]):	  
	  GregDay_all_jul[jj]=int(DT.datetime.strftime( num2date(GregDay_all[jj]), "%j" ))	
	#YearJulianDay.append(int(DT.datetime.strftime( num2date(greg_date), "%Y%j" )))

	#Core year with nan's
	Day_one_core=(str(year_number)+'01'+'01')
	Day_last_core=(str(year_number)+'12'+'31')	
	greg_date_core_start= date2num(DT.datetime.strptime(Day_one_core, "%Y%m%d"))
	greg_date_core_stop= date2num(DT.datetime.strptime(Day_last_core, "%Y%m%d"))
	Greg_day_core_year_nan=numpy.arange(greg_date_core_start,greg_date_core_stop+1,1) #see numpy.arange

	#create nan's array
	val_nir_nan=numpy.empty(GregDay_all.shape[0])
	val_nir_nan[:]=numpy.NAN


	
	#replace the nan's vector with observations from val_dias
	for jj in range (0,val_dias.shape[0]):
	  day=val_dias[jj]
	  index_date=numpy.where(GregDay_all==day)
	  val_nir_nan[index_date]=novoy4[jj]		
	



	#index for nan's and not nan's localization for plotting
	index_nan=numpy.isnan(val_nir_nan)
	
	index_nan_val_nir_nan=numpy.where(index_nan==True)
	index_not_nan_val_nir_nan=numpy.where(index_nan==False)
		
	###############MAXIMUM######################
	#novoy7,total3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(val_nir,5,4,0)#Moving window for  time series (maximum)
	#novoy9,total3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy7,5,0,0)#Moving window for maxima time series (minimum)--see pseudo code g)
	#novoy9,total3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy9,5,0,0)#Moving window for maxima time series (minimum)--see pseudo code g)
	#novoy9,total3=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window(novoy9,5,4,0)#Moving window for maxima time series (maximum)--see pseudo code g)
	
		
	#########################parameters##########
	
	density_alldata=float(val_dias_year_core.shape[0])/(val_dias_year_core[-1]-val_dias_year_core[0]) #time series density in the core year-number of valids divided by number of days
	
	'''
	#number of valids in the window will depend of the density of the time series
	if(density_alldata<=0.25):
	  num_valids_inside_window=4
	elif ((density_alldata>0.25) &(density_alldata<=0.333)):
	  num_valids_inside_window=5
	elif ((density_alldata>0.333) &(density_alldata<=0.416)):
	  num_valids_inside_window=6

	elif ((density_alldata>0.416) &(density_alldata<=0.5)):
	  num_valids_inside_window=7

	elif (density_alldata>0.5):
	  num_valids_inside_window=8	
	'''
	

	num_valids_inside_window=6
	#density_window=0.25


	window_size=16
	coef=2
	
	group=1 #0 with groups values from NB
	power_refl=4 #power in calculating score segments=(pre_obs_vector_mean*num_days)/pow(refl_vector_mean,power_refl)
	Num_classes_nb=3 #number of classes of natural breaks 
	#percent_parameter=10
	
	##########################################################


	#Calculate the median of the left window for the minimum time series.Predict the next value
	
	#novoy_median,total,window=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left(novoy3,10,'median',2)
	
	novoy_median,total,window=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left(val_nir_nan,window_size,num_valids_inside_window,'nanmedian',coef)
	#print"pre",novoy_median[107]
	#print"val",val_nir_nan[107]
	
	
		
		
	####################check how many novox elements are in the core year. Create the indexes###############
	#common_days_novox_year_core_index = numpy.in1d(val_dias,Greg_day_core_year)
	common_days_novox_year_core_index = numpy.in1d(GregDay_all,Greg_day_core_year_nan)

	#common_days_novox_year_core = val_dias[common_days_novox_year_core_index]
	common_days_novox_year_core = GregDay_all[common_days_novox_year_core_index]

	#common_days_novoy_year_core = novoy3[common_days_novox_year_core_index]
	common_days_novoy_year_core = val_nir_nan[common_days_novox_year_core_index]	# from the minimum time series
	#########################################################################################################

	
	

	
	
	###############################reflectances original differences##########################################
		
	#the calculations will be based on the differences only when we have data (so taking into account the difference in days between observations different from nan)
	#reflectances original min/max/max/min
	#diff_reflect2=numpy.zeros((novoy4.shape[0]), numpy.float)
	diff_reflect2=numpy.empty(novoy4.shape[0])
	diff_reflect2[:]=numpy.NAN	
	for z in range(0, len(novoy4)-1):
	  diff_reflect2[z+1]=(novoy4[z+1]-novoy4[z])/(val_dias[z+1]-val_dias[z])	  
	

	
	diff_reflect2_nan=numpy.empty(GregDay_all.shape[0])
	diff_reflect2_nan[:]=numpy.NAN      
	#for z in range(0, len(val_nir_nan)-1):
	  #diff_reflect2_nan[z+1]=(val_nir_nan[z+1]-val_nir_nan[z])/(GregDay_all[z+1]-GregDay_all[z])	  
  
	#replace the nan's vector with observations from val_dias
	for jj in range (0,val_dias.shape[0]):
	  day=val_dias[jj]
	  index_date=numpy.where(GregDay_all==day)
	  diff_reflect2_nan[index_date]=diff_reflect2[jj]	  
	#taking the core year
	  
	diff_reflected2_core_year=diff_reflect2_nan[common_days_novox_year_core_index]
	
	
	
	################calculate the modified trimmed median for the daily differences time series#####################
      
	
	diff_reflect2_predicted,diff_reflect2_trimmed,window=GTr1AlgFunc_vs_3_3_wings_vs2.moving_window_left(diff_reflect2_nan,window_size,num_valids_inside_window,'nanmedian',coef)
	#print"tri",diff_reflect2_trimmed
	
	######################################Final predicted###########################################################
	#final_predicted=numpy.zeros((novoy3.shape[0]), numpy.float)
	#for z in range (0,novoy3.shape[0]):
	  #final_predicted[z]=novoy_median[z]+ ((window/2)*diff_reflect2_trimmed[z])

	final_predicted=numpy.zeros((val_nir_nan.shape[0]), numpy.float)
	for z in range (0,val_nir_nan.shape[0]):
	  final_predicted[z]=novoy_median[z]+ ((window/2)*diff_reflect2_trimmed[z])	  
	  
	#print"sum",final_predicted 
	#print "ori",novoy8
	
	#taking the core year
	
	final_predicted_core_year=final_predicted[common_days_novox_year_core_index]
	
	#index for nan's localization for plotting
	index_nan=numpy.isnan(final_predicted)
		
	index_nan_final_predicted_nan=numpy.where(index_nan==True)	
	index_not_nan_final_predicted_nan=numpy.where(index_nan==False)	
	
	
	#####################predicted-observed################################################
	
	#diff_pre_obs=final_predicted-novoy3
	diff_pre_obs=final_predicted-val_nir_nan
	#print "dif",diff_pre_obs
	
	#taking the core year
		
	diff_pre_obs_core_year=diff_pre_obs[common_days_novox_year_core_index]	
	
	#index for nan's localization for plotting
	index_nan=numpy.isnan(diff_pre_obs)
			
	index_nan_diff_pre_obs_nan=numpy.where(index_nan==True)		
	index_not_nan_diff_pre_obs_nan=numpy.where(index_nan==False)
	
	

###########################Finite values#######################################################
	index_predicted=numpy.isfinite(diff_pre_obs_core_year)
	index_predicted_index=numpy.where(index_predicted==True)
	val_nir_nan_finite=val_nir_nan[index_predicted_index]
	
	#rearrange new variables
	common_days_novoy_year_core=common_days_novoy_year_core[index_predicted_index] #create new vector based on finite values of the difference between predicted and observed
	diff_reflected2_core_year=diff_reflected2_core_year[index_predicted_index]
	diff_pre_obs_core_year=diff_pre_obs_core_year[index_predicted_index]
	common_days_novox_year_core=common_days_novox_year_core[index_predicted_index]
	final_predicted_core_year=final_predicted_core_year[index_predicted_index]
#################################################################################################################################


####################calculate the values below certain percentage (parameter)and the index###############
	###################################NATURAL BREAKS#############################################################
	common_days_novoy_year_core_jenks=numpy.array(common_days_novoy_year_core)
	kclass=GTr1AlgFunc_vs_3_3_wings_vs2.getJenksBreaks(common_days_novoy_year_core_jenks,Num_classes_nb) #makes the sort
	cut_value=kclass[1]

	
	group_first=numpy.where(common_days_novoy_year_core<=kclass[1]) #first quartil
	group_second=numpy.where(common_days_novoy_year_core<=kclass[2])	
	first_nb=group_first[0].shape[0]/float(common_days_novoy_year_core.shape[0])
	second_nb=group_second[0].shape[0]/float(common_days_novoy_year_core.shape[0])
	

	##########################################plot###################################################################
	
	values,base=numpy.histogram(common_days_novoy_year_core,bins=numpy.arange(min(numpy.rint(common_days_novoy_year_core)), max(numpy.rint(common_days_novoy_year_core)) + 2, 1))
	values_normed=values/float(common_days_novoy_year_core.shape[0])
	if (plot_nb==1):
	  plt.figure()
	  fig=plt.gcf()
	  fig.set_size_inches(15.5,9.5)
  
	  ax1=plt.subplot(1,1,1)
	    
	  
	  histogram=plt.hist(common_days_novoy_year_core, bins=numpy.arange(min(numpy.rint(common_days_novoy_year_core)), max(numpy.rint(common_days_novoy_year_core)) + 2, 1), normed=False) #
	  ax1.set_xticks(histogram[1])
  
	  ax1.legend(loc='upper left',prop={'size':8})	    
	  #ax1.grid(True)
		  
	  ax2 = ax1.twinx()
	  ax2.grid(True)
	  cumulative = numpy.cumsum(values_normed)
	  ax2.plot(base[1:,], cumulative, c='green',label='F(x)')
	  ax2.legend(loc='upper right',prop={'size':8})
	  ax1.set_xlabel('NIR')
		      
	  ax2.set_ylabel('F(x)')
	  ax1.set_ylabel('Frequency')
	  
	  #put lines s of the natural breaks
	  plt.axvline(kclass[1],color='r')
	  plt.axvline(kclass[2],color='r')
	  #plt.axvline(kclass[3],color='r')
	  plt.text(kclass[1], 0.5, '%:' + str(numpy.rint(first_nb*100))+ '/' + str(round(kclass[1],2)),color='r')
	  plt.text(kclass[2], 0.5, '%:' + str(numpy.rint(second_nb*100)),color='r')
	  plt.title('Tile: LL' + str(numlinhas) + '_' +'CC' + str(numcolunas) + '_row: ' +str(i) + ', col:' + str(j) + '; ' + 'N.classes NB:' + '' + str(Num_classes_nb))    
	  
	
	#################################################DEFINE GROUP OF VALUES########################################################################	
	percent_parameter=numpy.rint(first_nb*100)
	
	#the maximum value for natural breaks is 33%
	if percent_parameter>25:
	  percent_parameter=25
	
	
	refl_values,index_refl,value_ref=GTr1AlgFunc_vs_3_3_wings_vs2.cal_percentage(common_days_novoy_year_core,percent_parameter)
	refl_values_group=common_days_novoy_year_core[index_refl]

	#DETERMINE THE GROUP WHERE ARE THE DATA TO BE ANALYZED
	#DAYS
	days_group=common_days_novox_year_core[index_refl]

	
	#DAILY DIFFERENCE GROUP
	daily_diff_group=diff_reflected2_core_year[index_refl] 
	
	#PREDICTED
	
	final_predicted_group=final_predicted_core_year[index_refl] 
	#DIFFERENCE BETWEEN PREDICTED AND OBSERVED
		
	diff_pre_obs_group=diff_pre_obs_core_year[index_refl] #index_refl- index calculated from the "refl_values,index_refl=GTr1AlgFunc_vs_3_3_wings_vs2.cal_percentage(common_days_novoy_year_core,10)"
	
	
	#a)choose the biggest negative value from daily differences and where is it-observed in the group of values
	
	daily_diff_group_max=numpy.nanmin(daily_diff_group) #minimum valuewithout nan!
	daily_diff_group_max_index=daily_diff_group_max==daily_diff_group#index of the maximum (True)

	
	
	#b)choose the highest difference from diff_pre_obs_group and where is it-observed in the group of values
	#put with two decimal - NO EFFECT#############################################################
	diff_pre_obs_group=numpy.round(diff_pre_obs_group,2)
	
	
	'''
	diff_pre_obs_group_max=numpy.nanmax(diff_pre_obs_group) #maximum valuewithout nan!
	 
	diff_pre_obs_group_max_index=diff_pre_obs_group_max==diff_pre_obs_group#index of the maximum (True)
	'''
	
	#########################SCORING##################################################################3
	
	#(a) Choose the values in the group where predicted-observed is positive and daily difference is negative(to define the area)
	diff_pre_obs_group_positive=numpy.where(diff_pre_obs_group>0)
	
	##b) Check where is the first negative----diff_pre_obs_group_negative[1][0]
	#diff_pre_obs_group_negative=numpy.where(diff_pre_obs_group<0)
	
	
	##same size of diff_pre_obs_group_positive
	diff_pre_obs_group_sub_group=diff_pre_obs_group[0][diff_pre_obs_group_positive[1]]
	daily_diff_group_sub_group=daily_diff_group[0][diff_pre_obs_group_positive[1]]
	days_group_subgroup=days_group[0][diff_pre_obs_group_positive[1]]
	refl_values_group_subgroup=refl_values_group[0][diff_pre_obs_group_positive[1]]
	
	predicted_group_subgroup=final_predicted_group[0][diff_pre_obs_group_positive[1]]
	
	
	
	##choose the values in the previous group (a) where predicted-observed is positive (diff_pre_obs_group>0) and daily difference is negative (daily_diff_group<0)
	daily_diff_sub_group_negative=numpy.where(daily_diff_group_sub_group<0) 
	
	if (daily_diff_sub_group_negative[0].shape[0]==0): #empty- none of areas have candidates
	    print "no candidates dates"
	    #sys.exit() #exit the script
	
	##POTENTIAL CANDIDATES TO SCORING-same size of daily_diff_sub_group_negative
	diff_pre_obs_group_sub_group_candidates=diff_pre_obs_group_sub_group[daily_diff_sub_group_negative]
	daily_diff_group_sub_group_candidates=daily_diff_group_sub_group[daily_diff_sub_group_negative]
	days_group_subgroup_candidates=days_group_subgroup[daily_diff_sub_group_negative]
	refl_values_group_subgroup_candidates=refl_values_group_subgroup[daily_diff_sub_group_negative]
	
	predicted_values_group_subgroup_candidates=predicted_group_subgroup[daily_diff_sub_group_negative]
	
	#################################################################################################
	
	area_segment=[]
	score_individual_area=[]	
	
	predicted_all=[]	 
	pre_obs_all=[]		
	refl_all=[]	
	days_all=[]
	daily_diff_all=[]
	
	ll=0
	
	while ll<diff_pre_obs_core_year.shape[0]-1:
	#for ll in range(diff_pre_obs_group.shape[1]):
	  
	  predicted=[]
	  pre_obs=[]
	  refl=[]
	  days=[]
	  daily_diff=[]
	  area_segment=[]
	  score_individual_area=[]		  
	  #date_ini=days_group[0][ll]
	  
	  #loop for each value
	  if (diff_pre_obs_core_year[ll]<=0): #or (days_group[0][jj]-days_group[0][index_new]>window_size):
	    ll=ll+1
	    continue
	  
	  else:
	    
	    #loop for each area
	    # first positive (to append)
	    
	    predicted.append(final_predicted_core_year[ll])
	    pre_obs.append(diff_pre_obs_core_year[ll])
	    refl.append(common_days_novoy_year_core[ll])
	    days.append(common_days_novox_year_core[ll]) 
	    daily_diff.append(diff_reflected2_core_year[ll])
	    
	    ll=ll+1
	    for jj in range(ll,diff_pre_obs_core_year.shape[0]):
	      
	      if (diff_pre_obs_core_year[jj]<=0)or (common_days_novox_year_core[jj]-common_days_novox_year_core[jj-1]>window_size):
		if ll==jj:
		  jj=jj+1
		
		break
	      else:
	      
		
		diffpre_obs_orginal=diff_pre_obs_core_year[jj]
		refl_original=common_days_novoy_year_core[jj]
		days_original=common_days_novox_year_core[jj]
		predicted_original=final_predicted_core_year[jj]
		daily_diff_original=diff_reflected2_core_year[jj]
			
			
		predicted.append(predicted_original)
		pre_obs.append(diffpre_obs_orginal)
		refl.append(refl_original)
		days.append(days_original) 
		daily_diff.append(daily_diff_original)
	    
	    ll=jj
	  
	  predicted_all.append(predicted)  
	  pre_obs_all.append(pre_obs)
	  refl_all.append(refl)
	  days_all.append(days) 
	  daily_diff_all.append(daily_diff)
	  
	    
	    
	    
	predicted_all=numpy.array(predicted_all) #convert list to array
	pre_obs_all=numpy.array(pre_obs_all)#convert list to array
	refl_all=numpy.array(refl_all) #convert list to array
	days_all=numpy.array(days_all)#convert list to array	
	daily_diff_all=numpy.array(daily_diff_all)#convert list to array
	    
	    
	    
	    
	#Score each area (the score of each area is the sum of each segment of each area)####################
	
	
	score_area=numpy.zeros((days_all.shape[0]), numpy.float)
	for uu in range(0,days_all.shape[0]): #for each area
	  days=numpy.array(days_all[uu])
	  pre_obs=numpy.array(pre_obs_all[uu])
	  refl=numpy.array(refl_all[uu])
	  predicted=numpy.array(predicted_all[uu])
	  daily_diff=numpy.array(daily_diff_all[uu])
	  
	  ##############score segments###########################################
	  if (days.shape[0]<2): #only one day
	    if (group==1):
	      #score_segments=(pre_obs.mean()*1)/refl.mean()
	      score_segments=(pre_obs.mean()*1)
	    else:
	      score_segments=(pre_obs.mean()*1)/pow(refl.mean(),power_refl)
	      
	  else:
	    
	    num_days=numpy.zeros((days.shape[0]-1), numpy.int32)
	    pre_obs_vector_mean=numpy.zeros((pre_obs.shape[0]-1), numpy.float)
	    refl_vector_mean=numpy.zeros((refl.shape[0]-1), numpy.float)	    
	    for kk in range(0,days.shape[0]-1):
	      num_days[kk]=days[kk+1]-days[kk]
	      pre_obs_values=pre_obs[kk:kk+1+1]
	      refl_values=refl[kk:kk+1+1]
	      
	      pre_obs_vector_mean[kk]=pre_obs_values.mean() #mean predicted - observed
	      refl_vector_mean[kk]=refl_values.mean() #mean reflectance
	
	    #area_total_segments=pre_obs_vector*num_days
	    if (group==1):
	      #score_segments=(pre_obs_vector_mean*num_days)/refl_vector_mean
	      score_segments=(pre_obs_vector_mean*num_days)
	    else:
	      score_segments=(pre_obs_vector_mean*num_days)/pow(refl_vector_mean,power_refl)	    
	    
	   
	    
	
	    #area total is the sum of the area of the segments
	  score_area[uu]=score_segments.sum()
	  
	if (group==1): #with groups ############################################################################3
	  #check if the area selected are elements of the group (through the dates using )
	  area_included=numpy.zeros(days_all.shape[0], numpy.int)
	  #to identify which area correspond to each value
	    
	  area_id_days=[]	
	  area_id_daily=[]
	  area_id_refl=[]
	  score_id_area=[]
	  area_id_pre_obsv=[]
	  area_id_predicted=[]
	  for tt in range (0,days_all.shape[0]):
	    #common=numpy.in1d(days_all[tt],days_group) #intersect between the two dates arrays (one area by one area)
	    common=numpy.in1d(days_all[tt],days_group_subgroup_candidates) #intersect between the two dates arrays (one area by one area)
	    common_2=numpy.in1d(days_group_subgroup_candidates,days_all[tt]) #to link each date to each area
	    common_index=numpy.where(common==True)
	    common2_index=numpy.where(common_2==True)
	    
	    if (common_index[0].shape[0]==0): #doesnt have any date group elements included
	      area_included[tt]=0
	      area_id_days.append(area_included[tt])
	      area_id_daily.append(area_included[tt])
	      area_id_refl.append(area_included[tt])
	      area_id_pre_obsv.append(area_included[tt])
	      area_id_predicted.append(area_included[tt])
	      score_id_area.append(area_included[tt])
	    else: #some date group elements included
	      area_included[tt]=1
	      area_id_days.append(days_group_subgroup_candidates[common2_index[0]])
	      
	      area_id_daily.append(daily_diff_group_sub_group_candidates[common2_index[0]])
	      area_id_refl.append(refl_values_group_subgroup_candidates[common2_index[0]])
	      area_id_pre_obsv.append(diff_pre_obs_group_sub_group_candidates[common2_index[0]])
	      area_id_predicted.append(predicted_values_group_subgroup_candidates[common2_index[0]])
	      
	      score_id_area.append(score_area[tt])
	  #score only with those areas with date values belonging to group
	  
	  score_area_group=score_area[area_included==1]
	  
	  #area_id_days=numpy.array(area_id_days)
	  #area_id_daily=numpy.array(area_id_daily)
	  #area_id_refl=numpy.array(area_id_refl)
	  #score_id_area=numpy.array(score_id_area)
	 
	  if (score_area_group.shape[0]==0): #empty- none of areas have candidates
	    print "none of areas have candidates"
	    sys.exit() #exit the script
	  else:
	  
	    refl_selected_group=refl_all[area_included==1]
	    days_selected_group=days_all[area_included==1]
	    predicted_selected_group=predicted_all[area_included==1]
	    pre_obs_selected_group=pre_obs_all[area_included==1]
	    daily_diff_selected_group=daily_diff_all[area_included==1]
	else:
	
	
	  #without groups
	  score_area_group=score_area
	  
	  refl_selected_group=refl_all
	  days_selected_group=days_all
	  predicted_selected_group=predicted_all
	  pre_obs_selected_group=pre_obs_all
	  daily_diff_selected_group=daily_diff_all	
	
	
	
	  ############################################################################################################### 
	
	score_area_index=numpy.where(area_included==1)
	
	
	
	
	score_final_all=[]
	max_score_index=numpy.zeros(len(score_id_area),numpy.float)
	max_score=numpy.zeros(len(score_id_area),numpy.float)
	for yy in range(0,len(score_id_area)):
	
	  score_final_all.append(pow(score_id_area[yy],1.0/2)-area_id_daily[yy]-area_id_refl[yy]) #score
	 
	
	#Where is the maximum score for date (maximum of a list)
	
	max_score=numpy.zeros(score_area_index[0].shape[0],numpy.float)
	for yy in range(0,score_area_index[0].shape[0]):
	  
	  #score_final_all.append(pow(score_area[score_area_index[0][yy]],1.0/2)-area_id_daily[score_area_index[0][yy]]-area_id_refl[score_area_index[0][yy]])
	  max_score[yy]=score_final_all[score_area_index[0][yy]].max()
      
      
      
	index_max_score_date_area=numpy.where(max_score==max_score.max()) 
	
	value_refl=area_id_refl[score_area_index[0][index_max_score_date_area[0]]]
	value_daily=area_id_daily[score_area_index[0][index_max_score_date_area[0]]]
	value_days=area_id_days[score_area_index[0][index_max_score_date_area[0]]]
	value_pre_obsv=area_id_pre_obsv[score_area_index[0][index_max_score_date_area[0]]]
	value_predicted=area_id_predicted[score_area_index[0][index_max_score_date_area[0]]]
	score=score_final_all[score_area_index[0][index_max_score_date_area[0]]]
	
	#area where is the maximum
	
	value_score_area=score_final_all[score_area_index[0][index_max_score_date_area[0]]] #values of scores of observations within the area
	
	index_maximum=numpy.where(max_score[index_max_score_date_area[0]]==value_score_area) #where is the maximum score ?
	score_date=score[index_maximum[0][0]]
	  
	'''
	#CHOOSE THE BEST SCORING (HIGHEST AREA SCORE WITHIN AREA GROUP)
		  
	
	selected_candidate=numpy.where(score_area_group==score_area_group.max())	
	'''
	
	score_area_selected=score_area_group[index_max_score_date_area[0]] #score of the area selected
	selected_candidate_area=index_max_score_date_area
	
	
	
	
	
	if(len(score_final_all)==1): #just one area choosen
	 
	  index_max_score_date_area=numpy.where(score_final_all==score_final_all.max())
	  scores_candidates=score_final_all[0][index_max_score_date_area[1][0]]
	  index_max_score_date=index_max_score_date_area
	  score_date=scores_candidates
	  
	  selected_candidate_area=index_max_score_date_area #area selected
	  
	  #VARIABLES OF THE AREA=REMAINS
	  refl_area_selected=numpy.array(refl_selected_group[selected_candidate_area[0][0]])
	  days_group_area_selected=numpy.array(days_selected_group[selected_candidate_area[0][0]])
	  diff_pre_obs_area_selected=numpy.array(pre_obs_selected_group[selected_candidate_area[0][0]])
	  predicted_area_selected=numpy.array(predicted_selected_group[selected_candidate_area[0][0]])
	  daily_diff_area_selected=numpy.array(daily_diff_selected_group[selected_candidate_area[0][0]])	  
	  
	  #VARIABLES FOR THE DATE
	  refl_selected=value_refl[index_maximum[1][0]]
	  day_selected=value_days[index_maximum[1][0]]
	  diff_pre_obs_selected=value_pre_obsv[index_maximum[1][0]]
	  daily_diff_selected=value_daily[index_maximum[1][0]]
	  predicted_selected=value_predicted[index_maximum[1][0]]	  
	  
	else:
	  
	
	  ##choose the highest
	  #index_max_score_date_area=numpy.where(max_score==max_score.max()) #which area has the date with the highest score
	  
	  ##in that area where is the maximum date?
	  
	  #scores_candidates=score_final_all[index_max_score_date_area[0][0]]
	  
	  #index_max_score_date=numpy.where(scores_candidates==scores_candidates.max())
	  #score_date=scores_candidates[index_max_score_date[0][0]]
	  
	  
	  #selected_candidate_area=index_max_score_date_area #area selected
	
	
	
	
	  #VARIABLES for the area selected############################################################################
	  
	  refl_area_selected=numpy.array(refl_selected_group[selected_candidate_area[0][0]])
	  days_group_area_selected=numpy.array(days_selected_group[selected_candidate_area[0][0]])
	  diff_pre_obs_area_selected=numpy.array(pre_obs_selected_group[selected_candidate_area[0][0]])
	  predicted_area_selected=numpy.array(predicted_selected_group[selected_candidate_area[0][0]])
	  daily_diff_area_selected=numpy.array(daily_diff_selected_group[selected_candidate_area[0][0]])	
	  #score_area_selected=score_area_group[selected_candidate[0][0]]
	
	
	
	  ##calculate the date inside the area choosen
	  
	  #score_final_all=diff_pre_obs_area_selected-daily_diff_area_selected-refl_area_selected
	  
		  
	  #selected_candidate=numpy.where(score_final_all==score_final_all.max())
	  
	  
	  #VARIABLES FOR THE DATE############################################################################
	  '''
	  refl_selected=refl_area_selected[selected_candidate[0][0]]
	  day_selected=days_group_area_selected[selected_candidate[0][0]]
	  diff_pre_obs_selected=diff_pre_obs_area_selected[selected_candidate[0][0]]
	  daily_diff_selected=daily_diff_area_selected[selected_candidate[0][0]]
	  predicted_selected=predicted_area_selected[selected_candidate[0][0]]
	  score_date=score_final_all[selected_candidate[0][0]]
	  '''
	  
	  refl_selected=value_refl[index_maximum[0][0]]
	  day_selected=value_days[index_maximum[0][0]]
	  diff_pre_obs_selected=value_pre_obsv[index_maximum[0][0]]
	  daily_diff_selected=value_daily[index_maximum[0][0]]
	  predicted_selected=value_predicted[index_maximum[0][0]]
	
	
	
	#score_date=score_final_all[selected_candidate[0][0]]	
	
	
	##b) Check where is the first negative----diff_pre_obs_group_negative[1][0]
	#diff_pre_obs_group_negative=numpy.where(diff_pre_obs_group<0)
	
	
	##same size of diff_pre_obs_group_positive
	#diff_pre_obs_group_sub_group=diff_pre_obs_group[0][diff_pre_obs_group_positive[1]]
	#daily_diff_group_sub_group=daily_diff_group[0][diff_pre_obs_group_positive[1]]
	#days_group_subgroup=days_group[0][diff_pre_obs_group_positive[1]]
	#refl_values_group_subgroup=refl_values_group[0][diff_pre_obs_group_positive[1]]
	#final_predicted_group_subgroup=final_predicted_group[0][diff_pre_obs_group_positive[1]]
	
	
	#diff_pre_obs_group_sub_group_candidates=diff_pre_obs_group_sub_group
	#daily_diff_group_sub_group_candidates=daily_diff_group_sub_group
	#days_group_subgroup_candidates=days_group_subgroup
	#refl_values_group_subgroup_candidates=refl_values_group_subgroup	
	#final_predicted_group_subgroup_candidates=final_predicted_group_subgroup
	
	#'''
	
	##choose the values in the previous group (a) where predicted-observed is positive (diff_pre_obs_group>0) and daily difference is negative (daily_diff_group<0)
	#daily_diff_sub_group_negative=numpy.where(daily_diff_group_sub_group<0) 
	
	##POTENTIAL CANDIDATES TO SCORING-same size of daily_diff_sub_group_negative
	#diff_pre_obs_group_sub_group_candidates=diff_pre_obs_group_sub_group[daily_diff_sub_group_negative]
	#daily_diff_group_sub_group_candidates=daily_diff_group_sub_group[daily_diff_sub_group_negative]
	#days_group_subgroup_candidates=days_group_subgroup[daily_diff_sub_group_negative]
	#refl_values_group_subgroup_candidates=refl_values_group_subgroup[daily_diff_sub_group_negative]
	#final_predicted_group_subgroup_candidates=final_predicted_group_subgroup[daily_diff_sub_group_negative]
	#'''
	
	##CALCULATE THE AREA UNTIL THE PREDICTED-OBSERVED IS NEGATIVE (where predicted-observed is positive)
	##area_individual=numpy.zeros((daily_diff_group_sub_group_candidates.shape[0]), numpy.float) #create a vector with the same size of the number of candidates
	
	#areas_pre_obs=[]
	#areas_refl=[] 
	#areas_days=[]	
	#areas_predicted=[]
	#score_area=numpy.zeros((diff_pre_obs_group_sub_group_candidates.shape[0]), numpy.float)
	
	##values_predicted_m_observed=[]
	
	
	########CONTEMPLAR CASOS ONDE NAO HA CANDIDATOS--SO PLOTAR###################
	
	#if (diff_pre_obs_group_sub_group_candidates.shape[0]==0):#no candidates
	  
	  #reflectance=[]
	  #percentage_minimum=[]
	  #date=[]
	  #pred_minus_observ=[]
	  #daily_diff_selected=[]  
	  #temporal_unc=[]
	  
	  #print "no candidates in this time series"
	#else:
	  
	
	  #for z in range(0,diff_pre_obs_group_sub_group_candidates.shape[0]):
	   
		
		
	    ##values_predicted_m_observed=diff_pre_obs_group[0][diff_pre_obs_group_positive[1][daily_diff_sub_group_negative[0][i]]]
	    ###index_new=daily_diff_sub_group_negative[0][z]
	    #index_new=diff_pre_obs_group_positive[1][z]
	    #date_ini=days_group[0][diff_pre_obs_group_positive[1][index_new]] #date of the first selected candidate
	    #predicted=[]
	    #pre_obs=[]
	    #refl=[]
	    #days=[]
	    #area_segment=[]
	    #score_individual_area=[]
	    ##for jj in range(diff_pre_obs_group_positive[1][index_new],diff_pre_obs_group.shape[1]): #all grou
	    #for jj in range(index_new,diff_pre_obs_group.shape[1]): #all group
	      
	      ##If the value of difference between predicted or observed is negative or the time differencee is bigger than 16	      
  
	      ##if (diff_pre_obs_group[0][jj]<0) or (days_group[0][jj]-days_group[0][diff_pre_obs_group_positive[1][index_new]]>16):
	      #if (diff_pre_obs_group[0][jj]<0) or (days_group[0][jj]-days_group[0][index_new]>window_size):
		#break
	      #else:
	      
		#diffpre_obs_orginal=diff_pre_obs_group[0][jj]
		#refl_original=refl_values_group[0][jj]
		#days_original=days_group[0][jj]
		#predicted_original=final_predicted_group[0][jj]
		
		#predicted.append(predicted_original)
		#pre_obs.append(diffpre_obs_orginal)
		#refl.append(refl_original)
		#days.append(days_original)
		
		
		
	    ##################################################################
	    #pre_obs=numpy.array(pre_obs) #convert list to array
	    #refl=numpy.array(refl)#convert list to array
	    #days=numpy.array(days) #convert list to array
	    #predicted=numpy.array(predicted)#convert list to array
	    
	    #if (days.shape[0]<2): #only one day
	      #score_individual_area=(pre_obs.mean()*1)/refl.mean()
	      ##appending all the values inside each area
	      #areas_pre_obs.append(pre_obs)
	      #areas_refl.append(refl) 
	      #areas_days.append(days)
	      #areas_predicted.append(predicted)
	    #else:
	      
	    
	      #date_end=days_group[0][jj-1] #date of the final curve
	      #num_days=int(date_end-date_ini)+1#number of days of the integral	    
	      
	      ##for each area
	      #score_individual_area=(pre_obs.mean()*num_days)/refl.mean() #score of each area
	      
	      ##appending all the values inside each area
	      #areas_pre_obs.append(pre_obs)
	      #areas_refl.append(refl) 
	      #areas_days.append(days)
	      #areas_predicted.append(predicted)
	    
	    ##array with all the scores for each individual area
	    #score_area[z]=score_individual_area 
	    #'''
	    #if (days.shape[0]<2): #only one day
	      ##area_total_segments=pre_obs.mean()*1
	      #score_segments=(pre_obs.mean()*num_days)/refl.mean() #just one value
	    #else:
	      
	      ##calculate individual areas/scores per segment
	    
	     
	      #num_days=numpy.zeros((days.shape[0]-1), numpy.int32)
	      #pre_obs_vector_mean=numpy.zeros((pre_obs.shape[0]-1), numpy.float)
	      #refl_vector_mean=numpy.zeros((refl.shape[0]-1), numpy.float)
	      #for kk in range(0,days.shape[0]-1):
		#num_days[kk]=days[kk+1]-days[kk]
		#pre_obs_values=pre_obs[kk:kk+1+1]
		#refl_values=refl[kk:kk+1+1]
		#pre_obs_vector_mean[kk]=pre_obs_values.mean() #mean predicted - observed
		#refl_vector_mean[kk]=refl_values.mean() #mean reflectance
	      
	      ##area_total_segments=pre_obs_vector*num_days
	      #score_segments=(pre_obs_vector_mean*num_days)/refl_vector_mean
	      #'''
	      
	      ##date_end=days_group[0][jj] #date of the final curve
	      ##num_days=int(date_end-date_ini)+1#number of days of the integral
	      
	      ##NOTE : se nao houver pre-obs negativos a seguir ao candidato--fica o ultimo valor registado!!!!#######################################
	      
	     
	    ##area_individual[z]=pre_obs.mean()*num_days
	    ###score_area[z]=score_segments.sum()
	    ##area_individual[z]=area_total_segments.sum()
	    ##score_individual[z]=area_individual[z]/refl.mean() #SAME SIZE OF daily_diff_group_sub_group_candidates
	  
	  ##CHOOSE THE BEST SCORING (HIGHEST AREA SCORE)
	  
	  
	  #selected_candidate=numpy.where(score_area==score_area.max())
	  
	  
	  ##VARIABLES############################################################################
	  
	  #refl_selected=refl_values_group_subgroup_candidates[selected_candidate[0][0]]
	  #days_group_selected=days_group_subgroup_candidates[selected_candidate[0][0]]
	  #diff_pre_obs_selected=diff_pre_obs_group_sub_group_candidates[selected_candidate[0][0]]
	  #daily_diff_selected=daily_diff_group_sub_group_candidates[selected_candidate[0][0]]
	   
	   
	   
	  #'''
	  ########################CALCULATE VARIABLES##########################################################
	  ##based on the group of data and the diff_pre_obs_group_max_index or daily_diff_group_max_index
	
	  ##value of the reflectance selected-REFLECTANCE
	  ##index_true=numpy.where(diff_pre_obs_group_max_index==True)
	  #index_true=numpy.where(daily_diff_group_max_index==True)
	  
	  ##check if the index has more than one value (in case when we will have more than one maximum)-choose the observation with the earliest date
  
		  
	  #refl_selected=refl_values_group[0][index_true[1][0]] #even if its more than one element will choose always the the observations with the earliest date
	  #'''
	  ##percentage of all values in the core year (common_days_novoy_year_core) below refl_selected - PERCENTAGE_MINIMUM
	notnan_index=numpy.where(numpy.isfinite(common_days_novoy_year_core)==True)
	index_below_refl_selected=numpy.where(common_days_novoy_year_core[notnan_index]<refl_selected)
	per_below_refl_selected  =(index_below_refl_selected[0].shape[0]/float(common_days_novoy_year_core[notnan_index].shape[0]))*100.00
	  #'''
	  ##date selected - DATE
	  #days_group_selected=days_group[0][index_true[1][0]]
	  #days_group_selected_jul=int(DT.datetime.strftime( num2date(days_group_selected), "%j" ))
	  
	  ##difference between predicted and observed on the selected date - PRED_MINUS_OBSERV-based on the biggest negative daily difference index
	  #diff_pre_obs_selected=diff_pre_obs_group[0][index_true[1][0]]
  
  
	  ##daily difference value
	  
	  #daily_diff_selected=daily_diff_group[0][index_true[1][0]]
  
	  #'''
	  
	  
	#name of the variables
	reflectance=refl_selected
	percentage_minimum=per_below_refl_selected
	date=day_selected
	pred_minus_observ=diff_pre_obs_selected
	daily_diff_selected=daily_diff_selected
	predicted_selected=predicted_selected
	##################################################################################################
  
	########################################Temporal uncertainty#####################################
	'''
	#Identify where is the max value in the original daily differences
	daily_diff_original_max_index=daily_diff_group_max==diff_reflect2
	index_true_daily_diff=numpy.where(daily_diff_original_max_index==True)
	temporal_unc=val_dias[index_true_daily_diff[0]]-val_dias[index_true_daily_diff[0]-1]-1 #difference between days minus one
	temporal_unc=temporal_unc[0]
	
	'''
	
	date_original_index=date==val_dias
	index_true_daily_diff=numpy.where(date_original_index==True)
	temporal_unc=val_dias[index_true_daily_diff[0]]-val_dias[index_true_daily_diff[0]-1]-1 #difference between days minus one
	temporal_unc=temporal_unc[0]	  


	##############################write the file with the novoy4######################################
	WorkDir = '/media/HardDisk/Datasets/PIXELS_PROFILES/PROFILES_ALGORITHM_VERSION_4/'
	
	#files
	FirefileName= os.path.basename(File1)[0:-4] + '_' + str(i) + '_' + str(j) + '_profile_vs4.txt'
	FirefileName_mat = os.path.basename(File1)[0:-4] + '_' + str(i) + '_' + str(j) + '_profile_vs4.mat'
	FirefileName_png = os.path.basename(File1)[0:-4] + '_' + str(i) + '_' + str(j) + '_profile_vs4.png'
		  
	      
	
	profile_vector = numpy.vstack((GregDay_all,GregDay_all_jul, val_nir_nan,final_predicted,diff_pre_obs,numpy.tile(i,GregDay_all_jul.shape[0]),numpy.tile(j,GregDay_all_jul.shape[0]))).T
	
	#save
	numpy.savetxt(WorkDir + '/' + FirefileName,profile_vector,fmt='%d \t%d \t%5.3f \t%5.3f \t%5.3f \t%d \t%d')
		
	
	#save it to matfile
	#scipy.io.savemat(WorkDir+'/'+ FirefileName_mat, {'Day_ori':val_dias, 'Day_jul':val_dias_jul, 'Min':novoy4, 'Final_predicted':final_predicted, 'Diff_pre_obs':diff_pre_obs,'row':i, 'col':j})
	
	###############VARIABLES#############################################
	reflectance=reflectance
	percentage_minimum=percentage_minimum
	date= date
	pred_minus_observ=pred_minus_observ


	#####################################################################

	
	  		

	
	
	
	#return  days_group_subgroup_candidates,refl_values_group_subgroup_candidates,diff_pre_obs_group_sub_group_candidates,daily_diff_group_sub_group_candidates,daily_diff_selected,temporal_unc,num_valids_inside_window,percent_parameter,daily_diff_group,val_nir,density_alldata,window,days_group,diff_pre_obs_group,refl_values_group,year,numlinhas,numcolunas,GregDay_all,GregDay_all_jul,val_nir_nan,val_dias,novoy3,novoy4,novoy_median,diff_reflect2,diff_reflect2_nan,diff_reflect2_trimmed,final_predicted,diff_pre_obs,common_days_novox_year_core,common_days_novoy_year_core,diff_reflected2_core_year,final_predicted_core_year,diff_pre_obs_core_year,common_days_novox_year_core_index,reflectance,percentage_minimum,date,pred_minus_observ,value_ref,index_not_nan_diff_pre_obs_nan,index_not_nan_val_nir_nan,index_not_nan_final_predicted_nan
	return  novoy,novoy1,novoy2,novoy3_1,novoy3_2,novoy3_3,days_group_subgroup_candidates,refl_values_group_subgroup_candidates,diff_pre_obs_group_sub_group_candidates,daily_diff_group_sub_group_candidates,group,score_date,score_area_selected,daily_diff_area_selected,predicted_area_selected,diff_pre_obs_area_selected,days_group_area_selected,refl_area_selected,predicted_all,pre_obs_all,refl_all,days_all,daily_diff_all,daily_diff_selected,predicted_selected,temporal_unc,num_valids_inside_window,percent_parameter,daily_diff_group,val_nir,density_alldata,window,days_group,diff_pre_obs_group,refl_values_group,year,numlinhas,numcolunas,GregDay_all,GregDay_all_jul,val_nir_nan,val_dias,novoy3,novoy4,novoy_median,diff_reflect2,diff_reflect2_nan,diff_reflect2_trimmed,final_predicted,diff_pre_obs,common_days_novox_year_core,common_days_novoy_year_core,diff_reflected2_core_year,final_predicted_core_year,diff_pre_obs_core_year,common_days_novox_year_core_index,reflectance,percentage_minimum,date,pred_minus_observ,value_ref,index_not_nan_diff_pre_obs_nan,index_not_nan_val_nir_nan,index_not_nan_final_predicted_nan
	#return  novoy,novoy1,novoy2,days_group_subgroup_candidates,refl_values_group_subgroup_candidates,diff_pre_obs_group_sub_group_candidates,daily_diff_group_sub_group_candidates,group,score_date,score_area_selected,daily_diff_area_selected,predicted_area_selected,diff_pre_obs_area_selected,days_group_area_selected,refl_area_selected,predicted_all,pre_obs_all,refl_all,days_all,daily_diff_all,daily_diff_selected,predicted_selected,temporal_unc,num_valids_inside_window,percent_parameter,daily_diff_group,val_nir,density_alldata,window,days_group,diff_pre_obs_group,refl_values_group,year,numlinhas,numcolunas,GregDay_all,GregDay_all_jul,val_nir_nan,val_dias,novoy3,novoy4,novoy_median,diff_reflect2,diff_reflect2_nan,diff_reflect2_trimmed,final_predicted,diff_pre_obs,common_days_novox_year_core,common_days_novoy_year_core,diff_reflected2_core_year,final_predicted_core_year,diff_pre_obs_core_year,common_days_novox_year_core_index,reflectance,percentage_minimum,date,pred_minus_observ,value_ref,index_not_nan_diff_pre_obs_nan,index_not_nan_val_nir_nan,index_not_nan_final_predicted_nan
      
      
def verCPD_GT(data_nir,numleituras, diajul, cellref,numlinhas,numcolunas,dist,GregDay, YearJulianDay,ano,year,indices_core_year,Greg_day_core_year,plot_nb,i,j):
	
	#days_group_subgroup_candidates,refl_values_group_subgroup_candidates,diff_pre_obs_group_sub_group_candidates,daily_diff_group_sub_group_candidates,daily_diff_selected,temporal_unc,num_valids_inside_window,percent_parameter,daily_diff_group,val_nir,density_alldata,window,days_group,diff_pre_obs_group,refl_values_group,year,numlinhas,numcolunas,GregDay_all,GregDay_all_jul,val_nir_nan,val_dias,novoy3,novoy4,novoy_median,diff_reflect2,diff_reflect2_nan,diff_reflect2_trimmed,final_predicted,diff_pre_obs,common_days_novox_year_core,common_days_novoy_year_core,diff_reflected2_core_year,final_predicted_core_year,diff_pre_obs_core_year,common_days_novox_year_core_index,reflectance,percentage_minimum,date,pred_minus_observ,value_ref,index_not_nan_diff_pre_obs_nan,index_not_nan_val_nir_nan,index_not_nan_final_predicted_nan=calCPD_GT(data_nir,numleituras, diajul, cellref, dist,GregDay, YearJulianDay,ano,year,indices_core_year,Greg_day_core_year,i,j)
	novoy,novoy1,novoy2,novoy3_1,novoy3_2,novoy3_3,days_group_subgroup_candidates,refl_values_group_subgroup_candidates,diff_pre_obs_group_sub_group_candidates,daily_diff_group_sub_group_candidates,group,score_date,score_area_selected,daily_diff_area_selected,predicted_area_selected,diff_pre_obs_area_selected,days_group_area_selected,refl_area_selected,predicted_all,pre_obs_all,refl_all,days_all,daily_diff_all,daily_diff_selected,predicted_selected,temporal_unc,num_valids_inside_window,percent_parameter,daily_diff_group,val_nir,density_alldata,window,days_group,diff_pre_obs_group,refl_values_group,year,numlinhas,numcolunas,GregDay_all,GregDay_all_jul,val_nir_nan,val_dias,novoy3,novoy4,novoy_median,diff_reflect2,diff_reflect2_nan,diff_reflect2_trimmed,final_predicted,diff_pre_obs,common_days_novox_year_core,common_days_novoy_year_core,diff_reflected2_core_year,final_predicted_core_year,diff_pre_obs_core_year,common_days_novox_year_core_index,reflectance,percentage_minimum,date,pred_minus_observ,value_ref,index_not_nan_diff_pre_obs_nan,index_not_nan_val_nir_nan,index_not_nan_final_predicted_nan=calCPD_GT(plot_nb,data_nir,numleituras, diajul, cellref, dist,GregDay, YearJulianDay,ano,year,indices_core_year,Greg_day_core_year,i,j)
	#novoy,novoy1,novoy2,days_group_subgroup_candidates,refl_values_group_subgroup_candidates,diff_pre_obs_group_sub_group_candidates,daily_diff_group_sub_group_candidates,group,score_date,score_area_selected,daily_diff_area_selected,predicted_area_selected,diff_pre_obs_area_selected,days_group_area_selected,refl_area_selected,predicted_all,pre_obs_all,refl_all,days_all,daily_diff_all,daily_diff_selected,predicted_selected,temporal_unc,num_valids_inside_window,percent_parameter,daily_diff_group,val_nir,density_alldata,window,days_group,diff_pre_obs_group,refl_values_group,year,numlinhas,numcolunas,GregDay_all,GregDay_all_jul,val_nir_nan,val_dias,novoy3,novoy4,novoy_median,diff_reflect2,diff_reflect2_nan,diff_reflect2_trimmed,final_predicted,diff_pre_obs,common_days_novox_year_core,common_days_novoy_year_core,diff_reflected2_core_year,final_predicted_core_year,diff_pre_obs_core_year,common_days_novox_year_core_index,reflectance,percentage_minimum,date,pred_minus_observ,value_ref,index_not_nan_diff_pre_obs_nan,index_not_nan_val_nir_nan,index_not_nan_final_predicted_nan=calCPD_GT(plot_nb,data_nir,numleituras, diajul, cellref, dist,GregDay, YearJulianDay,ano,year,indices_core_year,Greg_day_core_year,i,j)
	
	#verSerie_GT_all_profile(days_group_subgroup_candidates,refl_values_group_subgroup_candidates,diff_pre_obs_group_sub_group_candidates,daily_diff_group_sub_group_candidates,daily_diff_selected,temporal_unc,num_valids_inside_window,percent_parameter,daily_diff_group,val_nir,density_alldata,window,days_group,diff_pre_obs_group,refl_values_group,year,numlinhas,numcolunas,GregDay_all,GregDay_all_jul,val_nir_nan,val_dias,novoy3,novoy4,novoy_median,diff_reflect2,diff_reflect2_nan,diff_reflect2_trimmed,final_predicted,diff_pre_obs,common_days_novox_year_core,common_days_novoy_year_core,diff_reflected2_core_year,final_predicted_core_year,diff_pre_obs_core_year,common_days_novox_year_core_index,reflectance,percentage_minimum,date,pred_minus_observ,value_ref,index_not_nan_diff_pre_obs_nan,index_not_nan_val_nir_nan,index_not_nan_final_predicted_nan,i,j)

	verSerie_GT_all_profile(novoy,novoy1,novoy2,novoy3_1,novoy3_2,novoy3_3,days_group_subgroup_candidates,refl_values_group_subgroup_candidates,diff_pre_obs_group_sub_group_candidates,daily_diff_group_sub_group_candidates,group,score_date,score_area_selected,daily_diff_area_selected,predicted_area_selected,diff_pre_obs_area_selected,days_group_area_selected,refl_area_selected,predicted_all,pre_obs_all,refl_all,days_all,daily_diff_all,daily_diff_selected,predicted_selected,temporal_unc,num_valids_inside_window,percent_parameter,daily_diff_group,val_nir,density_alldata,window,days_group,diff_pre_obs_group,refl_values_group,year,numlinhas,numcolunas,GregDay_all,GregDay_all_jul,val_nir_nan,val_dias,novoy3,novoy4,novoy_median,diff_reflect2,diff_reflect2_nan,diff_reflect2_trimmed,final_predicted,diff_pre_obs,common_days_novox_year_core,common_days_novoy_year_core,diff_reflected2_core_year,final_predicted_core_year,diff_pre_obs_core_year,common_days_novox_year_core_index,reflectance,percentage_minimum,date,pred_minus_observ,value_ref,index_not_nan_diff_pre_obs_nan,index_not_nan_val_nir_nan,index_not_nan_final_predicted_nan,i,j)
	#verSerie_GT_all_profile(novoy,novoy1,novoy2,days_group_subgroup_candidates,refl_values_group_subgroup_candidates,diff_pre_obs_group_sub_group_candidates,daily_diff_group_sub_group_candidates,group,score_date,score_area_selected,daily_diff_area_selected,predicted_area_selected,diff_pre_obs_area_selected,days_group_area_selected,refl_area_selected,predicted_all,pre_obs_all,refl_all,days_all,daily_diff_all,daily_diff_selected,predicted_selected,temporal_unc,num_valids_inside_window,percent_parameter,daily_diff_group,val_nir,density_alldata,window,days_group,diff_pre_obs_group,refl_values_group,year,numlinhas,numcolunas,GregDay_all,GregDay_all_jul,val_nir_nan,val_dias,novoy3,novoy4,novoy_median,diff_reflect2,diff_reflect2_nan,diff_reflect2_trimmed,final_predicted,diff_pre_obs,common_days_novox_year_core,common_days_novoy_year_core,diff_reflected2_core_year,final_predicted_core_year,diff_pre_obs_core_year,common_days_novox_year_core_index,reflectance,percentage_minimum,date,pred_minus_observ,value_ref,index_not_nan_diff_pre_obs_nan,index_not_nan_val_nir_nan,index_not_nan_final_predicted_nan,i,j)
	  
	  
	  
#def verSerie_GT_all_profile(days_group_subgroup_candidates,refl_values_group_subgroup_candidates,diff_pre_obs_group_sub_group_candidates,daily_diff_group_sub_group_candidates,daily_diff_selected,temporal_unc,num_valids_inside_window,percent_parameter,daily_diff_group,val_nir,density_alldata,window,days_group,diff_pre_obs_group,refl_values_group,year,numlinhas,numcolunas,GregDay_all,GregDay_all_jul,val_nir_nan,val_dias,novoy3,novoy4,novoy_median,diff_reflect2,diff_reflect2_nan,diff_reflect2_trimmed,final_predicted,diff_pre_obs,common_days_novox_year_core,common_days_novoy_year_core,diff_reflected2_core_year,final_predicted_core_year,diff_pre_obs_core_year,common_days_novox_year_core_index,reflectance,percentage_minimum,date,pred_minus_observ,value_ref,index_not_nan_diff_pre_obs_nan,index_not_nan_val_nir_nan,index_not_nan_final_predicted_nan,i,j):
def verSerie_GT_all_profile(novoy,novoy1,novoy2,novoy3_1,novoy3_2,novoy3_3,days_group_subgroup_candidates,refl_values_group_subgroup_candidates,diff_pre_obs_group_sub_group_candidates,daily_diff_group_sub_group_candidates,group,score_date,score_area_selected,daily_diff_area_selected,predicted_area_selected,diff_pre_obs_area_selected,days_group_area_selected,refl_area_selected,predicted_all,pre_obs_all,refl_all,days_all,daily_diff_all,daily_diff_selected,predicted_selected,temporal_unc,num_valids_inside_window,percent_parameter,daily_diff_group,val_nir,density_alldata,window,days_group,diff_pre_obs_group,refl_values_group,year,numlinhas,numcolunas,GregDay_all,GregDay_all_jul,val_nir_nan,val_dias,novoy3,novoy4,novoy_median,diff_reflect2,diff_reflect2_nan,diff_reflect2_trimmed,final_predicted,diff_pre_obs,common_days_novox_year_core,common_days_novoy_year_core,diff_reflected2_core_year,final_predicted_core_year,diff_pre_obs_core_year,common_days_novox_year_core_index,reflectance,percentage_minimum,date,pred_minus_observ,value_ref,index_not_nan_diff_pre_obs_nan,index_not_nan_val_nir_nan,index_not_nan_final_predicted_nan,i,j):
#def verSerie_GT_all_profile(novoy,novoy1,novoy2,days_group_subgroup_candidates,refl_values_group_subgroup_candidates,diff_pre_obs_group_sub_group_candidates,daily_diff_group_sub_group_candidates,group,score_date,score_area_selected,daily_diff_area_selected,predicted_area_selected,diff_pre_obs_area_selected,days_group_area_selected,refl_area_selected,predicted_all,pre_obs_all,refl_all,days_all,daily_diff_all,daily_diff_selected,predicted_selected,temporal_unc,num_valids_inside_window,percent_parameter,daily_diff_group,val_nir,density_alldata,window,days_group,diff_pre_obs_group,refl_values_group,year,numlinhas,numcolunas,GregDay_all,GregDay_all_jul,val_nir_nan,val_dias,novoy3,novoy4,novoy_median,diff_reflect2,diff_reflect2_nan,diff_reflect2_trimmed,final_predicted,diff_pre_obs,common_days_novox_year_core,common_days_novoy_year_core,diff_reflected2_core_year,final_predicted_core_year,diff_pre_obs_core_year,common_days_novox_year_core_index,reflectance,percentage_minimum,date,pred_minus_observ,value_ref,index_not_nan_diff_pre_obs_nan,index_not_nan_val_nir_nan,index_not_nan_final_predicted_nan,i,j):

	  
	  
	  import matplotlib.pyplot as plt
	  import GTr1AlgFunc_vs_3_3_wings_vs2
	  import numpy
	  import pylab as p
	  from matplotlib.dates import date2num
	  from matplotlib.dates import num2date
	  from scipy.stats import gaussian_kde
	  from matplotlib.patches import Rectangle
	  import matplotlib.collections as collections
	  import matplotlib.mlab as mlab
	  import math
	  from scipy import stats
	  #VARIABLES
	  previous_year  = int(year) - 1
	  following_year = int(year) + 1
	  previous_year = str (previous_year)
	  following_year = str (following_year)		  
	 	  
	  
	  ########PLOT FIGURE WITH THE STEPS OF C-O OPERATIONS########################
	  
	  
	  
	  
		
	   #find the index day of the selected point 
	   
	   
	#preparing the ticks
	
	  x_ticks = [335 , 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 1]
	  x_ticks_2 = []
	  x_ticks_2.append (previous_year + GetStrDoY(x_ticks[0]))
	   
	  for l in range(1,len(x_ticks)-1):
		       
	    x_ticks_2.append (year + GetStrDoY(x_ticks[l]))
	     
	  x_ticks_2.append (following_year + GetStrDoY(x_ticks[-1])) 
		       
	  x_ticks_2 = numpy.array(x_ticks_2)    
		   
	  
	  #convert to greg days
	  x_ticks_greg = numpy.zeros(x_ticks_2.shape,numpy.int32)
	  for m in range(0, x_ticks_2.shape[0]):
		       
	    x_ticks_greg[m] = date2num(DT.datetime.strptime(x_ticks_2[m], "%Y%j"))	    
	   
		   
	  Months = ['dec' + previous_year[2:4],'jan' + year[2:4],'feb','mar','apr','mai','jun' + year[2:4],'jul','ago','sep','oct','nov','dec' + year[2:4],'jan' + following_year[2:4]]#14 months
	  Months = numpy.array(Months)	  
	  
	  
	  ########PLOT FIGURE WITH THE STEPS OF C-O OPERATIONS########################
	  plt.figure()
		    
	  fig=plt.gcf()#DO E RENATA
	  fig.set_size_inches(15.5,9.5)	
		   
	  if (daily_diff_group_sub_group_candidates.shape[0]==0): #NO CANDIDATES
	    ax1=plt.subplot(2,1,1)###################################################
	   
	    
	    ax1.plot(val_dias, val_nir, color=[0.72, 0.72, 0.72],linewidth=0.5,marker='.')  #original data
	    
	    
	    ax1.plot(val_dias, novoy3, color=[0.72, 0.72, 0.72],linewidth=2,marker='.',markersize=8)  #just to have all the information behind
	    
	    
	    
	    
	    #ax1.plot(GregDay_all[index_not_nan_val_nir_nan], val_nir_nan[index_not_nan_val_nir_nan], color='k',marker='.',markersize=4,label='original')
	    #ax1.plot(GregDay_all[index_not_nan_final_predicted_nan], final_predicted[index_not_nan_final_predicted_nan], color='r',marker='.',markersize=4,label='predicted')
	    #ax1.plot(GregDay_all,novoy_median, color='m',marker='.',markersize=4,label='median')
	     
  
  
	    ax1.plot(GregDay_all,val_nir_nan, color='k',marker='.',markersize=5,label='observed')
	    
	    ax1.plot(GregDay_all,final_predicted, color='r',marker='.',markersize=4,label='predicted')
		     
	    ax1.plot(days_group,refl_values_group, color='b',marker='o',markersize=4)
	    
	   
	    
	    
	    
	  
	    #set limits
	    plt.xlim([x_ticks_greg[0], x_ticks_greg[-1]+30])
	    minobserv=val_nir_nan[index_not_nan_val_nir_nan].min()
	    #ax1.set_ylim([minobserv-11, val_nir_nan[index_not_nan_val_nir_nan].max()])
	    
	     
		    
	    
	    #ticks
	    plt.xticks(x_ticks_greg,Months)
	
	
	    #put lines limiting the year
	    ax1.plot((x_ticks_greg[1],x_ticks_greg[1]),(minobserv-11,val_nir_nan[index_not_nan_val_nir_nan].max()) , color='b',label='limit year')#january year
	    
	    
	    ax1.plot((x_ticks_greg[-1],x_ticks_greg[-1]),(minobserv-11,val_nir_nan[index_not_nan_val_nir_nan].max()) , color='b') #january following year
	  
	    #lines
	    ax1.plot((GregDay_all[0], GregDay_all[-1]), (30, 30), color='r', linestyle='dashdot',label='Limit 30')
	    ax1.plot((GregDay_all[0], GregDay_all[-1]), (value_ref, value_ref), color='g', linestyle='dashdot',label='Limit' + str(percent_parameter) + '%')
	    
	    
	    
	    
  
	    ax1.grid(True)
	    
	    ax1.set_xlabel('Julian dates')
	    ax1.set_ylabel('NIR reflectance')	  
	    ax1.legend(loc='upper right',prop={'size':8})
  
	    
	    plt.title('Vs4_Four_variables -: NO CANDIDATES Tile: LL' + str(numlinhas) + '_' +'CC' + str(numcolunas) + '_row: ' +str(i) + ', col:' + str(j)+ ' ' + str(window) + 'days window'+ '_' + str(num_valids_inside_window) + 'val_window')    
	    
	    ax2=plt.subplot(2,1,2,sharex=ax1)    #link axes
	    ax2.plot(GregDay_all, diff_pre_obs,marker='.',markersize=4,color='k')  
	    ax2.plot(GregDay_all[index_not_nan_diff_pre_obs_nan], diff_pre_obs[index_not_nan_diff_pre_obs_nan], color=[0.72, 0.72, 0.72],marker='.',markersize=7,label='Pred-Obser')
	    ax2.plot(GregDay_all, diff_pre_obs,marker='.',markersize=4,color='k')
	    ax2.plot(GregDay_all, diff_reflect2_trimmed,marker='.',markersize=4,color='m')#trimmed median
	    
	    ax2.plot(GregDay_all, diff_reflect2_nan,marker='.',markersize=4,color='c',label='Daily diff.') #daily differences
	    
  
	    ax2.plot(days_group,diff_pre_obs_group, color='b',marker='o',markersize=4)
	    ax2.plot(days_group,daily_diff_group, color='r',marker='o',markersize=4)
	    
  
	    minobserv_2=numpy.nanmin(daily_diff_group)
	    
	    
	      
	    #ticks
	    plt.xticks(x_ticks_greg,Months)	  
  
	    ax2.grid(True)
	    ax2.set_xlabel('Julian dates')
	    ax2.set_ylabel('Daily  differences')	  
	    ax2.legend(loc='upper right',prop={'size':8})	     
		
	       
	    #plt.show()
	    plt.savefig('New_Plot_tiles_wings_vs4/VS4_1/Vs4_Four_variables_LL' + str(numlinhas) + '_'+ 'CC' + str(numcolunas) + '_row_' + str(i) + '_col_' + str(j) + '_' + str(year) + '_' + str(window)+ '_days' + '_' + str(num_valids_inside_window)+ '_val' + '.pdf',dpi=500)
	      
	      
	    
	  else:	  
	    ax1=plt.subplot(4,1,1)###################################################
		       
			
	    ax1.plot(val_dias, val_nir, color=[0.72, 0.72, 0.72],linewidth=0.5,marker='.')  #original data
	    ax1.plot(val_dias, novoy3, color='r',linewidth=1,marker='.',markersize=5,label='pass_2')  #just to have all the information behind
	    plt.title('Vs4_Four_variables_PAR -: Tile: LL' + str(numlinhas) + '_' +'CC' + str(numcolunas) + '_row: ' +str(i) + ', col:' + str(j)+ ' ' + str(window) + 'days window'+ '_' + str(num_valids_inside_window) + 'val_window' + ' ' + 'Minimum_values=' + str(group))    
	    
	    ax2=plt.subplot(4,1,2,sharex=ax1)
	    ax2.plot(val_dias, val_nir, color=[0.72, 0.72, 0.72],linewidth=0.5,marker='.')  #original data
	    ax2.plot(val_dias, novoy3_1, color='b',linewidth=1,marker='.',markersize=5,label='pass_3')  #just to have all the information behind
	    
	    ax3=plt.subplot(4,1,3,sharex=ax1)
	    ax3.plot(val_dias, val_nir, color=[0.72, 0.72, 0.72],linewidth=0.5,marker='.')  #original data
	    ax3.plot(val_dias, novoy3_2, color='g',linewidth=1,marker='.',markersize=5,label='pass_4')  #just to have all the information behind
	    
	    ax4=plt.subplot(4,1,4,sharex=ax1)
	    ax4.plot(val_dias, val_nir, color=[0.72, 0.72, 0.72],linewidth=0.5,marker='.')  #original data
	    ax4.plot(val_dias, novoy3_3, color='y',linewidth=2,marker='.',markersize=8,label='pass_5')  #just to have all the information behind	
	
	
	    #set limits
	    plt.xlim([x_ticks_greg[0], x_ticks_greg[-1]+30])
	    minobserv=val_nir_nan[index_not_nan_val_nir_nan].min()
	    #ax1.set_ylim([minobserv-11, val_nir_nan[index_not_nan_val_nir_nan].max()])
	    
	    #ticks
	    plt.xticks(x_ticks_greg,Months)
	
	
	    #put lines limiting the year
	    ax4.plot((x_ticks_greg[1],x_ticks_greg[1]),(minobserv-11,val_nir_nan[index_not_nan_val_nir_nan].max()) , color='b',label='limit year')#january year
	    
	    
	    ax4.plot((x_ticks_greg[-1],x_ticks_greg[-1]),(minobserv-11,val_nir_nan[index_not_nan_val_nir_nan].max()) , color='b') #january following year
	  
	    #lines
	    ax4.plot((GregDay_all[0], GregDay_all[-1]), (30, 30), color='r', linestyle='dashdot',label='Limit 30')
	    ax4.plot((GregDay_all[0], GregDay_all[-1]), (value_ref, value_ref), color='g', linestyle='dashdot',label='Limit' + str(percent_parameter) + '%')
	    plt.axvline(date,color='r',label='Selected ')
	    
	    
	    #candidates
	    for z in range(0,days_group_subgroup_candidates.shape[0]):
	      plt.axvline(days_group_subgroup_candidates[z],color='b',linestyle='dashdot')
	   
	    #days_group_subgroup_candidates,
	   # refl_values_group_subgroup_candidates,
	    #diff_pre_obs_group_sub_group_candidates,
	    #refl_values_group_subgroup_candidates,
	    
	    ax1.grid(True)
	    ax2.grid(True)
	    ax3.grid(True)
	    ax4.grid(True)
	    
	    ax4.set_xlabel('Julian dates')
	    
	    ax1.set_ylabel('NIR reflectance')	
	    ax2.set_ylabel('NIR reflectance')
	    ax3.set_ylabel('NIR reflectance')
	    ax4.set_ylabel('NIR reflectance')
	    
	    ax1.legend(loc='upper right',prop={'size':8})	    
	    ax2.legend(loc='upper right',prop={'size':8})	
	    ax3.legend(loc='upper right',prop={'size':8})	
	    ax4.legend(loc='upper right',prop={'size':8})	
		
		   
	    plt.savefig('New_Plot_tiles_wings_vs4/VS4_1/Vs4_3_par_step_Four_variables_LL' + str(numlinhas) + '_'+ 'CC' + str(numcolunas) + '_row_' + str(i) + '_col_' + str(j) + '_' + str(year) + '_' + str(window)+ '_days' + '_' + str(num_valids_inside_window)+ '_val_' +'min_' + str(group) + '.pdf',dpi=500)
	   
		   
	##################################################################################################################################################	   
	  plt.figure()
	
	  fig=plt.gcf()#DO E RENATA
	  fig.set_size_inches(15.5,9.5)		  
	  
	  
	  
	  
	  
	  #SEPARATE PIXELS WITH CANDIDATES AND PIXELS WITHOUT CANDIDATES
	  
	  if (daily_diff_group_sub_group_candidates.shape[0]==0): #NO CANDIDATES
	    ax1=plt.subplot(2,1,1)###################################################
	   
	    
	    ax1.plot(val_dias, val_nir, color=[0.72, 0.72, 0.72],linewidth=0.5,marker='.')  #original data
	    
	    
	    ax1.plot(val_dias, novoy3, color=[0.72, 0.72, 0.72],linewidth=2,marker='.',markersize=8)  #just to have all the information behind
	    
	    
	    
	    
	    #ax1.plot(GregDay_all[index_not_nan_val_nir_nan], val_nir_nan[index_not_nan_val_nir_nan], color='k',marker='.',markersize=4,label='original')
	    #ax1.plot(GregDay_all[index_not_nan_final_predicted_nan], final_predicted[index_not_nan_final_predicted_nan], color='r',marker='.',markersize=4,label='predicted')
	    #ax1.plot(GregDay_all,novoy_median, color='m',marker='.',markersize=4,label='median')
	     
  
  
	    ax1.plot(GregDay_all,val_nir_nan, color='k',marker='.',markersize=5,label='observed')
	    
	    ax1.plot(GregDay_all,final_predicted, color='r',marker='.',markersize=4,label='predicted')
		     
	    ax1.plot(days_group,refl_values_group, color='b',marker='o',markersize=4)
	    
	   
	    
	    
	    
	  
	    #set limits
	    plt.xlim([x_ticks_greg[0], x_ticks_greg[-1]+30])
	    minobserv=val_nir_nan[index_not_nan_val_nir_nan].min()
	    #ax1.set_ylim([minobserv-11, val_nir_nan[index_not_nan_val_nir_nan].max()])
	    
	     
		    
	    
	    #ticks
	    plt.xticks(x_ticks_greg,Months)
	
	
	    #put lines limiting the year
	    ax1.plot((x_ticks_greg[1],x_ticks_greg[1]),(minobserv-11,val_nir_nan[index_not_nan_val_nir_nan].max()) , color='b',label='limit year')#january year
	    
	    
	    ax1.plot((x_ticks_greg[-1],x_ticks_greg[-1]),(minobserv-11,val_nir_nan[index_not_nan_val_nir_nan].max()) , color='b') #january following year
	  
	    #lines
	    ax1.plot((GregDay_all[0], GregDay_all[-1]), (30, 30), color='r', linestyle='dashdot',label='Limit 30')
	    ax1.plot((GregDay_all[0], GregDay_all[-1]), (value_ref, value_ref), color='g', linestyle='dashdot',label='Limit' + str(percent_parameter) + '%')
	    
	    
	    
	    
  
	    ax1.grid(True)
	    
	    ax1.set_xlabel('Julian dates')
	    ax1.set_ylabel('NIR reflectance')	  
	    ax1.legend(loc='upper right',prop={'size':8})
  
	    
	    plt.title('Vs4_Four_variables -: NO CANDIDATES Tile: LL' + str(numlinhas) + '_' +'CC' + str(numcolunas) + '_row: ' +str(i) + ', col:' + str(j)+ ' ' + str(window) + 'days window'+ '_' + str(num_valids_inside_window) + 'val_window')    
	    
	    ax2=plt.subplot(2,1,2,sharex=ax1)    #link axes
	    ax2.plot(GregDay_all, diff_pre_obs,marker='.',markersize=4,color='k')  
	    ax2.plot(GregDay_all[index_not_nan_diff_pre_obs_nan], diff_pre_obs[index_not_nan_diff_pre_obs_nan], color=[0.72, 0.72, 0.72],marker='.',markersize=7,label='Pred-Obser')
	    ax2.plot(GregDay_all, diff_pre_obs,marker='.',markersize=4,color='k')
	    ax2.plot(GregDay_all, diff_reflect2_trimmed,marker='.',markersize=4,color='m')#trimmed median
	    
	    ax2.plot(GregDay_all, diff_reflect2_nan,marker='.',markersize=4,color='c',label='Daily diff.') #daily differences
	    
  
	    ax2.plot(days_group,diff_pre_obs_group, color='b',marker='o',markersize=4)
	    ax2.plot(days_group,daily_diff_group, color='r',marker='o',markersize=4)
	    
  
	    minobserv_2=numpy.nanmin(daily_diff_group)
	    
	    
	      
	    #ticks
	    plt.xticks(x_ticks_greg,Months)	  
  
	    ax2.grid(True)
	    ax2.set_xlabel('Julian dates')
	    ax2.set_ylabel('Daily  differences')	  
	    ax2.legend(loc='upper right',prop={'size':8})	     
		
	       
	    #plt.show()
	    plt.savefig('New_Plot_tiles_wings_vs4/VS4_1/Vs4_Four_variables_LL' + str(numlinhas) + '_'+ 'CC' + str(numcolunas) + '_row_' + str(i) + '_col_' + str(j) + '_' + str(year) + '_' + str(window)+ '_days' + '_' + str(num_valids_inside_window)+ '_val' + '.pdf',dpi=500)
	      
	    
	  
	  else:
	   
	      
	      
	      
	  
	  
	  ############################################################AXE1#####################################################################
	    ax1=plt.subplot(3,1,1)###################################################
	   
	    
	    ax1.plot(val_dias, val_nir, color=[0.72, 0.72, 0.72],linewidth=0.5,marker='.')  #original data
	    
	    
	    ax1.plot(val_dias, novoy3, color='r',linewidth=1,marker='.',markersize=5,label='pass_2')  #just to have all the information behind
	    ax1.plot(val_dias, novoy3_1, color='b',linewidth=1,marker='.',markersize=5,label='pass_3')  #just to have all the information behind
	    ax1.plot(val_dias, novoy3_2, color='g',linewidth=1,marker='.',markersize=5,label='pass_4')  #just to have all the information behind
	    ax1.plot(val_dias, novoy3_3, color='y',linewidth=2,marker='.',markersize=8,label='pass_5')  #just to have all the information behind	    
	    
	    
	    
	    
	    #ax1.plot(val_dias, novoy3, color=[0.72, 0.72, 0.72],linewidth=2,marker='.',markersize=8)  #just to have all the information behind
	    #ax1.plot(val_dias, novoy3_3, color='y',linewidth=2,marker='.',markersize=8)  #just to have all the information behind
	    
	    #ax1.plot(GregDay_all[index_not_nan_val_nir_nan], val_nir_nan[index_not_nan_val_nir_nan], color='k',marker='.',markersize=4,label='original')
	    #ax1.plot(GregDay_all[index_not_nan_final_predicted_nan], final_predicted[index_not_nan_final_predicted_nan], color='r',marker='.',markersize=4,label='predicted')
	    #ax1.plot(GregDay_all,novoy_median, color='m',marker='.',markersize=4,label='median')
	     
  
  
	    ax1.plot(GregDay_all,val_nir_nan, color='k',marker='.',markersize=5,label='observed')
	    
	    ax1.plot(GregDay_all,final_predicted, color='r',marker='.',markersize=4,label='predicted')
	    
	    
	    #to plot the region between reflectance and predicted curve for all the areas
	    for pp in range (0,days_all.shape[0]):
	      
	      ax1.fill_between(days_all[pp],refl_all[pp],predicted_all[pp], color='c', alpha=0.1)
	      
	     # to plot the selected area
	     
	    ax1.fill_between(days_group_area_selected,refl_area_selected,predicted_area_selected, color='y', alpha=0.5) 
	     
	    
		     
	    ax1.plot(days_group,refl_values_group, color='b',marker='o',markersize=4)
	    
	   
	    
	    
	    
	  
	    #set limits
	    plt.xlim([x_ticks_greg[0], x_ticks_greg[-1]+30])
	    minobserv=val_nir_nan[index_not_nan_val_nir_nan].min()
	    #ax1.set_ylim([minobserv-11, val_nir_nan[index_not_nan_val_nir_nan].max()])
	    
	     
		    
	    
	    #ticks
	    plt.xticks(x_ticks_greg,Months)
	
	
	    #put lines limiting the year
	    ax1.plot((x_ticks_greg[1],x_ticks_greg[1]),(minobserv-11,val_nir_nan[index_not_nan_val_nir_nan].max()) , color='b',label='limit year')#january year
	    
	    
	    ax1.plot((x_ticks_greg[-1],x_ticks_greg[-1]),(minobserv-11,val_nir_nan[index_not_nan_val_nir_nan].max()) , color='b') #january following year
	  
	    #lines
	    ax1.plot((GregDay_all[0], GregDay_all[-1]), (30, 30), color='r', linestyle='dashdot',label='Limit 30')
	    ax1.plot((GregDay_all[0], GregDay_all[-1]), (value_ref, value_ref), color='g', linestyle='dashdot',label='Limit' + str(percent_parameter) + '%')
	    plt.axvline(date,color='r',label='Selected ')
	    
	    
	    #candidates
	    for z in range(0,days_group_subgroup_candidates.shape[0]):
	      plt.axvline(days_group_subgroup_candidates[z],color='b',linestyle='dashdot')
	   
	    #days_group_subgroup_candidates,
	   # refl_values_group_subgroup_candidates,
	    #diff_pre_obs_group_sub_group_candidates,
	    #refl_values_group_subgroup_candidates,
	    
	    
	    
  
	    ax1.grid(True)
	    
	    ax1.set_xlabel('Julian dates')
	    ax1.set_ylabel('NIR reflectance')	  
	    ax1.legend(loc='upper right',prop={'size':8})
  
  
	    '''
	    #put text
	    ax1.text(date+2, minobserv-4, str(round(reflectance,2)), fontsize=8, color='b')  #reflectance value
	    ax1.text(date+2, minobserv-6, str(round(pred_minus_observ,2)), fontsize=8, color='c')#Predicted-observed
	    ax1.text(date+2, minobserv-8, str(round(percentage_minimum,2)), fontsize=8, color='k')#Percentage below the reflectance selected
	    ax1.text(date+2, minobserv-10, str(int(DT.datetime.strftime(num2date(date), "%j" ))), fontsize=8, color='m')#Date 
	    
	   #ax1.text(date+2, minobserv-12, str(round(density_alldata,2)), fontsize=8, color='r')#Density of the time series
	    ax1.text(date+2, minobserv-12, str(round(score_date,3)), fontsize=8, color='r')#Score of the area
	    ax1.text(date+2, minobserv-14, str(round(daily_diff_selected,2)), fontsize=8, color='k')#Daily difference
	    
	    
	    
	    ax1.text(GregDay_all[0], minobserv-4, 'Reflectance', fontsize=8, color='b')
	    ax1.text(GregDay_all[0], minobserv-6, 'Pre-Obs', fontsize=8, color='c')
	    ax1.text(GregDay_all[0], minobserv-8, '%below', fontsize=8, color='k')	      	 
	    ax1.text(GregDay_all[0], minobserv-10, 'Date', fontsize=8, color='m')   
   
	    #ax1.text(GregDay_all[0], minobserv-12, 'Time series dens.', fontsize=8, color='r')  
	    ax1.text(GregDay_all[0], minobserv-12, 'Score date.', fontsize=8, color='r')
	    ax1.text(GregDay_all[0], minobserv-14, 'Daily diff.', fontsize=8, color='k') 
	  '''
	    plt.title('Vs4_Four_variables_PAR -: Tile: LL' + str(numlinhas) + '_' +'CC' + str(numcolunas) + '_row: ' +str(i) + ', col:' + str(j)+ ' ' + str(window) + 'days window'+ '_' + str(num_valids_inside_window) + 'val_window' + ' ' + 'Minimum_values=' + str(group))    
	    
	    
	    
	    ############################################################AXE2#####################################################################
	    ax2=plt.subplot(3,1,2,sharex=ax1)###################################################
	   
	    
	    ax2.plot(val_dias, val_nir, color=[0.72, 0.72, 0.72],linewidth=0.5,marker='.')  #original data
	    
	    
	    #ax2.plot(val_dias, novoy3, color='m',linewidth=2,marker='.',markersize=8,label='pass_2')  #just to have all the information behind
	    #ax2.plot(val_dias, novoy3_1, color='b',linewidth=2,marker='.',markersize=8,label='pass_3')  #just to have all the information behind
	    #ax2.plot(val_dias, novoy3_2, color='g',linewidth=2,marker='.',markersize=8,label='pass_4')  #just to have all the information behind
	    #ax1.plot(val_dias, novoy3_3, color='y',linewidth=2,marker='.',markersize=8,label='pass_3')  #just to have all the information behind	    
	    
	    
	    
	    
	    #ax2.plot(val_dias, novoy3, color=[0.72, 0.72, 0.72],linewidth=2,marker='.',markersize=8)  #just to have all the information behind
	    ax2.plot(val_dias, novoy3_3, color=[0.72, 0.72, 0.72],linewidth=2,marker='.',markersize=8)  #just to have all the information behind
	    
	    #ax1.plot(GregDay_all[index_not_nan_val_nir_nan], val_nir_nan[index_not_nan_val_nir_nan], color='k',marker='.',markersize=4,label='original')
	    #ax1.plot(GregDay_all[index_not_nan_final_predicted_nan], final_predicted[index_not_nan_final_predicted_nan], color='r',marker='.',markersize=4,label='predicted')
	    #ax1.plot(GregDay_all,novoy_median, color='m',marker='.',markersize=4,label='median')
	     
  
  
	    ax2.plot(GregDay_all,val_nir_nan, color='k',marker='.',markersize=5,label='observed')
	    
	    ax2.plot(GregDay_all,final_predicted, color='r',marker='.',markersize=4,label='predicted')
	    
	    
	    #to plot the region between reflectance and predicted curve for all the areas
	    for pp in range (0,days_all.shape[0]):
	      
	      ax2.fill_between(days_all[pp],refl_all[pp],predicted_all[pp], color='c', alpha=0.1)
	      
	     # to plot the selected area
	     
	    ax2.fill_between(days_group_area_selected,refl_area_selected,predicted_area_selected, color='y', alpha=0.5) 
	     
	    
		     
	    ax2.plot(days_group,refl_values_group, color='b',marker='o',markersize=4)
	    
	   
	    
	    
	    
	  
	    #set limits
	    plt.xlim([x_ticks_greg[0], x_ticks_greg[-1]+30])
	    minobserv=val_nir_nan[index_not_nan_val_nir_nan].min()
	    #ax1.set_ylim([minobserv-11, val_nir_nan[index_not_nan_val_nir_nan].max()])
	    
	     
		    
	    
	    #ticks
	    plt.xticks(x_ticks_greg,Months)
	
	
	    #put lines limiting the year
	    ax2.plot((x_ticks_greg[1],x_ticks_greg[1]),(minobserv-11,val_nir_nan[index_not_nan_val_nir_nan].max()) , color='b',label='limit year')#january year
	    
	    
	    ax2.plot((x_ticks_greg[-1],x_ticks_greg[-1]),(minobserv-11,val_nir_nan[index_not_nan_val_nir_nan].max()) , color='b') #january following year
	  
	    #lines
	    ax2.plot((GregDay_all[0], GregDay_all[-1]), (30, 30), color='r', linestyle='dashdot',label='Limit 30')
	    ax2.plot((GregDay_all[0], GregDay_all[-1]), (value_ref, value_ref), color='g', linestyle='dashdot',label='Limit' + str(percent_parameter) + '%')
	    plt.axvline(date,color='r',label='Selected ')
	    
	    
	    #candidates
	    for z in range(0,days_group_subgroup_candidates.shape[0]):
	      plt.axvline(days_group_subgroup_candidates[z],color='b',linestyle='dashdot')
	   
	    #days_group_subgroup_candidates,
	   # refl_values_group_subgroup_candidates,
	    #diff_pre_obs_group_sub_group_candidates,
	    #refl_values_group_subgroup_candidates,
	    
	    
	    
  
	    ax2.grid(True)
	    
	    ax2.set_xlabel('Julian dates')
	    ax2.set_ylabel('NIR reflectance')	  
	    ax2.legend(loc='upper right',prop={'size':8})
  
	    #put text
	    ax2.text(date+2, minobserv-4, str(round(reflectance,2)), fontsize=8, color='b')  #reflectance value
	    ax2.text(date+2, minobserv-6, str(round(pred_minus_observ,2)), fontsize=8, color='c')#Predicted-observed
	    ax2.text(date+2, minobserv-8, str(round(percentage_minimum,2)), fontsize=8, color='k')#Percentage below the reflectance selected
	    ax2.text(date+2, minobserv-10, str(int(DT.datetime.strftime(num2date(date), "%j" ))), fontsize=8, color='m')#Date 
	    
	   #ax1.text(date+2, minobserv-12, str(round(density_alldata,2)), fontsize=8, color='r')#Density of the time series
	    ax2.text(date+2, minobserv-12, str(round(score_date,3)), fontsize=8, color='r')#Score of the area
	    ax2.text(date+2, minobserv-14, str(round(daily_diff_selected,2)), fontsize=8, color='k')#Daily difference
	    
	    
	    
	    ax2.text(GregDay_all[0], minobserv-4, 'Reflectance', fontsize=8, color='b')
	    ax2.text(GregDay_all[0], minobserv-6, 'Pre-Obs', fontsize=8, color='c')
	    ax2.text(GregDay_all[0], minobserv-8, '%below', fontsize=8, color='k')	      	 
	    ax2.text(GregDay_all[0], minobserv-10, 'Date', fontsize=8, color='m')   
   
	    #ax1.text(GregDay_all[0], minobserv-12, 'Time series dens.', fontsize=8, color='r')  
	    ax2.text(GregDay_all[0], minobserv-12, 'Score date.', fontsize=8, color='r')
	    ax2.text(GregDay_all[0], minobserv-14, 'Daily diff.', fontsize=8, color='k') 	    
	    
	    
	    
	    
	    '''
	    ax2=plt.subplot(2,1,2,sharex=ax1)    #link axes
	    ax2.plot(GregDay_all, diff_pre_obs,marker='.',markersize=4,color='k')  
	    ax2.plot(GregDay_all[index_not_nan_diff_pre_obs_nan], diff_pre_obs[index_not_nan_diff_pre_obs_nan], color=[0.72, 0.72, 0.72],marker='.',markersize=7,label='Pred-Obser')
	    ax2.plot(GregDay_all, diff_pre_obs,marker='.',markersize=4,color='k')
	    ax2.plot(GregDay_all, diff_reflect2_trimmed,marker='.',markersize=4,color='m')#trimmed median
	    
	    ax2.plot(GregDay_all, diff_reflect2_nan,marker='.',markersize=4,color='c',label='Daily diff.') #daily differences
	    
  
	    ax2.plot(days_group,diff_pre_obs_group, color='b',marker='o',markersize=4)
	    ax2.plot(days_group,daily_diff_group, color='r',marker='o',markersize=4)
	    
  
	    minobserv_2=numpy.nanmin(daily_diff_group)
	    ax2.text(date+2, 2.5, str(temporal_unc), fontsize=8, color='k')#Temporal uncertainty
	    ax2.text(GregDay_all[0], 2.5, 'Temporal unc.', fontsize=8, color='k')
	    '''
      ######################################AXE 3###################################################
	    ax3=plt.subplot(3,1,3,sharex=ax1)    #link axes
	    ax3.plot(GregDay_all, diff_pre_obs,marker='.',markersize=4,color='k')  
	    ax3.plot(GregDay_all[index_not_nan_diff_pre_obs_nan], diff_pre_obs[index_not_nan_diff_pre_obs_nan], color=[0.72, 0.72, 0.72],marker='.',markersize=7,label='Pred-Obser')
	    ax3.plot(GregDay_all, diff_pre_obs,marker='.',markersize=4,color='k')
	    ax3.plot(GregDay_all, diff_reflect2_trimmed,marker='.',markersize=4,color='m')#trimmed median
	    
	    ax3.plot(GregDay_all, diff_reflect2_nan,marker='.',markersize=4,color='c',label='Daily diff.') #daily differences
	    
      
	    ax3.plot(days_group,diff_pre_obs_group, color='b',marker='o',markersize=4)
	    ax3.plot(days_group,daily_diff_group, color='r',marker='o',markersize=4)
	    
      
	    minobserv_2=numpy.nanmin(daily_diff_group)
	    ax3.text(date+2, 2.5, str(temporal_unc), fontsize=8, color='k')#Temporal uncertainty
	    ax3.text(GregDay_all[0], 2.5, 'Temporal unc.', fontsize=8, color='k')	    
	    
	    
	    
	    
	    
	    
	    
	   
	    plt.axvline(date,color='r',label='Selected ')
	    
	    ##candidates
	    for z in range(0,days_group_subgroup_candidates.shape[0]):
	      plt.axvline(days_group_subgroup_candidates[z],color='b',linestyle='dashdot')	  
	    #set limits
	    plt.xlim([x_ticks_greg[0], x_ticks_greg[-1]+30])
	    minobserv=diff_pre_obs[index_not_nan_diff_pre_obs_nan].min()
	    #ax2.set_ylim([minobserv-11, diff_pre_obs[index_not_nan_diff_pre_obs_nan].max()])	  
	    #ticks
	    plt.xticks(x_ticks_greg,Months)	  
  
  
	    '''
	    ax2.grid(True)
	    ax2.set_xlabel('Julian dates')
	    ax2.set_ylabel('Daily  differences')	  
	    ax2.legend(loc='upper right',prop={'size':8})	
	    '''
	    
	    ax3.grid(True)
	    ax3.set_xlabel('Julian dates')
	    ax3.set_ylabel('Daily  differences')	  
	    ax3.legend(loc='upper right',prop={'size':8})	    
		
	       
	    #plt.show()
	    plt.savefig('New_Plot_tiles_wings_vs4/VS4_1/Vs4_3_par_Four_variables_LL' + str(numlinhas) + '_'+ 'CC' + str(numcolunas) + '_row_' + str(i) + '_col_' + str(j) + '_' + str(year) + '_' + str(window)+ '_days' + '_' + str(num_valids_inside_window)+ '_val_' +'min_' + str(group) + '.pdf',dpi=500)
	      
