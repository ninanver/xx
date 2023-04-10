#Script to compute the Discrete Cross-Correlation Function (DCF, Edelson & Krolik 1988) and/or the Interpolated Cross-Correlation Function (ICF, White & Peterson 1994).


#######################################
########### USER-INTERFACE ############
#######################################

#Input path of the first light curve (7mm)
input_path_first = '/enter/the/path/to/Mrk421_7mm.csv' #enter the path to Mrk421_7mm.csv or OQ334_7mm.csv

#Input path of the second light curve (14mm)
input_path_second = '/enter/the/path/to/Mrk421_14mm.csv' #enter the path to Mrk421_14mm.csv or OQ334_14mm.csv

#Output path of the folder where the results are stored
output_path = '/enter/the/path/to/folder/where/results/are/stored/' #enter the path to the folder where the results are stored

#Source name
source_name = 'Mrk421' #set to Mrk421 or OQ334, names the output files by source_name

#DCF 
comp_DCF = True #if comp_DCF == True, the Discrete Cross-Correlation Function will be computed
binsize =  #the bin size should be set to the average sampling time across BOTH light curves

#ICF
comp_ICF = True #if comp_ICF == True, the Interpolated Cross-Correlation Function will be computed
interp_unit =  #the interpolation unit should be set to 1/2 of the average sampling time across BOTH light curves
interp_first = True #if interp_first == True, the first light curve will be interpolated
interp_second = True #if interp_second == True, the second light curve will be interpolated



#######################################
########### PROGRAM-CORE ##############
#######################################

import numpy as np
import pandas as pd

#Function to compute the DCF
def dcf_funct(first, second, binsize, startfactor):

	#first: first light curve (LC) (as Pandas Dataframe with columns: mjd, flux density, uncertainty of flux density)
	#second: second LC (as Pandas Dataframe with columns: mjd, flux density, uncertainty of flux density)
	#startfactor: 1/3 of duration where LCs overlap in sampling / binsize

	first_mjd = first.iloc[:, 0].values
	first_flux = first.iloc[:, 1].values
	second_mjd = second.iloc[:, 0].values
	second_flux = second.iloc[:, 1].values

	lag_array = np.array([])
	dcf_array = np.array([])
	dcf_err_array = np.array([])

	#Number of different time lags
	interval = 2. * startfactor + 1.

	#Compute DCF
	for i in range(int(interval)):
		
		#Compute time lag
		if i > 0:
			lag = lag_array[i - 1] + binsize
		else:
			lag = -startfactor * binsize
		lag_array = np.append(lag_array, lag)

		#Compute correlation coefficient
		corr_flux1 = np.array([])
		corr_flux2 = np.array([])
		for j in range(len(first_flux)):
			flux2 = second_flux[np.where(second_mjd >= first_mjd[j] + lag - binsize / 2.)]
			mjd2 = second_mjd[np.where(second_mjd >= first_mjd[j] + lag - binsize / 2.)]
			flux2 = flux2[np.where(mjd2 < first_mjd[j] + lag + binsize / 2.)]
			
			corr_flux2 = np.append(corr_flux2, flux2)

			if len(flux2) > 0:
				corr_flux1 = np.append(corr_flux1, np.full_like(flux2, first_flux[j]))

		if len(corr_flux2) >= 2:		 
			udcf = (corr_flux1 - np.mean(corr_flux1)) * (corr_flux2 - np.mean(corr_flux2)) / np.sqrt(np.sum((corr_flux1 - np.mean(corr_flux1)) ** 2.) / len(corr_flux1) * np.sum((corr_flux2 - np.mean(corr_flux2)) ** 2.) / len(corr_flux2))
			dcf_value = np.mean(udcf)
			dcf_array = np.append(dcf_array, dcf_value)
			dcf_err = np.sqrt(np.sum((udcf - dcf_value) ** 2.)/((len(udcf) - 1) * (len(udcf))))
			dcf_err_array = np.append(dcf_err_array, dcf_err)
		else:
			dcf_array = np.append(dcf_array, np.nan)
			dcf_err_array = np.append(dcf_err_array, np.nan)

	dcf = pd.DataFrame()
	dcf["Lag"] = lag_array
	dcf["Correlation_Coefficient"] = dcf_array
	dcf["Error"] = dcf_err_array
	dcf.dropna(inplace = True)

	return dcf

#Function to compute the ICF
def icf_funct(first, second, interp_unit, startfactor, interp_first, interp_second):

	#first: first light curve (LC) (as Pandas Dataframe with columns: mjd, flux density, uncertainty of flux density)
	#second: second LC (as Pandas Dataframe with columns: mjd, flux density, uncertainty of flux density)
	#interp_unit: interpolation unit where LCs will be interpolated
	#startfactor: 1/3 of duration where LCs overlap in sampling / interpolation unit
	#interp_first: if True, the first LC will be interpolated
	#interp_second: if True, the second LC will be interpolated

	first_mjd = first.iloc[:, 0].values
	first_flux = first.iloc[:, 1].values
	second_mjd = second.iloc[:, 0].values
	second_flux = second.iloc[:, 1].values

	lag_array = np.array([])
	icf1_array = np.array([])
	icf2_array = np.array([])
	weight1_array = np.array([])
	weight2_array = np.array([])

	#Number of different time lags
	interval = 2. * startfactor + 1.

	#Compute ICF
	if interp_first == False and interp_second == False:
		print("No light curve is chosen to be interpolated!")
		print("Please choose which light curve should be interpolated by setting interp_first and/or interp_second to True!")
		icf = pd.DataFrame()
	else:
		for i in range(int(interval)):
		
			#Compute time lag
			if i > 0:
				lag = lag_array[i - 1] + interp_unit
			else:
				lag = -startfactor * interp_unit
			lag_array = np.append(lag_array, lag)

			#Compute correlation coefficient
			if interp_second == True:
				flux1 = first_flux[np.where(first_mjd <= np.max(second_mjd) - lag)]
				mjd1 = first_mjd[np.where(first_mjd <= np.max(second_mjd) - lag)]
				flux1 = flux1[np.where(mjd1 >= np.min(second_mjd) - lag)]
				mjd1 = mjd1[np.where(mjd1 >= np.min(second_mjd) - lag)]

				flux2_interp = np.interp(mjd1 + lag, second_mjd, second_flux)
			
				weight1_array = np.append(weight1_array, len(mjd1))

				icf1_value = np.mean((flux1 - np.mean(flux1)) * (flux2_interp - np.mean(flux2_interp)) / (np.std(flux1) * np.std(flux2_interp)))
				icf1_array = np.append(icf1_array, icf1_value)

			if interp_first == True:
				flux2 = second_flux[np.where(second_mjd <= np.max(first_mjd) + lag)]
				mjd2 = second_mjd[np.where(second_mjd <= np.max(first_mjd) + lag)]
				flux2 = flux2[np.where(mjd2 >= np.min(first_mjd) + lag)]
				mjd2 = mjd2[np.where(mjd2 >= np.min(first_mjd) + lag)]

				flux1_interp = np.interp(mjd2 - lag, first_mjd, first_flux)
			
				weight2_array = np.append(weight2_array, len(mjd2))

				icf2_value = np.mean((flux1_interp - np.mean(flux1_interp)) * (flux2 - np.mean(flux2)) / (np.std(flux1_interp) * np.std(flux2)))
				icf2_array = np.append(icf2_array, icf2_value)

		if interp_second == True and interp_first == False:
			icf_array = icf1_array

		elif interp_second == False and interp_first == True:
			icf_array = icf2_array

		elif interp_second == True and interp_first == True:
			icf_array = (icf1_array * weight1_array + icf2_array * weight2_array) / (weight1_array + weight2_array)

		icf = pd.DataFrame()
		icf["Lag"] = lag_array
		icf["Correlation_Coefficient"] = icf_array

	return icf



#Load first light curve (LC)
first = pd.read_csv(input_path_first)
first.astype(float)
first.columns = ['MJD', 'Flux Density [Jy]', 'Uncertainty [Jy]']
first.sort_values(by = ['MJD'], inplace = True)

#Load second LC
second = pd.read_csv(input_path_second)
second.astype(float)
second.columns = ['MJD', 'Flux Density [Jy]', 'Uncertainty [Jy]']
second.sort_values(by = ['MJD'], inplace = True)

#Compute duration where the LCs overlap in sampling
first_mjd = first.iloc[:, 0].values
second_mjd = second.iloc[:, 0].values

if np.min(first_mjd) < np.min(second_mjd):
    t_min = np.min(second_mjd)
else:
    t_min = np.min(first_mjd)

if np.max(first_mjd) < np.max(second_mjd):
    t_max = np.max(first_mjd)
else:
    t_max = np.max(second_mjd)

overlap_duration = t_max - t_min

#Compute DCF
if comp_DCF == True:
	print("Compute DCF")
	print(f"Binsize = {binsize:.3f} days")

	#Compute startfactor
	dcf_startfactor = np.around((1. / 3.) * overlap_duration / binsize)

	#Compute DCF
	DCF = dcf_funct(first = first, second = second, binsize = binsize, startfactor = dcf_startfactor)

	DCF.to_csv(output_path + source_name + "_dcf.csv", index = False)

	DCF_lag = DCF.iloc[:, 0].values
	DCF_coefficient = DCF.iloc[:, 1].values
	DCF_error = DCF.iloc[:, 2].values

	DCF_coefficient_max_total = np.max(DCF_coefficient)
	DCF_error_max_total = DCF_error[np.where(DCF_coefficient == DCF_coefficient_max_total)]
	DCF_lag_max_total = DCF_lag[np.where(DCF_coefficient == DCF_coefficient_max_total)]
	
	print(f"Peak correlation coefficient = {DCF_coefficient_max_total:.3f} +/- {DCF_error_max_total[0]:.3f}")
	print(f"Peak time lag = {DCF_lag_max_total[0]:.3f} days")

else:
	print("DCF is not computed!")

#Compute ICF
if comp_ICF == True:
	print("\nCompute ICF")
	print(f"Interpolation unit = {interp_unit:.3f} days\n")

	#Compute startfactor
	icf_startfactor = np.around((1. / 3.) * overlap_duration / interp_unit)

	#Compute ICF
	if interp_first == True and interp_second == False:
		print("Only first light curve is interpolated!\n")

	elif interp_first == False and interp_second == True:
		print("Only second light curve is interpolated!\n")

	elif interp_first == True and interp_second == True:
		print("Both light curve are interpolated!")
		print("The ICF is computed twice.")
		print("In the first pass, the second light curve is interpolated.")
		print("In the second pass, the first light curve is interpolated.")
		print("The results from the two passes are averaged.\n")


	ICF = icf_funct(first = first, second = second, interp_unit = interp_unit, startfactor = icf_startfactor, interp_first = interp_first, interp_second = interp_second)

	ICF_empty = ICF.empty
	if ICF_empty == False:
		ICF.to_csv(output_path + source_name + "_icf.csv", index = False)

		ICF_lag = ICF.iloc[:, 0].values
		ICF_coefficient = ICF.iloc[:, 1].values

		ICF_coefficient_max_total = np.max(ICF_coefficient)
		ICF_lag_max_total = ICF_lag[np.where(ICF_coefficient == ICF_coefficient_max_total)]

		print(f"Peak correlation coefficient = {ICF_coefficient_max_total:.3f}")
		print(f"Peak time lag = {ICF_lag_max_total[0]:.3f} days")

else:
	print("ICF is not computed!")


