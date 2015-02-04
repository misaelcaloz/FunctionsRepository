						#
						#
						#
		#################################################################################
		#																				#
		#		python functions are defined in this file								#
		#																				#
		#################################################################################
						#
						#
						#
						#
# --------packages 

import sys
sys.path.append("/usr/local/lib/root")
import ROOT
from ROOT import gROOT, TCanvas, TF1, TGraph, TLegend, TMath, TMultiGraph, TFile, TH2F, TH1F, TDirectory,TH2D, TH1D, TLatex, gPad, TH2I, TLine, TAxis, kTRUE, gStyle
from array import array
import numpy
import math
import pylab
import os
from decimal import Decimal






# -------------- THRESHOLD SCAN ------------------------------
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# -------------- Functions to analyse threshold scan .txt file 


def AnalyseThresholdScan(FileThresholdScan):

	#-------- conversion factor injection -> electrons

	injToElectrons=1660./.39	# WARNING: different depending the version of the chip: v2: 1660/0.39 v4:1660/0.25

	#------------ create canvas for Scurves

#	CanvasScurves = TCanvas("S-curves")
#	CanvasScurves.SetGrid()
#	CanvasScurves.SetFillColor(0)
#	CanvasScurves.cd()

	TDAC_value_dic = {}
	Threshold_value_dic = {}
	Sigma_value_dic = {}
	Chi2_value_dic = {}
	Scurve_plot_dic = {}
	Data_pointsX_dic = {}
	Data_pointsY_dic = {}
	
	lineNum=0
	grList = []
	TDAC_list=[]

	DataThresholdScan = open(FileThresholdScan)

	for line in DataThresholdScan:
		n=line.split()
		lineNum +=1
		TDAC_value=int(n[3])	
		TDAC_value_dic["r"+n[1]+"_c"+n[2]+""] = TDAC_value
				
		#------------ create odd and even list
		#------------ in order to select x and y values from data file
		
		lenght_nn=len(n)
		lenght_n=lenght_nn/2.
		seq = range(int(lenght_n))
		def add(x, y): 
			return x+y
		list_pair=map(add, seq, seq)
		list_pair.append(lenght_nn-2)
		list_pair.reverse()
		list_pair.pop()
		list_pair.pop()
		list_pair.pop()
		list_pair.reverse()
		list_impair=[]
		for i in list_pair:
			list_impair.append(i+1)
		
		#---------- create list of 'real' values
		
		x1=[]
		for r in list_pair:
			x1.append(float(n[r])*injToElectrons)
		y1=[]
		for r in list_impair:
			y1.append(float(n[r]))
		
		Data_pointsX_dic["r"+n[1]+"_c"+n[2]+""] = x1
		Data_pointsY_dic["r"+n[1]+"_c"+n[2]+""] = y1
		#---------- necessary to plot x and y
		
		x=array("d",x1)
		y=array("d",y1)
		
		#------------ create fitting function ERF

		myfit_ErrFunc=TF1("myfit_ErrFunc", "[0]+0.5*TMath::Erf([1]*([2]+x))",0.,3000.)
		
		#---------- create the plot
		
		g=TGraph(len(x),x,y)
		g.GetXaxis().SetTitle("injection [e]")
		g.GetYaxis().SetTitle("Probe")
		g.SetFillColor(33)
		g.SetMarkerColor(1)
		g.SetMarkerSize(1)
		g.SetLineColor(1)
		g.SetLineWidth(1)
		
		#----- set fitting parameters
		
		myfit_ErrFunc.SetParameter(0,0.5);
		myfit_ErrFunc.SetParameter(1,20/injToElectrons);
		myfit_ErrFunc.SetParameter(2,-0.02*injToElectrons)
		
		#------ fit
		
		g.Fit("myfit_ErrFunc")
		
		#--------- keep all s-curves into the same canvas 
		
#		if lineNum ==1:	
#			g.Draw("AC")
#		else:
#			g.Draw("C")		

		
		
		
		#-------- get the threshold
		
		thresh=-1*myfit_ErrFunc.GetParameter(2)
		sigma=math.sqrt(2.)/(myfit_ErrFunc.GetParameter(1)*2)
		thresh_int=int(thresh)
		sigma_int=int(sigma)
		
		Threshold_value_dic["r"+n[1]+"_c"+n[2]+""] = thresh
		Sigma_value_dic["r"+n[1]+"_c"+n[2]+""] = sigma
		Chi2_value_dic["r"+n[1]+"_c"+n[2]+""] = myfit_ErrFunc.GetChisquare()
		
		#--------- Set titles
		
		g.SetNameTitle("r"+n[1]+"_c"+n[2]+"_thr="+str(thresh_int)+"_TDAC="+str(int(n[3]))+" ",n[2]+"row_"+n[1]+"col"+n[2]+"_thr"+str(thresh_int)+"_TDAC."+str(int(n[3]))+"")
		
		g.SetTitle("r"+n[1]+"_c"+n[2]+" , thr = "+str(thresh_int)+" , TDAC = "+str(int(n[3]))+" , sigma = "+str(sigma_int)+" , Chi2 = "+str(myfit_ErrFunc.GetChisquare())+"")
		
		Scurve_plot_dic["r"+n[1]+"_c"+n[2]+""] = g
		grList.append(g)
		

	return [TDAC_value_dic, Threshold_value_dic, Sigma_value_dic, Chi2_value_dic, Scurve_plot_dic, Data_pointsX_dic, Data_pointsY_dic]
	
	
	
# --------------------------------------------------------------------------------------------------

# END -----------  THRESHOLD SCAN ---------------------------------------------------------------- END
























# -------------- ANALYSE TDAC THRESHOLD SCANS ------------
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# -------------- Functions to analyse several threshold scans with a different TDAC
#				and predict a TDAC for each pixel using the FITTING (not interpolation)

# FileTDAC_list is a list containing threshold scan file name
# TDAC_value_list is a list containing TDAC value for each threshold scans (in the same order) 



def AnalyseTDACThresholdScans(FileTDAC_list,TDAC_value_list, TargetThreshold):
	
	g_Thr_vs_TDAC_List = []

	
	#-------- conversion factor injection -> electrons

	injToElectrons=1660./.39			# WARNING: different depending the version of the chip: v2: 1660/0.39 v4: 1660/0.25



	AllThresh_dic = {}
	AllSigma_dic = {}
	AllChi2_dic = {}
	AllData_pointsX_dic = {}
	AllData_pointsY_dic = {}
	AllSuggTDAC_dic = {}
	AllExpectThresh_dic = {}
	AllThr_vs_TDAC_graphs_dic = {}
	DIC_AnalyseThresholdScan = {}
	
	for n, TDACcnt in zip(range(2,8), TDAC_value_list):
		FileDataThrScan = FileTDAC_list[n-2]
			
		# -------- Reset dictionaries and lists

		TDAC_value_dic = {}
		Threshold_value_dic = {}
		Sigma_value_dic = {}
		Chi2_value_dic = {}
		Scurve_plot_dic = {}
		Data_pointsX_dic = {}
		Data_pointsY_dic = {}
		graphList = []


		# -------- Fill dictionaries with function AnalyseThresholdScan()

		DIC_AnalyseThresholdScan = AnalyseThresholdScan(FileDataThrScan)

		TDAC_value_dic = DIC_AnalyseThresholdScan[0]
		Threshold_value_dic = DIC_AnalyseThresholdScan[1]
		Sigma_value_dic = DIC_AnalyseThresholdScan[2]
		Chi2_value_dic = DIC_AnalyseThresholdScan[3]
		Scurve_plot_dic = DIC_AnalyseThresholdScan[4]
		Data_pointsX_dic = DIC_AnalyseThresholdScan[5]
		Data_pointsY_dic = DIC_AnalyseThresholdScan[6]



		for r in range(24):
			for c in range(12,48):
				TDAC_value = TDAC_value_dic["r"+str(r)+"_c"+str(c)+""]
				Threshold_value = Threshold_value_dic["r"+str(r)+"_c"+str(c)+""]
				Sigma_value = Sigma_value_dic["r"+str(r)+"_c"+str(c)+""]
				Chi2_value = Chi2_value_dic["r"+str(r)+"_c"+str(c)+""]
				Data_pointsX_list = Data_pointsX_dic["r"+str(r)+"_c"+str(c)+""]
				Data_pointsY_list = Data_pointsY_dic["r"+str(r)+"_c"+str(c)+""]


				AllThresh_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Threshold_value
				AllSigma_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Sigma_value
				AllChi2_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Chi2_value
				AllData_pointsX_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Data_pointsX_list
				AllData_pointsY_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Data_pointsY_list
		
	# ------- Find suggested TDAC and expected threshold with sugg TDAC	
		
	for r in range(24):
		for c in range(12,48):
		
			x = []
			y = []			
									
			for TDACcnt in TDAC_value_list:				
				if AllThresh_dic["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] > 300:				
					x.append(TDACcnt)
					y.append(AllThresh_dic["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""])
					
					
			# - necessary to plot x and y
		
			x_array=array("d",x)
			y_array=array("d",y)
				
			if len(x) > 0: 	
				g_Thr_vs_TDAC = TGraph(len(x_array),x_array,y_array)
			else: 
				x = []
				y = []
				
				for TDACcnt in TDAC_value_list:				
					x.append(TDACcnt)
					y.append(AllThresh_dic["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""])
				
				# - necessary to plot x and y
				
				x_array=array("d",x)
				y_array=array("d",y)
				
				g_Thr_vs_TDAC = TGraph(len(x_array),x_array,y_array)


			myfit_Thr_vs_TDAC=TF1("myfit_Thr_vs_TDAC", "pol1")
			g_Thr_vs_TDAC.Fit("myfit_Thr_vs_TDAC")

			p0 = myfit_Thr_vs_TDAC.GetParameter(0)
			p1 = myfit_Thr_vs_TDAC.GetParameter(1)			
			
			if p1 != 0:
				suggTDAC = (TargetThreshold-p0)/p1
			else:
				suggTDAC = 0

			suggTDAC_arrondi = int(Decimal(str(round(suggTDAC,0))))	# arrondi 
			
			if suggTDAC_arrondi < min(x):
				suggTDAC_arrondi = min(x)
			if suggTDAC_arrondi > 15:
				suggTDAC_arrondi = 15
			
			expected_Threshold = p1*suggTDAC_arrondi+p0



			# - Customize graph
			
			g_Thr_vs_TDAC.GetXaxis().SetTitle("TDAC")
			g_Thr_vs_TDAC.GetYaxis().SetTitle("Threshold [e]")
			Xaxis = g_Thr_vs_TDAC.GetXaxis()
			Xaxis.SetLimits(0.,16.)
			g_Thr_vs_TDAC.GetHistogram().SetMaximum(2000.)
			g_Thr_vs_TDAC.GetHistogram().SetMinimum(0.)


			AllExpectThresh_dic ["r"+str(r)+"_c"+str(c)+""] = expected_Threshold
			AllSuggTDAC_dic ["r"+str(r)+"_c"+str(c)+""] = suggTDAC_arrondi
			
			
			g_Thr_vs_TDAC.SetNameTitle("r"+str(r)+"_c"+str(c)+"_suggTDAC "+str(suggTDAC_arrondi)+"","r"+str(r)+"_c"+str(c)+" ")
			g_Thr_vs_TDAC.SetTitle("Thr_vs_TDAC for  r"+str(r)+"_c"+str(c)+",  suggTDAC "+str(suggTDAC_arrondi)+", expected threshold = "+str(int(expected_Threshold))+"")
		

			AllThr_vs_TDAC_graphs_dic ["r"+str(r)+"_c"+str(c)+""] = g_Thr_vs_TDAC



			g_Thr_vs_TDAC_List.append(g_Thr_vs_TDAC)




				
	return [AllThresh_dic, AllSigma_dic, AllChi2_dic, AllData_pointsX_dic, AllData_pointsY_dic, AllSuggTDAC_dic, AllExpectThresh_dic, AllThr_vs_TDAC_graphs_dic]



# --------------------------------------------------------------------------------------------------

# END -----------  ANLYSE TDAC THRESHOLD SCANS ---------------------------------------------------------------- END






















# -------------- ANALYSE TDAC THRESHOLD SCANS VERSION 2 ------------
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# -------------- Functions to analyse several threshold scans with a different TDAC
#				and predict a TDAC for each pixel by INTERPOLATION

# FileTDAC_list is a list containing threshold scan file name
# TDAC_value_list is a list containing TDAC value for each threshold scans (in the same order)
# TargetThreshold is the target threshold in [e]


def AnalyseTDACThresholdScansV2(FileTDAC_list,TDAC_value_list, TargetThreshold):
			
	#-------- conversion factor injection -> electrons

	injToElectrons=1660./.39			# WARNING: different depending the version of the chip: v2: 1660/0.39 v4: 1660/0.25

	# --- Initialize dictionaries and lists just to be sure

	AllThresh_dic = {}
	AllSigma_dic = {}
	AllChi2_dic = {}
	AllData_pointsX_dic = {}
	AllData_pointsY_dic = {}
	AllSuggTDAC_dic = {}
	AllExpectThresh_dic = {}
	AllThr_vs_TDAC_graphs_dic = {}
	AllThresh_simpl_dic = {}
	DIC_AnalyseThresholdScan = {}
	g_Thr_vs_TDAC_List = []
	
	for n, TDACcnt in zip(range(2,8), TDAC_value_list):
		FileDataThrScan = FileTDAC_list[n-2]
			
		# -------- Reset dictionaries and lists

		TDAC_value_dic = {}
		Threshold_value_dic = {}
		Sigma_value_dic = {}
		Chi2_value_dic = {}
		Scurve_plot_dic = {}
		Data_pointsX_dic = {}
		Data_pointsY_dic = {}
		graphList = []

		# -------- Fill dictionaries with function AnalyseThresholdScan()

		DIC_AnalyseThresholdScan = AnalyseThresholdScan(FileDataThrScan)

		TDAC_value_dic = DIC_AnalyseThresholdScan[0]
		Threshold_value_dic = DIC_AnalyseThresholdScan[1]
		Sigma_value_dic = DIC_AnalyseThresholdScan[2]
		Chi2_value_dic = DIC_AnalyseThresholdScan[3]
		Scurve_plot_dic = DIC_AnalyseThresholdScan[4]
		Data_pointsX_dic = DIC_AnalyseThresholdScan[5]
		Data_pointsY_dic = DIC_AnalyseThresholdScan[6]

		for r in range(24):
			for c in range(12,48):
				TDAC_value = TDAC_value_dic["r"+str(r)+"_c"+str(c)+""]
				Threshold_value = Threshold_value_dic["r"+str(r)+"_c"+str(c)+""]
				Sigma_value = Sigma_value_dic["r"+str(r)+"_c"+str(c)+""]
				Chi2_value = Chi2_value_dic["r"+str(r)+"_c"+str(c)+""]
				Data_pointsX_list = Data_pointsX_dic["r"+str(r)+"_c"+str(c)+""]
				Data_pointsY_list = Data_pointsY_dic["r"+str(r)+"_c"+str(c)+""]

				AllThresh_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Threshold_value
				AllSigma_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Sigma_value
				AllChi2_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Chi2_value
				AllData_pointsX_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Data_pointsX_list
				AllData_pointsY_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Data_pointsY_list
		
	# ------- Find suggested TDAC and expected threshold with sugg TDAC	
		
	for r in range(24):
		for c in range(12,48):
		
			x = []
			y = []
			Diff_Abs_list = []
			Diff_list = []
			Thresh_list = []
			TDAC_list = []
			Diff_Abs_dic = {}			
									
			for TDACcnt in TDAC_value_list:				
				
				Thresh_value = AllThresh_dic["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""]
				Chi2_val = AllChi2_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""]
				
				if len(AllData_pointsY_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""]) > 2:
					First_Y_point_scurve = AllData_pointsY_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""][0]
					Second_Y_point_scurve = AllData_pointsY_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""][1]
					Third_Y_point_scurve = AllData_pointsY_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""][2]

				else:
					First_Y_point_scurve = 2
					Second_Y_point_scurve = 2
					Third_Y_point_scurve = 2				

				# --- Selection of good S-curves: good Chi2 and 3 first points of the S-curve under 
				# ||||||||||||||||||||||||   WARNING VERY IMPORTANT SELECTIONS ||||||||||||||||||||||||||||||
				# TODO : verify selections are OK 
				
				if Chi2_val > 0.01 or Chi2_val == 0:
					Thresh_value = 0
				if First_Y_point_scurve > 0.1 and Second_Y_point_scurve > 0.1 and Third_Y_point_scurve > 0.1:
					Thresh_value = 0
				if Thresh_value < 0:
					Thresh_value = 0
													
				# --- Create 2 lists of data points and compute the difference between target threshold and data points
				
				if Thresh_value > 0:

					Diff_Abs_list.append(abs(TargetThreshold - Thresh_value))
					Diff_list.append(TargetThreshold - Thresh_value)
					Thresh_list.append(Thresh_value)
					TDAC_list.append(TDACcnt)
					
				# --- Necessary to plot TGraph

				x.append(TDACcnt)
				y.append(Thresh_value)		
					
				AllThresh_simpl_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Thresh_value					
																	
			# --- If lengh of Diff
			
			if len(Diff_Abs_list) == 0:
				Diff_Abs_list.append(1)
				Diff_list.append(1)
				Thresh_list.append(0)
				TDAC_list.append(15)
							
															
			# --- Take the closer point to the target threshold	

			Min_Diff_Abs_list = min(Diff_Abs_list)
			Min_index = Diff_Abs_list.index(Min_Diff_Abs_list)
			
			first_point_TDAC = TDAC_list[Min_index]
			first_point_thresh = Thresh_list[Min_index]			
				
			# --- Define which second point has to be taken in order to interpolate
			if len(Diff_Abs_list) < 2:
				second_point_TDAC = -1										# to identify it easely later
				second_point_thresh = -1				
															
			elif Diff_list[Min_index] < 0 and Min_index > 0:					# first selection: first point
				second_point_TDAC = TDAC_list[Min_index-1]
				second_point_thresh = Thresh_list[Min_index-1]
			elif Diff_list[Min_index] > 0 and Min_index < 5:				# second selection: last point
				second_point_TDAC = TDAC_list[Min_index+1]
				second_point_thresh = Thresh_list[Min_index+1]
			else:															# else put a value of -1 in order
				second_point_TDAC = -1										# to identify it easely later
				second_point_thresh = -1
				
			if second_point_thresh < 300 or second_point_thresh > 2000:		# third selection: only good data points
				second_point_thresh = -1									# can be selected as a second points
							
			# --- Suggest a TDAC value 	
			
			#---
			if second_point_thresh == -1:							# if the second point has not be found 
				suggTDAC = first_point_TDAC							# or doesn't match the criteria
				suggTDAC_arrondi = suggTDAC							# suggested TDAC is the FIRST POINT (the point
																	# which is the closer to the target threshold)
				expected_Threshold = first_point_thresh
				p = -1
			
			#---
			elif second_point_TDAC != -1: 
				if first_point_TDAC - second_point_TDAC < 0:		# define first point and second point
					x1 = first_point_TDAC							# in order to calculate the coefficient p 
					x2 = second_point_TDAC							# in y = p * x + b
					y1 = first_point_thresh
					y2 = second_point_thresh
					orientation = 1
															
				else:												# orientation 2 is when the threshold of the 
					x2 = first_point_TDAC							# second point is higher than the threshold
					x1 = second_point_TDAC							# of the first point 
					y2 = first_point_thresh							# THIS SHOULD ALMOST NEVER HAPPEN NORMALLY
					y1 = second_point_thresh
					orientation = 2
				
				p = (y2-y1)/(x2-x1)								# calculate coeff p
				b = (y1 - p * x1)								# calculate coeff b
				suggTDAC = (TargetThreshold - b)/p				# suggest a TDAC 
				
				if p > 0:
					suggTDAC_arrondi = int(Decimal(str(round(suggTDAC,0))))	# arrondi 
					expected_Threshold = p * suggTDAC_arrondi + b
				else:
					if orientation == 1:											# everything is fine
						suggTDAC_arrondi = int(Decimal(str(round(suggTDAC,0))))		# define suggTDAC_arrondi
						expected_Threshold = p * suggTDAC_arrondi + b
						
					elif orientation == 2:											# bad case: go back to 
						suggTDAC_arrondi = first_point_TDAC							# the first point
						expected_Threshold = first_point_thresh
						
					else:								# else condition: IS NORMALLY IMPOSSIBLE
						suggTDAC_arrondi = -2			# -2 in order to identify easily
						expected_Threshold = -2			# THIS CONDITION SHOULD NOT HAPPEN
					
					
			#---
			else:										# else condition: IS NORMALLY IMPOSSIBLE
				suggTDAC = -2							# -2 in order to identify easily 
				expected_Threshold = -2					# THIS CONDITION SHOULD NOT HAPPEN
				p = -2

			
			# - necessary to plot x and y
		
			x_array=array("d",x)
			y_array=array("d",y)
								
			g_Thr_vs_TDAC = TGraph(len(x_array),x_array,y_array)

			# - Customize graph
			
			g_Thr_vs_TDAC.GetXaxis().SetTitle("TDAC")
			g_Thr_vs_TDAC.GetYaxis().SetTitle("Threshold [e]")
			Xaxis = g_Thr_vs_TDAC.GetXaxis()
			Xaxis.SetLimits(0.,16.)
			g_Thr_vs_TDAC.GetHistogram().SetMaximum(2000.)
			g_Thr_vs_TDAC.GetHistogram().SetMinimum(0.)

			AllExpectThresh_dic ["r"+str(r)+"_c"+str(c)+""] = expected_Threshold
			AllSuggTDAC_dic ["r"+str(r)+"_c"+str(c)+""] = suggTDAC_arrondi			
			
			g_Thr_vs_TDAC.SetNameTitle("r"+str(r)+"_c"+str(c)+"_suggTDAC "+str(suggTDAC_arrondi)+"","r"+str(r)+"_c"+str(c)+" ")
			g_Thr_vs_TDAC.SetTitle("r"+str(r)+"_c"+str(c)+",  suggTDAC "+str(suggTDAC_arrondi)+", exp th = "+str(int(expected_Threshold))+", first TDAC = "+str(first_point_TDAC)+" , second TDAC = "+str(second_point_TDAC)+" , p = "+str(p)+"")		

			AllThr_vs_TDAC_graphs_dic ["r"+str(r)+"_c"+str(c)+""] = g_Thr_vs_TDAC

			g_Thr_vs_TDAC_List.append(g_Thr_vs_TDAC)




				
	return [AllThresh_dic, AllSigma_dic, AllChi2_dic, AllData_pointsX_dic, AllData_pointsY_dic, AllSuggTDAC_dic, AllExpectThresh_dic, AllThr_vs_TDAC_graphs_dic, AllThresh_simpl_dic]



# --------------------------------------------------------------------------------------------------

# END -----------  ANLYSE TDAC THRESHOLD SCANS ---------------------------------------------------------------- END
































## -------------- SET TDAC VALUE ------------
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## -------------- Functions to set TDAC value into a .txt file
#
############
## WARNING # : to call this function file_TDAC_list must be defined before in the code
## ------- #
############
#
##file_TDAC = open("suggested_TDACfileMathieuTuning.txt", "w")  # "r" instead of "w" in order to take the original suggested_TDACfile as an input instead of create TDACfile with all 6
##
##file_TDAC_list = []
##lineNum_TDAC=0
##
##for counter in range(1,1441): # old way: create a TDAC all 6 file and modify it
##	if counter < 289:
##		file_TDAC_list.append(15)
##	elif 288 < counter < 1153:
##		file_TDAC_list.append(6)
##	elif counter > 1152:
##		file_TDAC_list.append(15)
#
#def setTDAC(CMOSrow,CMOScol,newTDAC_value):
#	cnt=0
#	for i in range(60) :
#		for j in range(24):
#			if j==CMOSrow and i==CMOScol:
#				file_TDAC_list[cnt]=newTDAC_value
#			cnt += 1
#
#
## --------------------------------------------------------------------------------------------------
#
## END -----------  SET TDAC VALUE ---------------------------------------------------------------- END







