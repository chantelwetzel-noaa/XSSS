#V3.30.11.00-safe;_2018_03_07;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_11.6
#This is a work of the U.S. Government and is not subject to copyright protection in the United States.
#Foreign copyrights may apply. See copyright.txt for more information.
#_user_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_user_info_available_at:https://vlab.ncep.noaa.gov/group/stock-synthesis
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Morph_between/within_stdev_ratio (no read if N_morphs=1)
#_Cond  1 #vector_Morphdist_(-1_in_first_val_gives_normal_approx)
#
2 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
0 #_Nblock_Patterns
# begin and end years of blocks
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#  autogen
1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
# 
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement 
#
0 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
  #_no additional input for selected M option; read 1P per morph
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K; 4=not implemented
1 #_Age(post-settlement)_for_L1;linear growth below this
30 #_Growth_Age_for_L2 (999 to use as Linf)
0.055 #_exponential decay for growth above maxage (fixed at 0.2 in 3.24; value should approx initial Z; -999 replicates 3.24)
0  #_placeholder for future growth feature
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
6 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
 6.8e-006 1.39e-005 2.85e-005 5.83e-005 0.000119295 0.000244164 0.000499664 0.00102224 0.00209011 0.00426839 0.00869559 0.0176272 0.0353807 0.0696495 0.132205 0.235559 0.381109 0.545888 0.692066 0.796225 0.859413 0.894075 0.912046 0.92109 0.925574 0.92778 0.928862 0.929391 0.92965 0.929776
2 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_	LO	HI	INIT		PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn		
0.001	0.70	0.088		-2.43	0.5323	3	-99	0	0	0	0	0	0	0	#	NatM_p_1_Fem_GP_1
2	15	7.63127		4	50	0	-99	0	0	0	0	0	0	0	#	L_at_Amin_Fem_GP_1
50	70	59.8661		60	50	0	-99	0	0	0	0	0	0	0	#	L_at_Amax_Fem_GP_1
0.02	0.21	0.132332	0.14	50	0	-99	0	0	0	0	0	0	0	#	VonBert_K_Fem_GP_1
0.02	0.21	0.106538	0.15	50	0	-99	0	0	0	0	0	0	0	#	CV_young_Fem_GP_1
0.01	0.20	0.029	        0.029	50	0	-99	0	0	0	0	0	0	0	#	CV_old_Fem_GP_1
0	1	1.18E-05	1.55E-05 50	6	-99	0	0	0	0	0	0	0	#	Wtlen_1_Fem
2	4	3.094		3.03	50	6	-99	0	0	0	0	0	0	0	#	Wtlen_2_Fem
40	41	40.5		40.5	50	6	-99	0	0	0	0	0	0	0	#	Mat50%_Fem
-3	3	-2.50E-01	-0.25	50	6	-99	0	0	0	0	0	0	0	#	Mat_slope_Fem
-3	3	0.2619		1	50	6	-99	0	0	0	0	0	0	0	#	Eggs/kg_inter_Fem
-1	1	0.0217		0	50	6	-99	0	0	0	0	0	0	0	#	Eggs/kg_slope_wt_Fem
0.001	0.70	0.088		-2.43	50	3	-99	0	0	0	0	0	0	0	#	NatM_p_1_Mal_GP_1
2	15	7.63127		4	50	0	-99	0	0	0	0	0	0	0	#	L_at_Amin_Mal_GP_1
50	70	53.4		60	50	0	-99	0	0	0	0	0	0	0	#	L_at_Amax_Mal_GP_1
0.02	0.21	0.165		0.14	50	0	-99	0	0	0	0	0	0	0	#	VonBert_K_Mal_GP_1
0.02	0.21	0.114		0.114	50	0	-99	0	0	0	0	0	0	0	#	CV_young_Mal_GP_1
0	0.2	0.042		0.042	50	0	-99	0	0	0	0	0	0	0	#	CV_old_Mal_GP_1
0	1	1.06E-05	1.55E-05 50	6	-99	0	0	0	0	0	0	0	#	Wtlen_1_Mal
2	4	3.123		3.03	50	6	-99	0	0	0	0	0	0	0	#	Wtlen_2_Mal
0	2	1		1	50	0	-99	0	0	0	0	0	0	0	#	RecrDist_GP_1
0	2	1		1	50	0	-99	0	0	0	0	0	0	0	#	RecrDist_Area_1
0	2	1		1	50	0	-99	0	0	0	0	0	0	0	#	RecrDist_Bseas_1
-1	1	1		1	50	0	-99	0	0	0	0	0	0	0	#	CohortGrowDev
0.01	0.99	0.5		0.5	0.5	0	-99	0	0	0	0	0	0	0	#	FracFemale_GP_1
#
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
#_Spawner-Recruitment
3 #_SR_function: 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepard_3Parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
             5            20         7.926         7.926             5             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
           0.2             1         0.773         0.773          0.15             2         -2          0          0          0          0          0          0          0 # SR_BH_steep
           0.5           1.2           0.7           0.7            99             0         -6          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0            99             0        -99          0          0          0          0          0          0          0 # SR_regime
             0             2             0             1            99             0        -99          0          0          0          0          0          0          0 # SR_autocorr
0 #do_recdev:  0=none; 1=devvector; 2=simple deviations
2014 # first year of main recr_devs; early devs can preceed this era
2014 # last year of main recr_devs; forecast devs start in following year
-1 #_recdev phase 
1 # (0/1) to read 13 advanced options
 0 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 -3 #_recdev_early_phase
 -5 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1955 #_last_yr_nobias_adj_in_MPD; begin of ramp
 1975 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 2012 #_last_yr_fullbias_adj_in_MPD
 2014 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
 1.0 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
 0 #_period of cycles in recruitment (N parms read below)
 -6 #min rec_dev
 6 #max rec_dev
 0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
# implementation error by year in forecast:  0 0 0 0 0 0 0 0 0 0 0 0
#
#Fishing Mortality info 
0.03 # F ballpark
-1999 # F ballpark year (neg value to disable)
1 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
0.9 # max F or harvest rate, depends on F_Method
# no additional F input needed for Fmethod 1
#
#_Q_setup for fleets with cpue or survey data
#_1:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm)
#_2:  extra input for link, i.e. mirror fleet
#_3:  0/1 to select extra sd parameter
#_4:  0/1 for biasadj or not
#_5:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
         2         1         0         0         0         1  #  Survey 1
         3         1         0         0         0         1  #  Survey 2
         4         1         0         0         0         1  #  Survey 3
         5         1         2         0         0         1  #  Depl Survey
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -15            15      -2.54651             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_Survey1
#             0           0.5      	0              0             1             0          2          0          0          0          0          0          0          0  #  ExtraSD_base_Survey1
           -15            15      -2.54651             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_Survey2
#             0           0.5      	0              0             1             0          2          0          0          0          0          0          0          0  #  ExtraSD_base_Survey2
           -15            15      -2.54651             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_Survey3
#             0           0.5      	0              0             1             0          2          0          0          0          0          0          0          0  #  ExtraSD_base_Survey3
           -15            15      -2.54651             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_Depl
#_no timevary Q parameters
#
#_size_selex_patterns
#_Pattern Discard Male Special
 1 0 0 0 # 1 Fishery
 15 0 0 1 # 2 Survey1
 15 0 0 1 # 3 Survey2
 15 0 0 1 # 4 Survey3
 10 0 0 0 # 5 Depl
#
#_age_selex_types
#_Pattern Discard Male Special
 10 0 0 0 # 1 Fishery
 10 0 0 0 # 2 Survey1
 10 0 0 0 # 3 Survey2
 10 0 0 0 # 4 Survey3
 10 0 0 0 # 5 Depl
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
            20            70          40.5          40.5            10             0          -1          0          0          0          0          0          0          0  #  SizeSel_P1
         0.001            50             7            15             5             0          -3          0          0          0          0          0          0          0  #  SizeSel_P2
# timevary selex parameters 
# info on dev vectors created for selex parms are reported with other devs after tag parameter section 
#
0   #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# Input variance adjustments factors: 
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
 -9999   1    0  # terminator
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 12 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark
#like_comp fleet  phase  value  sizefreq_method
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  0 # F_ballpark_lambda
0 # (0/1) read specs for more stddev reporting 
999

