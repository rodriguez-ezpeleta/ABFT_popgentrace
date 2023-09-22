import sys
import os
import numpy
import dadi
import pylab
from datetime import datetime
import Optimize_Functions
import Models_3D

'''
Usage: python dadi_Run_3D_Set.py

This is a modified version of the 'dadi_Run_Optimizations.py' script in which
we run optimizations for 3D comparisons for a large set of models that have been
made available as part of published works. These models are stored in the
Models_3D.py script, and will be called directly here. The user can delete or
comment out models to analyze a subset of the models available. 

This script must be in the same working directory as Optimize_Functions.py, which
contains all the functions necessary, as well as the  Models_3D.py script, which
has all the model definitions.


General workflow:
 The optimization routine runs a user-defined number of rounds, each with a user-defined
 or predefined number of replicates. The starting parameters are initially random, but after
 each round is complete the parameters of the best scoring replicate from that round are
 used to generate perturbed starting parameters for the replicates of the subsequent round.
 The arguments controlling steps of the optimization algorithm (maxiter) and perturbation
 of starting parameters (fold) can be supplied by the user for more control across rounds.
 The user can also supply their own set of initial parameters, or set custom bounds on the
 parameters (upper_bound and lower_bound) to meet specific model needs. This flexibility
 should allow these scripts to be generally useful for model-fitting with any data set.

 
Outputs:
 For each model run, there will be a log file showing the optimization steps per replicate
 and a summary file that has all the important information. Here is an example of the output
 from a summary file, which will be in tab-delimited format:
 
 Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T)
 sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
 sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
 sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
 sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
 sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688


Notes/Caveats:
 The likelihood and AIC returned represent the true likelihood only if the SNPs are
 unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this
 is true, but if SNPs are linked across loci then the likelihood is actually a composite
 likelihood and using something like AIC is no longer appropriate for model comparisons.
 See the discussion group for more information on this subject. 

Citations:
 If you use these scripts or the sets of diversification models for your work, please
 cite the following publications:
    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 5245-5263.
    doi: 10.1111/mec.14266
    
    Barratt, C.D., Bwong, B.A., Jehle, R., Liedtke, H.C., Nagel, P., Onstein, R.E., 
	Portik, D.M., Streicher, J.W., and S.P. Loader. Vanishing refuge: testing the 
	forest refuge hypothesis in coastal East Africa using genome-wide sequence data 
	for five co-distributed amphibians. In Review, Molecular Ecology.

 If you are interesting in contributing your 3D models to this workflow, please email me!

-------------------------
Written for Python 2.7
Python modules required:
-Numpy
-Scipy
-dadi
-------------------------

Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated July 2018
'''

#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================

#**************
snps = "/home/natalia/bin/dadi_pipeline/Three_Population_Pipeline/ABFT_test/input_file.txt"

#Create python dictionary from snps file
dd = dadi.Misc.make_data_dict(snps)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=['MED','GOM','SLOPE']

#**************
#projection sizes, in ALLELES not individuals
proj = [20,20,20]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = True)
fs.mask [0,1,0] = True 
fs.mask [1,0,0] = True 
fs.mask [0,0,1] = True

#print some useful information about the afs or jsfs
print "\n\n============================================================================\nData for site frequency spectrum\n============================================================================\n"
print "projection", proj
print "sample sizes", fs.sample_sizes
sfs_sum = numpy.around(fs.S(), 2)
print "Sum of SFS = ", sfs_sum, '\n', '\n'

#================================================================================
# Calling external 3D models from the Models_3D.py script
#================================================================================
'''
 We will use a function from the Optimize_Functions.py script for our optimization routines:
 
 Optimize_Routine(fs, pts, outfile, model_name, func, rounds, param_number, fs_folded=True, reps=None, maxiters=None, folds=None, in_params=None, in_upper=None, in_lower=None, param_labels=" ")
 
   Mandatory Arguments =
    fs:  spectrum object name
    pts: grid size for extrapolation, list of three values
    outfile:  prefix for output naming
    model_name: a label to help label the output files; ex. "no_mig"
    func: access the model function from within 'moments_Run_Optimizations.py' or from a separate python model script, ex. after importing Models_3D, calling Models_3D.split_nomig
    rounds: number of optimization rounds to perform
    param_number: number of parameters in the model selected (can count in params line for the model)q
    fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False).

   Optional Arguments =
     reps: a list of integers controlling the number of replicates in each of the optimization rounds
     maxiters: a list of integers controlling the maxiter argument in each of the optimization rounds
     folds: a list of integers controlling the fold argument when perturbing input parameter values
     in_params: a list of parameter values 
     in_upper: a list of upper bound values
     in_lower: a list of lower bound values
     param_labels: list of labels for parameters that will be written to the output file to keep track of their order

Below, I give all the necessary information to call each model available in the
Models_2D.py script. I have set the optimization routine to be the same for each
model using the optional lists below, which are included as optional arguments for
each model. This particular configuration will run 4 rounds as follows:
Round1 - 10 replicates, maxiter = 3, fold = 3
Round2 - 20 replicates, maxiter = 5, fold = 2
Round3 - 30 replicates, maxiter = 10, fold = 2
Round4 - 40 replicates, maxiter = 15, fold = 1

If this script was run as is, each model would be called and optimized sequentially;
this could take a very long time for 3D models. For your actual analyses, I strongly 
recommend creating multiple scripts with only a few models each and running them
independently. 

'''


#create a prefix based on the population names to label the output files
#ex. Pop1_Pop2_Pop3
prefix = "_".join(pop_ids)

#**************
#make sure to define your extrapolation grid size (based on your projections)
pts = [20,30,40]

#**************
#Set the number of rounds here
rounds = 5
 

#define the lists for optional arguments
#you can change these to alter the settings of the optimization routine
reps = [10,20,30,40,50]
maxiters = [3,5,10,15,20]
folds = [3,2,2,1,1]


#**************
#Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = False


'''
Diversification Model Set

This first set of models come from the following publication:

    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 5245-5263.
    doi: 10.1111/mec.14266

'''
'''
# Split into three populations, no migration.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_nomig", Models_3D.split_nomig, rounds, 6, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, T1, T2")

# Split into three populations, symmetric migration between all populations (1<->2, 2<->3, and 1<->3).
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_symmig_all", Models_3D.split_symmig_all, rounds, 10, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2")

# Split into three populations, symmetric migration between 'adjacent' populations (1<->2, 2<->3, but not 1<->3).
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_symmig_adjacent", Models_3D.split_symmig_adjacent, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2")

# Adjacent Secondary contact, longest isolation - Split between pop 1 and (2,3) with no migration, then split between pop 2 and 3 with no migration. Period of symmetric secondary contact occurs between adjacent populations (ie 1<->2, 2<->3, but not 1<->3) after all splits are complete.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "refugia_adj_1", Models_3D.refugia_adj_1, rounds, 9, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, m1, m2, T1, T2, T3")

# Adjacent Secondary contact, shorter isolation - Split between pop 1 and (2,3), gene flow does not occur. Split between pop 2 and 3 occurs with gene flow. After appearance of 2 and 3, gene flow also occurs between 1 and 2.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "refugia_adj_2", Models_3D.refugia_adj_2, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, m1, m2, T1, T2")

# Adjacent Secondary contact, shortest isolation - Split between pop 1 and (2,3) with no migration. Split between pop 2 and 3 occurs with gene flow, and gene flow occurs between 1 and 2 as well.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "refugia_adj_3", Models_3D.refugia_adj_3, rounds, 10, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, mA, m1, m2, T1a, T1b, T2")

# Adjacent Ancient migration, longest isolation - Split between pop 1 and (2,3) with gene flow, which then stops. Split between pop 2 and 3, migration does not occur at all.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "ancmig_adj_3", Models_3D.ancmig_adj_3, rounds, 8, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, mA, T1a, T1b, T2")

# Adjacent Ancient migration, shorter isolation - Split between pop 1 and (2,3) with gene flow. Split between pop 2 and 3 with no migration between any populations.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "ancmig_adj_2", Models_3D.ancmig_adj_2, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, mA, T1, T2")

# Adjacent Ancient migration, shortest isolation - Split between pop 1 and (2,3) with gene flow. Split between pop 2 and 3 with gene flow, then all gene flow ceases.
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "ancmig_adj_1", Models_3D.ancmig_adj_1, rounds, 10, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2, T3")
'''

'''
This second set of models were developed for the following publication:

	Barratt, C.D., Bwong, B.A., Jehle, R., Liedtke, H.C., Nagel, P., Onstein, R.E., 
	Portik, D.M., Streicher, J.W., and S.P. Loader. Vanishing refuge: testing the 
	forest refuge hypothesis in coastal East Africa using genome-wide sequence data 
	for five co-distributed amphibians. In Review, Molecular Ecology.

'''

##########################################
##### HERE IS WHERE MY MODELS START ###### 
##########################################

'''
##### SCRIPT FOR GRID #######
argv = sys.argv[1:]
model = argv[0]
iter = arg[1]

if model==split_admix_nomig:
	prefix=model + "_" + iter
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_admix_nomig", Models_3D.split_admix_nomig , rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,0.5,0.5,0.95], in_upper=[20,20,20,10,5,0.99,0.995], in_lower=[0.01,0.01,0.01,0.01,0.01,0.01,0.9], param_labels = "nu1, nu2, nu3, T1, T2, f, O")
elif model==split_admix_mig:
        prefix=model + "_" + iter
 	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_admix_mig", Models_3D.split_admix_mig , rounds, 13, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,0.5,1,1,1,1,1,1,0.5,0.95], in_upper=[20,20,20,10,5,20,20,20,20,20,20,0.99,0.995], in_lower=[0.01,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0.01,0.9], param_labels = "nu1, nu2, nu3, T1, T2, m12, m21, m23, m32, m13, m31, f, O")
elif model==sim_split_nomig:
	prefix=model + "_" + iter
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sim_split_no_mig", Models_3D.sim_split_no_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,0.95], in_upper=[20,20,20,10,0.995], in_lower=[0.01,0.01,0.01,0.01,0.9], param_labels = "nu1, nu2, nu3, T1, O")
elif model==sim_split_mig:
	prefix=model + "_" + iter
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sim_split_mig", Models_3D.sim_split_mig, rounds, 11, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,1,1,1,1,1,1,0.95], in_upper=[20,20,20,10,20,20,20,20,20,20,0.995], in_lower=[0.01,0.01,0.01,0.01,0,0,0,0,0,0,0.9], param_labels = "nu1, nu2, nu3, T1, m12, m21, m23, m32, m13, m31, O")
elif model==split_nomig:
	prefix=model + "_" + iter
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_nomig", Models_3D.split_nomig, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,1,0.5,0.95], in_upper=[20,20,20,20,10,5,0.995], in_lower=[0.01,0.01,0.01,0.01,0.01,0.01,0.9], param_labels = "nu1, nuA, nu2, nu3, T1, T2, O")
elif model==split_mig:
	prefix=model + "_" + iter
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_mig", Models_3D.split_mig, rounds, 13, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,1,0.5,1,1,1,1,1,1,0.95], in_upper=[20,20,20,20,10,5,20,20,20,20,20,20,0.995], in_lower=[0.01,0.01,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0.9], param_labels = "nu1, nuA, nu2, nu3, T1, T2, m12, m21, m23, m32, m13, m31, O")
elif model==split_2_nomig:
	prefix=model + "_" + iter
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_2_nomig", Models_3D.split_2_nomig, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,1,0.5,0.95], in_upper=[20,20,20,20,10,5,0.995], in_lower=[0.01,0.01,0.01,0.01,0.01,0.01,0.9], param_labels = "nuA, nu2, nu1, nu3, T1, T2, O")
elif model==split_2_mig:
	prefix=model + "_" + iter
	Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_2_mig", Models_3D.split_2_mig, rounds, 13, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,1,0.5,1,1,1,1,1,1,0.95], in_upper=[20,20,20,20,10,5,20,20,20,20,20,20,0.995], in_lower=[0.01,0.01,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0.9], param_labels = "nuA, nu2, nu1, nu3, T1, T2, m12, m21, m23, m32, m13, m31, O")



#### Models with one split and the third population is generated as a hybrid population

# Split into two populations with posterior admixture event and generation of a thirds population

for i in range(1,11):
    prefix = "V10_Number_{}".format(i)
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_admix_nomig", Models_3D.split_admix_nomig , rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,0.5,0.5,0.95], in_upper=[20,20,20,10,5,0.99,0.995], in_lower=[0.01,0.01,0.01,0.01,0.01,0.01,0.9], param_labels = "nu1, nu2, nu3, T1, T2, f, O")

for i in range(1,11):
    prefix = "V10_Number_{}".format(i)
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_admix_mig", Models_3D.split_admix_mig , rounds, 13, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,0.5,1,1,1,1,1,1,0.5,0.95], in_upper=[20,20,20,10,5,20,20,20,20,20,20,0.99,0.995], in_lower=[0.01,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0.01,0.9], param_labels = "nu1, nu2, nu3, T1, T2, m12, m21, m23, m32, m13, m31, f, O") 



############# Models with simultaneous population splitting

# Simultaneous split into three populations

for i in range(1,11):
    prefix = "V10_Number_{}".format(i)
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sim_split_no_mig", Models_3D.sim_split_no_mig, rounds, 5, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,0.95], in_upper=[20,20,20,10,0.995], in_lower=[0.01,0.01,0.01,0.01,0.9], param_labels = "nu1, nu2, nu3, T1, O")

for i in range(1,11):
    prefix = "V10_Number_{}".format(i)
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "sim_split_mig", Models_3D.sim_split_mig, rounds, 11, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,1,1,1,1,1,1,0.95], in_upper=[20,20,20,10,20,20,20,20,20,20,0.995], in_lower=[0.01,0.01,0.01,0.01,0,0,0,0,0,0,0.9], param_labels = "nu1, nu2, nu3, T1, m12, m21, m23, m32, m13, m31, O")

############# Models with two splits

# Split 1Vs2, after split 2Vs3, no migration.
for i in range(1,11):
    prefix = "V10_Number_{}".format(i)
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_nomig", Models_3D.split_nomig, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,1,0.5,0.95], in_upper=[20,20,20,20,10,5,0.995], in_lower=[0.01,0.01,0.01,0.01,0.01,0.01,0.9], param_labels = "nu1, nuA, nu2, nu3, T1, T2, O")

# Split 1Vs2, after split 2Vs3, WITH migration.
for i in range(1,11):
    prefix = "V10_Number_{}".format(i)
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_mig", Models_3D.split_mig, rounds, 13, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,1,0.5,1,1,1,1,1,1,0.95], in_upper=[20,20,20,20,10,5,20,20,20,20,20,20,0.995], in_lower=[0.01,0.01,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0.9], param_labels = "nu1, nuA, nu2, nu3, T1, T2, m12, m21, m23, m32, m13, m31, O")


# Split 1Vs2, after split 1Vs3, no migration.
for i in range(1,11):
    prefix = "V10_Number_{}".format(i)
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_2_nomig", Models_3D.split_2_nomig, rounds, 7, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,1,0.5,0.95], in_upper=[20,20,20,20,10,5,0.995], in_lower=[0.01,0.01,0.01,0.01,0.01,0.01,0.9], param_labels = "nuA, nu2, nu1, nu3, T1, T2, O")

# Split 1Vs2, after split 1Vs3, WITH migration.
for i in range(1,11):
    prefix = "V10_Number_{}".format(i)
    Optimize_Functions.Optimize_Routine(fs, pts, prefix, "split_2_mig", Models_3D.split_2_mig, rounds, 13, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_params=[1,1,1,1,1,0.5,1,1,1,1,1,1,0.95], in_upper=[20,20,20,20,10,5,20,20,20,20,20,20,0.995], in_lower=[0.01,0.01,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0.9], param_labels = "nuA, nu2, nu1, nu3, T1, T2, m12, m21, m23, m32, m13, m31, O")
'''
