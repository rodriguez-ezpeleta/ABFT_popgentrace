import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

'''
Models for testing various three population scenarios.

############################################

Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated July 2018
'''

##########################################################################################
#Basic models of (no gene flow / gene flow) between (all / some) population pairs
##########################################################################################

##########################################################################################
#Various models based on forest refugia timing, all with symmetric gene flow estimates
##########################################################################################

def refugia_adj_1(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur. Split between pops
    2 and 3, gene flow does not occur. Period of symmetric secondary contact occurs between 
    adjacent populations (ie 1<->2, 2<->3, but not 1<->3) after all splits are complete. 
    'longest isolation'
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2: Migration rate between populations 2 and 3 
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    T3: The scaled time between the split of pops 2 and 3 and the present (in units of 2*Na generations).	
    """
    #9 parameters
    nu1, nuA, nu2, nu3, m1, m2, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    
    phi = Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


def refugia_adj_2(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur. Split between pops
    2 and 3, with gene flow. After appearance of 2 and 3, gene flow also occurs between 1 
    and 2.
    'shorter isolation'
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2: Migration rate between populations 2 and 3 
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #8 parameters
    nu1, nuA, nu2, nu3, m1, m2, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


def refugia_adj_3(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur, but then 
    secondary contact occurs. Split between pops 2 and 3 occurs with gene flow, and gene flow
    happens between 1 and 2 as well.
    'shortest isolation'
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    mA: Migration rate between population 1 and population (2,3)
    m1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2: Migration rate between populations 2 and 3 
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #10 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, T1a, T1b, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1a, nu1=nu1, nu2=nuA, m12=0, m21=0)
    
    phi = Integration.two_pops(phi, xx, T1b, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


##########################################################################################
#Various models based on ancient migration and contemporary isolation
##########################################################################################

def ancmig_adj_3(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), with gene flow, which then stops. Split 
    between pops 2 and 3, gene flow does not occur at all.
    'longest isolation'
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    mA: Migration rate between population 1 and population (2,3)
    m1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2: Migration rate between populations 2 and 3 
    m3: Migration rate between populations 1 and 3
    T1a: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T1b: The scaled time for no gene flow between pops 1 and (2, 3) (in units of 2*Na generations).
    T2: The scaled time between the cessation of gene flow and the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #8 parameters
    nu1, nuA, nu2, nu3, mA, T1a, T1b, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1a, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    
    phi = Integration.two_pops(phi, xx, T1b, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def ancmig_adj_2(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), with gene flow. Split 
    between pops 2 and 3, and all gene flow ceases.
    'shorter isolation'
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    mA: Migration rate between population 1 and population (2,3)
    m1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2: Migration rate between populations 2 and 3 
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #7 parameters
    nu1, nuA, nu2, nu3, mA, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def ancmig_adj_1(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), with gene flow. Split 
    between pops 2 and 3 with gene flow, then all gene flow ceases.
    'shortest isolation'
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    mA: Migration rate between population 1 and population (2,3)
    m1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2: Migration rate between populations 2 and 3 
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 and present (in units of 2*Na generations).
    """
    #10 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2, T3 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)
    
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)

    phi = Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


##########################################################################################
#Simultaneous split models (with/without migration/secondary contact and size changes)
##########################################################################################
'''
def sim_split_no_mig(params, ns, pts):
    """
    Model with simultaneous split between pop 1, 2 and 3, gene flow does not occur. 
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
	T1: The scaled time between the split and the present (in units of 2*Na generations).
	"""
    #4 parameters
    nu1, nu2, nu3, T1 = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T1, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs
'''	
def sim_split_no_mig_size(params, ns, pts):
    """
    Model with simultaneous split between pop 1, 2 and 3, gene flow does not occur, but size change does. 
    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu3a: Size of population 3 after split.
    nu1b: Size of population 1 after size change.
    nu2b: Size of population 2 after size change.
    nu3b: Size of population 3 after size change.
    T1: The scaled time between the split and the size change (in units of 2*Na generations).
    T2: The scaled time between the size change and present.
	"""
    #8 parameters
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, T1, T2 = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T1, nu1=nu1a, nu2=nu2a, nu3=nu3a, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1b, nu2=nu2b, nu3=nu3b, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs
	
def sim_split_sym_mig_all(params, ns, pts):
    """
    Model with simultaneous split between pop 1, 2 and 3
    Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2: Migration rate between populations 2 and 3 
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split and the present (in units of 2*Na generations).
	"""
	#7 parameters
    nu1, nu2, nu3, m1, m2, m3, T1 = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T1, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m3, m31=m3)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


def sim_split_sym_mig_adjacent(params, ns, pts):
    """
    Model with simultaneous split between pop 1, 2 and 3
    Migration is symmetrical between 'adjacent' population pairs (ie 1<->2, 2<->3, but not 1<->3).
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2: Migration rate between populations 2 and 3 
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split and the present (in units of 2*Na generations).
    """
    #6 parameters
    nu1, nu2, nu3, m1, m2, T1 = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T1, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs
  
def sim_split_refugia_sym_mig_all(params, ns, pts):
    """
    Model with simultaneous split between pop 1, 2 and 3 followed by isolation. Period of 
    symmetric secondary contact occurs between all populations after all splits are complete. 
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2: Migration rate between populations 2 and 3 
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split and the present (in units of 2*Na generations).
    T2: The scaled time between the migration (secondary contact) and the present (in units of 2*Na generations).
    """
    #8 parameters
    nu1, nu2, nu3, m1, m2, m3, T1, T2 = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T1, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m3, m31=m3)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


def sim_split_refugia_sym_mig_adjacent(params, ns, pts):
    """
    Model with simultaneous split between pop 1, 2 and 3 followed by isolation. Period of 
    symmetric secondary contact occurs between adjacent populations (ie 1<->2, 2<->3, but 
    not 1<->3) after all splits are complete. 
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2: Migration rate between populations 2 and 3 
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split and the present (in units of 2*Na generations).
    T2: The scaled time between the migration (secondary contact) and the present (in units of 2*Na generations).
	"""
    #7 parameters
    nu1, nu2, nu3, m1, m2, T1, T2 = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T1, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


###############################################################
#### Models with extra size change step (potential human impoact)
###############################################################

def split_nomig_size(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration does not occur between any population pair, size change
    nu1a: Size of population 1 after split.
    nuA: Size of ancestral population (pops 2 and 3)
	nu2a: Size of population 2 after split.
    nu3a: Size of population 3 after split.
    nu1b: Size of population 1 after size change.
	nu2b: Size of population 2 after size change.
    nu3b: Size of population 3 after size change.
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    T3: The scaled time between the size change and the present (in units of 2*Na generations).
    """
    #10 parameters	
    nu1a, nuA, nu2a, nu3a, nu1b, nu2b, nu3b, T1, T2, T3 = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1a, nu2a, m12=0, m21=0)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T2, nu1a, nu2a, nu3a, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T3, nu1b, nu2b, nu3b, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def ancmig_2_size(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), with gene flow. Split 
    between pops 2 and 3, and all gene flow ceases, then size change
    shorter isolation
    nu1a: Size of population 1 after split.
    nuA: Size of ancestral population (pops 2 and 3)
	nu2a: Size of population 2 after split.
    nu3a: Size of population 3 after split.
    nu1b: Size of population 1 after size change.
	nu2b: Size of population 2 after size change.
    nu3b: Size of population 3 after size change.
    mA: Migration rate between pops 1 and (2 and 3)
	T1: The scaled time between the split of pops 1 vs 2 and 3 (ancient migration)(in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    T3: The scaled time between the size change and the present (in units of 2*Na generations).
    """
    #11 parameters
    nu1a, nuA, nu2a, nu3a, nu1b, nu2b, nu3b, mA, T1, T2, T3 = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1a, nuA, m12=mA, m21=mA)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T2, nu1a, nu2a, nu3a, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    phi = Integration.three_pops(phi, xx, T3, nu1b, nu2b, nu3b, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


###############################################
###### HERE IS WHERE MY MODELS START ##########
###############################################

###############################################################
#### Models with one split and the third population is generated as a hybrid population
###############################################################


def split_admix_nomig(params, ns, pts):
    """
    Model with split between pop 1 and 3, then they admix to generate population 2.
    Migration does not occur between any population pair, except for the admixture event.
    nu1: Size of population 1.
    nu2: Size of population 2.
    nu3: Size of population 3 after admixture.
    T1: The scaled time between the split of pops 1 and 2 and admixture (in units of 2*Na generations).
    T2: The scaled time since the admixture event (in units of 2*Na generations).
    f: fraction of the gene influx coming from population 1. Influx coming from population 2 would be 1-f.
    O: proportion of correctly orientated sites
    """
    #7 parameters
    nu1,nu2,nu3,T1,T2,f,O = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nu2, m12=0, m21=0)
    phi = PhiManip.phi_2D_to_3D_admix(phi,f, xx, xx, xx)
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fsO = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    fsM = Numerics.reverse_array (fsO) 
    fs = O*fsO+(1-O)*fsM

    return fs


def split_admix_mig(params, ns, pts):
    """
    Model with split between pop 1 and 3, then they admix to generate population 2.
    Migration does not occur between any population pair, except for the admixture event.
    nu1: Size of population 1.
    nu2: Size of population 2.
    nu3: Size of population 3 after admixture.
    T1: The scaled time between the split of pops 1 and 2 and admixture (in units of 2*Na generations).
    T2: The scaled time since the admixture event (in units of 2*Na generations).
    m12: rate of migration from pop2 into pop1. It is equal to the fraction of individuals each generation in pop1 that are new migrants from pop 2, times the 2Nref.
    m21: rate of migration from pop1 into pop2.
    m23: rate of migration from pop3 into pop2.
    m32: rate of migration from pop2 into pop3.
    m13: rate of migration from pop3 into pop1.
    m31: rate of migration from pop1 into pop3.
    f: fraction of the gene influx coming from population 1. Influx coming from population 2 would be 1-f.
    O: proportion of correctly orientated sites
    """
    #13 parameters
    nu1,nu2,nu3,T1,T2,m12,m21,m23,m32,m13,m31,f,O = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nu2, m12=0, m21=0)
    phi = PhiManip.phi_2D_to_3D_admix(phi,f, xx, xx, xx)
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m12, m21=m21, m23=m23, m32=m32, m13=m13, m31=m31)
    fsO = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    fsM = Numerics.reverse_array (fsO)
    fs = O*fsO+(1-O)*fsM

    return fs


##########################################################################################
#Simultaneous split models
##########################################################################################

def sim_split_no_mig(params, ns, pts):
    """
    Model with simultaneous split between pop 1, 2 and 3, gene flow does not occur.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    T1: The scaled time between the split and the present (in units of 2*Na generations).
    O: proportion of correctly orientated sites. 
    """
    #5 parameters
    nu1, nu2, nu3, T1, O = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T1, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)
    fsO = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    fsM = Numerics.reverse_array (fsO)
    fs = O*fsO+(1-O)*fsM

    return fs



def sim_split_mig(params, ns, pts):
    """
    Model with simultaneous split between pop 1, 2 and 3, gene flow does occur.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    T1: The scaled time between the split and the present (in units of 2*Na generations).
    m12: rate of migration from pop2 into pop1. It is equal to the fraction of individuals each generation in pop1 that are new migrants from pop 2, times the 2Nref.
    m21: rate of migration from pop1 into pop2.
    m23: rate of migration from pop3 into pop2.
    m32: rate of migration from pop2 into pop3.
    m13: rate of migration from pop3 into pop1.
    m31: rate of migration from pop1 into pop3.
    O: proportion of correctly orientated sites.
    """
    #11 parameters
    nu1,nu2,nu3,T1,m12,m21,m23,m32,m13,m31,O = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T1, nu1=nu1, nu2=nu2, nu3=nu3, m12=m12, m21=m21, m23=m23, m32=m32, m13=m13, m31=m31)
    fsO = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    fsM = Numerics.reverse_array (fsO)
    fs = O*fsO+(1-O)*fsM

    return fs

##########################################################################################
# 2 Split Models
##########################################################################################

def split_nomig(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration does not occur between any population pair.
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    O: proportion of correctly orientated sites in the input_file
    """
    #7 parameters
    nu1, nuA, nu2, nu3, T1, T2, O = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=0, m21=0)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)

    fsO = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    fsM = Numerics.reverse_array (fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def split_mig(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration does not occur between any population pair.
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    m12: rate of migration from pop2 into pop1. It is equal to the fraction of individuals each generation in pop1 that are new migrants from pop 2, times the 2Nref.
    m21: rate of migration from pop1 into pop2.
    m23: rate of migration from pop3 into pop2.
    m32: rate of migration from pop2 into pop3.
    m13: rate of migration from pop3 into pop1.
    m31: rate of migration from pop1 into pop3.
    O: proportion of correctly orientated sites in the input_file
    """
    #13 parameters
    nu1, nuA, nu2, nu3, T1, T2, m12, m21, m23, m32, m13, m31, O = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=0, m21=0)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m12, m21=m21, m23=m23, m32=m32, m13=m13, m31=m31)

    fsO = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    fsM = Numerics.reverse_array (fsO)
    fs = O*fsO+(1-O)*fsM
    return fs



def split_2_nomig(params, ns, pts):
    """
    Model with split between pop (1,3) and 2, then split between 1 and 3.
    Migration does not occur between any population pair.
    nuA: Size of population (1,3) after split.
    nu2: Size of population 2 after split from (1,3).
    nu1: Size of population 1 after split.
    nu3: Size of population 3 after split.
    T1: The scaled time between the split of pops 1 and 3 vs 2 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 1 and 3 (in units of 2*Na generations).
    O: proportion of correctly orientated sites in the input_file
    """
    #7 parameters
    nuA, nu2, nu1, nu3, T1, T2, O = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=nuA, nu2=nu2, m12=0, m21=0)
    phi = PhiManip.phi_2D_to_3D_split_1(xx, phi)
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)

    fsO = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    fsM = Numerics.reverse_array (fsO)
    fs = O*fsO+(1-O)*fsM
    return fs

def split_2_mig(params, ns, pts):
    """
    Model with split between pop (1,3 and 2, then split between 1 and 3.
    Migration does occur between any population pair.
    nuA: Size of population (1,3) after split.    
    nu2: Size of population 2 after split from (1,3).
    nu1: Size of population 1 after split.
    nu3: Size of population 3 after split.
    T1: The scaled time between the split of pops 1 and 3 vs 2 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 1 and 3 (in units of 2*Na generations).
    m12: rate of migration from pop2 into pop1. It is equal to the fraction of individuals each generation in pop1 that are new migrants from pop 2, times the 2Nref.
    m21: rate of migration from pop1 into pop2.
    m23: rate of migration from pop3 into pop2.
    m32: rate of migration from pop2 into pop3.
    m13: rate of migration from pop3 into pop1.
    m31: rate of migration from pop1 into pop3.
    O: proportion of correctly orientated sites in the input_file
    """
    #13 parameters
    nuA, nu2, nu1, nu3, T1, T2, m12, m21, m23, m32, m13, m31, O = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=nuA, nu2=nu2, m12=0, m21=0)
    phi = PhiManip.phi_2D_to_3D_split_1(xx, phi)
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m12, m21=m21, m23=m23, m32=m32, m13=m13, m31=m31)

    fsO = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    fsM = Numerics.reverse_array (fsO)
    fs = O*fsO+(1-O)*fsM
    return fs


