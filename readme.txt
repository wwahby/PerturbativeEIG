This package contains matlab code for a standard EIG spectral partitioner, as well as a perturbative spectral partitioner.

## Check this file out first!
test_eig_partitoner.m is set up to run the standard solver once
for the original netlist chosen and once for the perturbed netlist, and
then it runs the perturbative solver once and generates a few plots to
compare cutsize and ratio cut.


## Warning about old versions of Matlab
Older versions of matlab are so bad at parsing text that it can actually
take way more time to parse the netlist than to actually calculate the eigenvectors!
I recommend using at least Matlab R2012a to avoid that situation, unless you really
need an excuse for a coffee break.


## Main routines
The main routines are eig_partitioner, for the standard solver, and
eig_partitioner_perturbed for the perturbative solver. The perturbative solver requires
a perturbed netlist. I've included premade perturbed netlists, with the naming convention
origfilename_addxyz.hgr where 0.xyz indicates the fraction of nets which were perturbed
(i.e. structP_add05.hgr has 5% of nets perturbed, structP_add005.hgr would have 0.5% perturbed, etc)


## Generating your own additive perturbations!
I've included a python script (Python 3.0 or better required) to generate random
additive perturbations to an arbitrary netlist

Usage: perturb_hgr <file to perturb> <output file name> <fraction of nets that
should get a random node added, i.e. 0.05 for 5% perturbation>
