% test eig_partitioner function
close all
clear all

%% Choose input file
% filename indicates original netlist
% filename_p indicates perturbed netlist

% filename = 'netlists/fract.hgr';
% filename_p = 'netlists/fract_add05.hgr';
% filename_p = 'netlists/fract_del05a.hgr';
% filename = 'netlists/p1.hgr';
% filename_p = 'netlists/p1_add05.hgr';
% filename_p = 'netlists/p1_del05a.hgr';
% filename = 'netlists/structP.hgr';
% filename_p = 'netlists/structP_add05.hgr';
% filename_p = 'netlists/structP_del05a.hgr';
% filename = 'netlists/p2.hgr';
% filename_p = 'netlists/p2_add05.hgr';
% filename_p = 'netlists/p2_del05a.hgr';
% filename = 'netlists/biomedP.hgr';
% filename_p = 'netlists/biomedP_add05.hgr';
% filename_p = 'netlists/biomedP_del05a.hgr';
% filename_p = 'netlists/biomedP_add005.hgr';
% filename = 'netlists/industry2.hgr';
% filename_p = 'netlists/industry2_add05.hgr';
% filename_p = 'netlists/industry2_del05a.hgr';
% filename = 'netlists/ibm01.hgr';
% filename_p = 'netlists/ibm01_add05.hgr';
% filename_p = 'netlists/ibm01_del05a.hgr';
% filename = 'netlists/industry3.hgr';
% filename_p = 'netlists/industry3_add05.hgr';
% filename_p = 'netlists/industry3_del05a.hgr';
%% Warning! Here be monsters. Have more than 3GB RAM free before you try to run these
% filename = 'netlists/ibm10.hgr';
% filename_p = 'netlists/ibm10_add01.hgr';
% filename_p = 'netlists/ibm10_add05.hgr';
% filename_p = 'netlists/ibm10_del05a.hgr';
% filename = 'netlists/ibm18.hgr';
% filename_p = 'netlists/ibm18_add05.hgr';
% filename_p = 'netlists/ibm18_del05a.hgr';


%% Parameters
num_eigs = 7; % number of eigenvalues and eigenvectors to use for the perturbed solution
node_areas = 1; % all nodes have unit area
area_constraint = 1e-5; % schmoo across all possible area splits (set this intentionally low so we have a minimum of 1e-3% of the nodes on one side and the rest on the other, and try all partitions between that and 1-1e5)



%% Run EIG to find exact partitioning results for original and perturbed systems
[metrics times matrices eigs] = eig_partitioner(filename,num_eigs,node_areas,area_constraint); % normal
[metrics_pe times_pe matrices_pe eigs_pe] = eig_partitioner(filename_p,2,node_areas,area_constraint); % exact perturbation, only need 2 eigenvalues here



%% Run the perturbed solver
% If you want to run the perturbed solver WITHOUT running the full
% partitioner for the perturbed solution, you'll need to read in the
% netlist with the line below
% [Qpe Dpe Ape] = parse_hgr_sparse_alt3(filename_p);

% Construct perturbation matrices (i.e. the difference between the new
% netlist and the old netlist)
Qp = matrices_pe.laplacian - matrices.laplacian;
Ape = matrices_pe.adjacency;
[metrics_p times_p matrices_p eigs_p] = eig_partitioner_perturbed(filename_p,matrices.laplacian,eigs.vals,eigs.vecs,node_areas,area_constraint); % approx perturbation




%% Full partitioning plots (here you can see the difference between a 40/60 split and a 60/50 split)
% figure(6)
% clf
% plot(metrics.partition_ratio, metrics.cutsize.vec,'k')
% hold on
% plot(metrics_pe.partition_ratio, metrics_pe.cutsize.vec,'b')
% plot(metrics_p.partition_ratio, metrics_p.cutsize.vec,'r')
% xlabel('partition ratio')
% ylabel('cutsize')
% title('Cutsize')
% fixfigs(6,3,14,12)
% 
% figure(7)
% clf
% plot(metrics.partition_ratio, metrics.ratio_cut.vec,'k')
% hold on
% plot(metrics_pe.partition_ratio, metrics_pe.ratio_cut.vec,'b')
% plot(metrics_p.partition_ratio, metrics_p.ratio_cut.vec,'r')
% xlabel('partition ratio')
% ylabel('ratio cut')
% title('Ratio cut')
% fixfigs(7,3,14,12)

%% Skew plots (just show best result from 0% skew up to full skew, i.e. compare a 40/60 split and a 60/40 split and show only the best result for skew of 0.1)
figure(8)
clf
plot(metrics.skew,metrics.cutsize.skew.vec,'k')
hold on
plot(metrics_pe.skew,metrics_pe.cutsize.skew.vec,'b')
plot(metrics_p.skew,metrics_p.cutsize.skew.vec,'r')
xlabel('skew')
ylabel('cutsize');
title('Cutsize')
xlim([0 0.5])
fixfigs(8,3,14,12)

figure(9)
clf
plot(metrics.skew,metrics.ratio_cut.skew.vec,'k')
hold on
plot(metrics_pe.skew,metrics_pe.ratio_cut.skew.vec,'b')
plot(metrics_p.skew,metrics_p.ratio_cut.skew.vec,'r')
xlabel('skew')
ylabel('ratio cut');
title('Ratio cut')
xlim([0 0.5])
fixfigs(9,3,14,12)