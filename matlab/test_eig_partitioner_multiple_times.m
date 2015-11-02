% test eig_partitioner function
close all
clear all

%% Choose input file
% filename indicates original netlist
% filename_p indicates perturbed netlist

% filename = '../netlists/fract.hgr';
% filename_p = '../netlists/fract_add05.hgr';
% filename_p = '../netlists/fract_del05a.hgr';
% filename = '../netlists/p1.hgr';
% filename_p = '../netlists/p1_add05.hgr';
% filename_p = '../netlists/p1_del05a.hgr';
% filename = '../netlists/structP.hgr';
% filename_p = '../netlists/structP_add05.hgr';
% filename_p = '../netlists/structP_del05a.hgr';
filename = '../netlists/p2.hgr';
filename_p = '../netlists/p2_add05.hgr';
% filename_p = '../netlists/p2_del05a.hgr';
% filename = '../netlists/biomedP.hgr';
% filename_p = '../netlists/biomedP_add05.hgr';
% filename_p = '../netlists/biomedP_del05a.hgr';
% filename_p = '../netlists/biomedP_add005.hgr';
% filename = '../netlists/industry2.hgr';
% filename_p = '../netlists/industry2_add05.hgr';
% filename_p = '../netlists/industry2_del05a.hgr';
% filename = '../netlists/ibm01.hgr';
% filename_p = '../netlists/ibm01_add05.hgr';
% filename_p = '../netlists/ibm01_del05a.hgr';
% filename = '../netlists/industry3.hgr';
% filename_p = '../netlists/industry3_add05.hgr';
% filename_p = '../netlists/industry3_del05a.hgr';
% filename = '../netlists/synnet_no1_2.hgr';
% filename_p = '../netlists/synnet_no1_mod_2.hgr';
%% Warning! Here be monsters. Have more than 3GB RAM free before you try to run these
% filename = '../netlists/ibm10.hgr';
% filename_p = '../netlists/ibm10_add01.hgr';
% filename_p = '../netlists/ibm10_add05.hgr';
% filename_p = '../netlists/ibm10_del05a.hgr';
% filename = '../netlists/ibm18.hgr';
% filename_p = '../netlists/ibm18_add05.hgr';
% filename_p = '../netlists/ibm18_del05a.hgr';


%% Parameters
num_eigs = 10; % number of eigenvalues and eigenvectors to use for the perturbed solution
node_areas = 1; % all nodes have unit area
area_constraint = 0.45; % schmoo across all possible area splits (set this intentionally low so we have a minimum of 1e-3% of the nodes on one side and the rest on the other, and try all partitions between that and 1-1e5)
num_runs = 5;

%% Rerun num_runs times and store results
metrics_cell = cell(3,num_runs);
time_cell = cell(3,num_runs);
time_vecs = zeros(3, num_runs);
time_vecs_no_parse = zeros(3, num_runs);
eig_cell = cell(3, num_runs);

for rind = 1:num_runs
    fprintf('Beginning iteration %d/%d...\n', rind, num_runs)
    % Run EIG to find exact partitioning results for original and perturbed systems
    fprintf('\tPartitioning original graph...\n')
    [metrics, times, matrices, eigs] = eig_partitioner(filename,num_eigs,node_areas,area_constraint); % normal
    fprintf('\tPartitioning perturbed graph (exact)...\n')
    [metrics_pe, times_pe, matrices_pe, eigs_pe] = eig_partitioner(filename_p,num_eigs,node_areas,area_constraint); % exact perturbation, only need 2 eigenvalues here



    % Run the perturbed solver
    eigs_to_correct = 2;
    fprintf('\tPartitioning perturbed graph (approximate)...\n')
    [metrics_p, times_p, matrices_p, eigs_p] = eig_partitioner_perturbed(filename_p,matrices.laplacian,eigs.vals,eigs.vecs,node_areas,area_constraint,eigs_to_correct); % approx perturbation
   
    metrics_cell{1,rind} = metrics;
    metrics_cell{2,rind} = metrics_pe;
    metrics_cell{3,rind} = metrics_p;
    
    time_cell{1,rind} = times;
    time_cell{2,rind} = times_pe;
    time_cell{3,rind} = times_p;
    
    time_vecs(1,rind) = times.total;
    time_vecs(2,rind) = times_pe.total;
    time_vecs(3,rind) = times_p.total;
    
    time_vecs_no_parse(1,rind) = times.total - times.parse;
    time_vecs_no_parse(2,rind) = times_pe.total - times.parse;
    time_vecs_no_parse(3,rind) = times_p.total - times.parse;
    
    eig_cell{1,rind} = eigs;
    eig_cell{2,rind} = eigs_pe;
    eig_cell{3,rind} = eigs_p;
end

%% Plot all cutsize vectors

for rind = 1:num_runs
    time_vecs_no_parse(1,rind) = times.total - times.parse;
    time_vecs_no_parse(2,rind) = times_pe.total - times.parse;
    time_vecs_no_parse(3,rind) = times_p.total - times.parse;
end
%%
% Normal cutsize
figure(1)
clf
hold on
for rind = 1:num_runs
    plot(metrics_cell{1,rind}.skew,metrics_cell{1,rind}.cutsize.skew.vec,'k')
end
xlabel('Skew')
ylabel('Cutsize')
title('Original Graph')

% PE cutsize
figure(2)
clf
hold on
for rind = 1:num_runs
    plot(metrics_cell{2,rind}.skew,metrics_cell{2,rind}.cutsize.skew.vec,'k')
end
xlabel('Skew')
ylabel('Cutsize')
title('Perturbed Graph - Exact Solution')

% P cutsize
figure(3)
clf
hold on
for rind = 1:num_runs
    plot(metrics_cell{3,rind}.skew,metrics_cell{3,rind}.cutsize.skew.vec,'k')
end
xlabel('Skew')
ylabel('Cutsize')
title('Perturbed Graph - Approx Solution')

fixfigs(1:3,3,14,12);

%% Plot all ratio cut vectors
% Normal 
figure(4)
clf
hold on
for rind = 1:num_runs
    plot(metrics_cell{1,rind}.skew,metrics_cell{1,rind}.ratio_cut.skew.vec,'k')
end
xlabel('Skew')
ylabel('Cutsize')
title('Original Graph')

% PE 
figure(5)
clf
hold on
for rind = 1:num_runs
    plot(metrics_cell{2,rind}.skew,metrics_cell{2,rind}.ratio_cut.skew.vec,'k')
end
xlabel('Skew')
ylabel('Cutsize')
title('Perturbed Graph - Exact Solution')

% P 
figure(6)
clf
hold on
for rind = 1:num_runs
    plot(metrics_cell{3,rind}.skew,metrics_cell{3,rind}.ratio_cut.skew.vec,'k')
end
xlabel('Skew')
ylabel('Cutsize')
title('Perturbed Graph - Approx Solution')

fixfigs(4:6,3,14,12)

%% All Cutsizes

figure(7)
clf
hold on
for rind = 1:num_runs
    plot(metrics_cell{1,rind}.skew,metrics_cell{1,rind}.cutsize.skew.vec,'k')
    plot(metrics_cell{2,rind}.skew,metrics_cell{2,rind}.cutsize.skew.vec,'b')
    plot(metrics_cell{3,rind}.skew,metrics_cell{3,rind}.cutsize.skew.vec,'r')
end
xlabel('Skew')
ylabel('Cutsize')
title('Cutsize Multiple Times')
xlim([0, 0.05])
fixfigs(7,3,14,12)


%% All Ratio Cuts

figure(8)
clf
hold on
for rind = 1:num_runs
    plot(metrics_cell{1,rind}.skew,metrics_cell{1,rind}.ratio_cut.skew.vec,'k')
    plot(metrics_cell{2,rind}.skew,metrics_cell{2,rind}.ratio_cut.skew.vec,'b')
    plot(metrics_cell{3,rind}.skew,metrics_cell{3,rind}.ratio_cut.skew.vec,'r')
end
xlabel('Skew')
ylabel('Ratio Cut')
title('Ratio Cut Multiple Times')
xlim([0, 0.05])
fixfigs(8,3,14,12)


%% Times

figure(9)
clf
h1 = cdfplot(time_vecs(1,:));
hold on
h2 = cdfplot(time_vecs(2,:));
h3 = cdfplot(time_vecs(3,:));
set(h1,'color','k')
set(h2,'color','b')
set(h3,'color','r')
xlabel('time (s)')
ylabel('Cumulative distribution')
title('Time distribution')
fixfigs(9,3,14,12)
set(gca,'xscale','log')

figure(10)
clf
h1 = cdfplot(time_vecs_no_parse(1,:));
hold on
h2 = cdfplot(time_vecs_no_parse(2,:));
h3 = cdfplot(time_vecs_no_parse(3,:));
set(h1,'color','k')
set(h2,'color','b')
set(h3,'color','r')
xlabel('time (s)')
ylabel('Cumulative distribution')
title('Time distribution (parse time excluded)')
fixfigs(10,3,14,12)
set(gca,'xscale','log')

%% All Eigs

figure(11)
clf
hold on
for rind = 1:num_runs
    plot(abs(eig_cell{1,rind}.vals),'k')
    plot(abs(eig_cell{2,rind}.vals),'b')
    plot(abs(eig_cell{3,rind}.vals),'r')
end
ylabel('Eigenvalues')
set(gca,'yscale','log')
title('Eigenvalues multiple times')
fixfigs(11,3,14,12)

