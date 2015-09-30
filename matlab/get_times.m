% test eig_partitioner function
filename = '../netlists/p2.hgr';
filename_p = '../netlists/p2_add05.hgr';
num_eigs = 7;
node_areas = 1; % all nodes have unit area
area_constraint = 0.55; % schmoo across all possible area splits (set this intentionally low so we have a minimum of 1e-3% of the nodes on one side and the rest on the other, and try all partitions between that and 1-1e5)


iterations = 20;

%% declarations
t_par = zeros(1,iterations);
t_place = zeros(1,iterations);
t_parse = zeros(1,iterations);
t_eigs = zeros(1,iterations);
t_tot = zeros(1,iterations);

tpe_par = zeros(1,iterations);
tpe_place = zeros(1,iterations);
tpe_parse = zeros(1,iterations);
tpe_eigs = zeros(1,iterations);
tpe_tot = zeros(1,iterations);

tp_par = zeros(1,iterations);
tp_place = zeros(1,iterations);
tp_parse = zeros(1,iterations);
tp_eigs = zeros(1,iterations);
tp_tot = zeros(1,iterations);
tp_tot_with_parse = zeros(1,iterations);

rc_vec_o = zeros(1,iterations);
mc_vec_o = zeros(1,iterations);
skew_rc_vec_o = zeros(1,iterations);
skew_mc_vec_o = zeros(1,iterations);

rc_vec_pe = zeros(1,iterations);
mc_vec_pe = zeros(1,iterations);
skew_rc_vec_pe = zeros(1,iterations);
skew_mc_vec_pe = zeros(1,iterations);

rc_vec_pa = zeros(1,iterations);
mc_vec_pa = zeros(1,iterations);
skew_rc_vec_pa = zeros(1,iterations);
skew_mc_vec_pa = zeros(1,iterations);

%% Run a bunch of times

for i= 1:iterations
    
    repstr = sprintf('\nIteration %d of %d\n===============',i,iterations);
    disp(repstr);
    %% Run EIG to find exact partitioning results for original and perturbed systems
    [metrics times matrices eigs] = eig_partitioner(filename,num_eigs,node_areas,area_constraint); % normal
    [metrics_pe times_pe matrices_pe eigs_pe] = eig_partitioner(filename_p,num_eigs,node_areas,area_constraint); % exact perturbation, only need 2 eigenvalues here


    %% Run the perturbed solver
    eigs_to_correct = 2;
    [metrics_p times_p matrices_p eigs_p] = eig_partitioner_perturbed(filename_p,matrices.laplacian,eigs.vals,eigs.vecs,node_areas,area_constraint,eigs_to_correct); % approx perturbation

    mc_vec_o(i) = metrics.cutsize.min;
    mc_vec_pe(i) = metrics_pe.cutsize.min;
    mc_vec_pa(i) = metrics_p.cutsize.min;
    
    rc_vec_o(i) = metrics.ratio_cut.min;
    rc_vec_pe(i) = metrics_pe.ratio_cut.min;
    rc_vec_pa(i) = metrics_p.ratio_cut.min;

    t_par(i) = times.partition;
    t_place(i) = times.placement;
    t_parse(i) = times.parse;
    t_eigs(i) = times.eigs;
    t_tot(i) = times.total;
    
    tpe_par(i) = times_pe.partition;
    tpe_place(i) = times_pe.placement;
    tpe_parse(i) = times_pe.parse;
    tpe_eigs(i) = times_pe.eigs;
    tpe_tot(i) = times_pe.total;
    
    tp_par(i) = times_p.partition;
    tp_place(i) = times_p.placement;
    tp_parse(i) = times_p.parse;
    tp_eigs(i) = times_p.eigs;
    tp_tot(i) = times_p.total;
    tp_tot_with_parse(i) = times_p.total + tpe_parse(i); % Perturbed version needs to parse new netlist
end

%% Plot!

figure(1)
clf
h1 = cdfplot(t_tot);
hold on
h2 = cdfplot(tpe_tot);
h3 = cdfplot(tp_tot);
set(h1,'color','k')
set(h2,'color','b')
set(h3,'color','r')
xlabel('time (s)')
ylabel('Cumulative distribution')
title('Time distribution (industry3, 20 runs, 5% skew)')
fixfigs(1,3,14,12)
set(gca,'xscale','log')

figure(2)
clf
h1 = cdfplot(mc_vec_o);
hold on
h2 = cdfplot(mc_vec_pe);
h3 = cdfplot(mc_vec_pa);
set(h1,'color','k')
set(h2,'color','b')
set(h3,'color','r')
xlabel('Mincut')
ylabel('Cumulative distribution')
fixfigs(2,3,14,12)
set(gca,'xscale','log')

figure(3)
clf
h1 = cdfplot(rc_vec_o);
hold on
h2 = cdfplot(rc_vec_pe);
h3 = cdfplot(rc_vec_pa);
set(h1,'color','k')
set(h2,'color','b')
set(h3,'color','r')
xlabel('Ratio Cut')
ylabel('Cumulative distribution')
fixfigs(3,3,14,12)
set(gca,'xscale','log')
