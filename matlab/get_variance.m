% test eig_partitioner function
%% Choose input file
% filename = 'ckt_f216.hgr';
% filename_p = 'ckt_f216_p.hgr';
% filename = 'fract.hgr';
% filename_p = 'fract_ptest.hgr';
% filename = 'p1.hgr';
% filename_p = 'p1_add05.hgr';
% filename = 'structP.hgr';
% filename_p = 'structP_add05.hgr';
% filename = 'p2.hgr';
% filename_p = 'p2_add05.hgr';
% filename = 'biomedP.hgr';
% filename_p = 'biomedP_add05.hgr';
% filename_p = 'biomedP_add005.hgr';
filename = 'industry2.hgr';
filename_p = 'industry2_add05.hgr';
% filename = 'ibm01.hgr';
% filename_p = 'ibm01_add05.hgr';
% filename = 'industry3.hgr';
% filename_p = 'industry3_add05.hgr';
% filename = 'ibm10.hgr';
% filename_p = 'ibm10_add05.hgr';
% filename = 'ibm18.hgr';
% filename_p = 'ibm18_add05.hgr';

%% Partitioning parameters
num_eigs = 7;
node_areas = 1; % all nodes have unit area
area_constraint = 1e-5; % schmoo across all possible area splits (set this intentionally low so we have a minimum of 1e-3% of the nodes on one side and the rest on the other, and try all partitions between that and 1-1e5)
iterations = 20;

%% Run a bunch of times

for i= 1:iterations
    
    repstr = sprintf('Iteration %d of %d',i,iterations);
    disp(repstr);
    [metrics times matrices eigs] = eig_partitioner(filename,num_eigs,node_areas,area_constraint); % normal
    [metrics_pe times_pe matrices_pe eigs_pe] = eig_partitioner(filename_p,2,node_areas,area_constraint); % exact perturbation

    Qp = matrices_pe.laplacian - matrices.laplacian;
    Ape = matrices_pe.adjacency;
    [metrics_p times_p matrices_p eigs_p] = eig_partitioner_perturbed(Qp,Ape,eigs.vals,eigs.vecs,node_areas,area_constraint); % approx perturbation
    
    if (i == 1) % initialize all our arrays
        rcv = zeros(iterations,length(metrics.ratio_cut.vec));
        cv = zeros(iterations,length(metrics.ratio_cut.vec));
        rcsv = zeros(iterations,length(metrics.ratio_cut.skew.vec));
        csv = zeros(iterations,length(metrics.cutsize.skew.vec));
        
        rcvp = zeros(iterations,length(metrics_p.ratio_cut.vec));
        cvp = zeros(iterations,length(metrics_p.ratio_cut.vec));
        rcsvp = zeros(iterations,length(metrics_p.ratio_cut.skew.vec));
        csvp = zeros(iterations,length(metrics_p.cutsize.skew.vec));
        
        rcvpe = zeros(iterations,length(metrics_pe.ratio_cut.vec));
        cvpe = zeros(iterations,length(metrics_pe.ratio_cut.vec));
        rcsvpe = zeros(iterations,length(metrics_pe.ratio_cut.skew.vec));
        csvpe = zeros(iterations,length(metrics_pe.cutsize.skew.vec));
    end
    
    rcv(i,:) = metrics.ratio_cut.vec;
    cv(i,:) = metrics.cutsize.vec;
    rcsv(i,:) = metrics.ratio_cut.skew.vec;
    csv(i,:) = metrics.cutsize.skew.vec;
    
    rcvp(i,:) = metrics_p.ratio_cut.vec;
    cvp(i,:) = metrics_p.cutsize.vec;
    rcsvp(i,:) = metrics_p.ratio_cut.skew.vec;
    csvp(i,:) = metrics_p.cutsize.skew.vec;
    
    rcvpe(i,:) = metrics_pe.ratio_cut.vec;
    cvpe(i,:) = metrics_pe.cutsize.vec;
    rcsvpe(i,:) = metrics_pe.ratio_cut.skew.vec;
    csvpe(i,:) = metrics_pe.cutsize.skew.vec;
end

%% Plot!

skew = metrics.skew;
figure(1)
clf
for i=1:iterations
    plot(skew,rcsvp(i,:))
    if i == 1
        hold on
    end
end
