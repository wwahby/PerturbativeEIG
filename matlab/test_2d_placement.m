% test eig_partitioner function
close all
% clear all

%% Choose input file
% filename indicates original netlist
% filename_p indicates perturbed netlist

% filename = '../netlists/fract.hgr';
% filename_p = '../netlists/fract_add05.hgr';
% filename_p = '../netlists/fract_del05a.hgr';
% filename = '../netlists/p1.hgr';
% filename_p = '../netlists/p1_block.hgr';
% filename_p = '../netlists/p1_add05.hgr';
% filename_p = '../netlists/p1_del05a.hgr';
% filename = '../netlists/structP.hgr';
% filename_p = '../netlists/structP_add05.hgr';
% filename_p = '../netlists/structP_del05a.hgr';
% filename = '../netlists/p2.hgr';
% filename_p = '../netlists/p2_block.hgr';
% filename_p = '../netlists/p2_add05.hgr';
% filename_p = '../netlists/p2_man.hgr';
% filename_p = '../netlists/p2_del05a.hgr';
% filename = '../netlists/biomedP.hgr';
% filename_p = '../netlists/biomedP_block.hgr';
% filename_p = '../netlists/biomedP_add05.hgr';
% filename_p = '../netlists/biomedP_man.hgr';
% filename_p = '../netlists/biomedP_del05a.hgr';
% filename_p = '../netlists/biomedP_add005.hgr';
% filename = '../netlists/industry2.hgr';
% filename_p = '../netlists/industry2_block.hgr';
% filename_p = '../netlists/industry2_add05.hgr';
% filename_p = '../netlists/industry2_del05a.hgr';
filename = '../netlists/ibm01.hgr';
filename_p = '../netlists/ibm01_block.hgr';
% filename_p = '../netlists/ibm01_add05.hgr';
% filename_p = '../netlists/ibm01_del05a.hgr';
% filename = '../netlists/industry3.hgr';
% filename_p = '../netlists/industry3_block.hgr';
% filename_p = '../netlists/industry3_add05.hgr';
% filename_p = '../netlists/industry3_del05a.hgr';
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



%% Run EIG to find exact partitioning results for original and perturbed systems
[metrics times matrices eigs] = eig_partitioner(filename,num_eigs,node_areas,area_constraint); % normal
[metrics_pe times_pe matrices_pe eigs_pe] = eig_partitioner(filename_p,num_eigs,node_areas,area_constraint); % exact perturbation, only need 2 eigenvalues here



%% Run the perturbed solver
eigs_to_correct = 2:3;
[metrics_p times_p matrices_p eigs_p] = eig_partitioner_perturbed(filename_p,matrices.laplacian,eigs.vals,eigs.vecs,node_areas,area_constraint,eigs_to_correct); % approx perturbation

%% 2D placement

%unperturbed
xnum = 2;
ynum = 3;
x = eigs.vecs(:,xnum);
y = eigs.vecs(:,ynum);

xpe = eigs_pe.vecs(:,xnum);
ype = eigs_pe.vecs(:,ynum);

xp = eigs_p.vecs(:,xnum);
yp = eigs_p.vecs(:,ynum);

figure(1)
clf
scatter(x,y,'k')
xlabel('x (au)')
ylabel('y (au)')
tstr = sprintf('%s - orig',filename);
title(tstr);

figure(2)
clf
scatter(xpe,ype,'b')
xlabel('x (au)')
ylabel('y (au)')
tstr = sprintf('%s - perturbed (exact)',filename);
title(tstr);

figure(3)
clf
scatter(xp,yp,'r')
xlabel('x (au)')
ylabel('y (au)')
tstr = sprintf('%s - perturbed (approx)',filename);
title(tstr);

fixfigs(1:3,3,14,12)

%% lengths

lsq_x = calc_l_squared(matrices.adjacency,eigs.vecs(:,xnum));
lsq_x_pe = calc_l_squared(matrices_pe.adjacency,eigs_pe.vecs(:,xnum));
lsq_x_p = calc_l_squared(matrices_p.adjacency,eigs_p.vecs(:,xnum));

lsq_y = calc_l_squared(matrices.adjacency,eigs.vecs(:,ynum));
lsq_y_pe = calc_l_squared(matrices_pe.adjacency,eigs_pe.vecs(:,ynum));
lsq_y_p = calc_l_squared(matrices_p.adjacency,eigs_p.vecs(:,ynum));

lsq = lsq_x + lsq_y;
lsq_pe = lsq_x_pe + lsq_y_pe;
lsq_p = lsq_x_p + lsq_y_p;

fprintf('\n%s\n',filename)
fprintf('Lsq:\t\tX\t\t\t\tY\t\t\t\tTOT\t\t\t\tL/orig\t\t\ttime\t\t\tt/t_orig \n======================================================================================================\n')
fprintf('Orig:\t\t%d\t%d\t%d\t    ---- \t\t%d\t    ----  \nPExact: \t%d\t%d\t%d\t%d\t%d\t%d \nPApprox: \t%d\t%d\t%d\t%d\t%d\t%d \n\n',lsq_x,lsq_y,lsq,times.total,lsq_x_pe,lsq_y_pe,lsq_pe,lsq_pe/lsq,times_pe.total,times_pe.total/times.total,lsq_x_p,lsq_y_p,lsq_p,lsq_p/lsq,times_p.total,times_p.total/times.total);
%% visualize eigs
% figure(3)
% clf
% hold all
% for i=1:length(eigs.vals)
%     plot(eigs.vecs(:,i));
% end
% 
% figure(4)
% clf
% hold all
% for i=1:length(eigs_pe.vals)
%     plot(eigs_pe.vecs(:,i));
% end



% figure(5)
% clf
% plot(eigs.vecs(:,2),'k');
% hold on
% plot(eigs_pe.vecs(:,2),'b');
% plot(eigs_p.vecs(:,2),'r--');
% 
% figure(6)
% clf
% plot(eigs.vecs(:,3),'k');
% hold on
% plot(eigs_pe.vecs(:,3),'b');
% plot(eigs_p.vecs(:,3),'r');


%% Skew plots (just show best result from 0% skew up to full skew, i.e. compare a 40/60 split and a 60/40 split and show only the best result for skew of 0.1)
% figure(8)
% clf
% plot(metrics.skew,metrics.cutsize.skew.vec,'k')
% hold on
% plot(metrics_pe.skew,metrics_pe.cutsize.skew.vec,'b')
% plot(metrics_p.skew,metrics_p.cutsize.skew.vec,'r')
% xlabel('skew')
% ylabel('cutsize');
% title('Cutsize')
% xlim([0 0.5])
% fixfigs(8,3,14,12)
% 
% figure(9)
% clf
% plot(metrics.skew,metrics.ratio_cut.skew.vec,'k')
% hold on
% plot(metrics_pe.skew,metrics_pe.ratio_cut.skew.vec,'b')
% plot(metrics_p.skew,metrics_p.ratio_cut.skew.vec,'r')
% xlabel('skew')
% ylabel('ratio cut');
% title('Ratio cut')
% xlim([0 0.5])
% fixfigs(9,3,14,12)