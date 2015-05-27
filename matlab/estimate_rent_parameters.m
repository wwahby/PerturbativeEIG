close all
clear all
%%
%filename = 'netlists/synnet_no1_2.hgr';
filename = 'netlists/mesh_4d_open_10.hgr';
area_constraint = 0.5;
max_partition_level = 8;

%%
disp('Recursively bipartitioning netlist...')
[terminals blocks cuts] = recursive_bipartition_eig(filename,area_constraint,max_partition_level);
disp('Bipartitioning complete!')

%%



med_num_nodes = zeros(1,length(blocks));
med_num_terminals_vec = zeros(1,length(blocks));
num_nodes_long = [];
num_terminals_long = [];

for bind = 2:length(blocks)
    
    num_nodes = zeros(1,length(blocks{bind}) );
    num_terminals = zeros(1,length(blocks{bind}) );
    
    
    
    for nind = 1:length(blocks{bind})
        
        num_nodes(nind) = length(blocks{bind}{nind});
        num_terminals(nind) = terminals{bind}(nind);
        med_num_terminals = median(terminals{bind}(nind));
        
        %scatter(num_nodes,num_terminals,'b','MarkerFaceColor','b')
        %scatter(num_nodes(nind),med_num_terminals,'r','MarkerFaceColor','r')
    end
    
    num_nodes_long = [num_nodes_long num_nodes];
    num_terminals_long = [num_terminals_long num_terminals];
    
    med_num_nodes(bind) = median(num_nodes);
    med_num_terminals_vec(bind) = median(num_terminals);
end

%% Estimate Rent Parameters

back_ind = 2;
[ks ps] = calc_rent_params_simple(med_num_terminals_vec(end),med_num_nodes(end),med_num_terminals_vec(end-back_ind),med_num_nodes(end-back_ind));
ks
ps

[kf pf kvec pvec] = fit_rent_params(med_num_nodes, med_num_terminals_vec);
kf
pf

trunc_min = num_nodes_long > 1*min(num_nodes_long);
trunc_max = num_nodes_long < 1/8*max(num_nodes_long);
trunc_cond = trunc_min & trunc_max;
num_nodes_long_trunc = num_nodes_long(trunc_cond);
num_terminals_long_trunc = num_terminals_long(trunc_cond);

[kt pt ktvec ptvec] = fit_rent_params_full(num_nodes_long, num_terminals_long);
kt
pt

k = kt;
p = pt;

%% estimate number of tsvs used for each partition level

%cuts_orig = cuts;
tsvs_actual{1} = 0;
tsvs_estimated{1} = 0;
tsvs_to{1} = 0;
tsvs_to_2{1} = 0;
tsvs_actual_tot = zeros(1,length(cuts));
tsvs_estimated_to = zeros(1,length(cuts));
tsvs_estimated_tot = zeros(1,length(cuts));
tsvs_actual_tot_2 = zeros(1,length(cuts));
tsvs_estimated_to_2 = zeros(1,length(cuts));
tsvs_estimated_tot_2 = zeros(1,length(cuts));
total_nodes = length(blocks{1}{1});
num_tiers_vec = ones(1,length(cuts));
for cind = 2:length(cuts)
    num_tiers = length(cuts{cind});
    num_tiers_vec(cind) = num_tiers;
    
%     % EIG partitioning to get 1D placement of tiers
%     tier_totals = sum(cuts{cind});
%     tier_weight = 1./(tier_totals-1);
%     
%     D = diag(tier_totals);
%     A = zeros(num_tiers,num_tiers);
%     for i =1:num_tiers
%         A(i,:) = cuts{cind}(i,:).*tier_weight;
%     end
%     Q = D - A;
%     [evecs evals] = eig(Q);
%     if(length(evals) > 1)
%         [y tier_order] = sort(evecs(:,2));
%         cuts{cind} = cuts{cind}(tier_order,tier_order); % reorder tiers according to 
%     end

    nt_tot = zeros(1,num_tiers);
    for i = 1:num_tiers
        nt_tot(i) = sum(sum( cuts{cind}(1:i,i:end) ));
    end
    tsvs_actual{cind} = nt_tot;
    tsvs_actual_tot(cind) = sum(nt_tot);
    
    [nt_max, nt_tot, nt_to, nt_through, Tacmat] = estimate_tsvs_required(total_nodes,num_tiers,kf,pf);
    [nt_max_2, nt_tot_2, nt_to_2, nt_through_2, Tacmat_2] = estimate_tsvs_required(total_nodes,num_tiers,kt,pt);
    tsvs_estimated{cind} = nt_tot;
    tsvs_to{cind} = nt_to;
    tsvs_estimated_tot(cind) = sum(nt_tot);
    tsvs_estimated_to(cind) = sum(nt_to);
    
    tsvs_estimated_2{cind} = nt_tot_2;
    tsvs_to_2{cind} = nt_to_2;
    tsvs_estimated_tot_2(cind) = sum(nt_tot_2);
    tsvs_estimated_to_2(cind) = sum(nt_to_2);
end




%% Plot the bipartitioning results and rent fit

% Create fit lines to illustrate the fits chosen
ts = ks .* med_num_nodes.^ps;
tf = kf .* med_num_nodes.^pf;
tt = kt .* med_num_nodes.^pt;

% Terminal/Node relation and fits
figure(1)
clf
hold on
scatter(num_nodes_long,num_terminals_long,'b','MarkerFaceColor','b')
scatter(med_num_nodes,med_num_terminals_vec,'r','MarkerFaceColor','r')
plot(med_num_nodes,ts,'k-')
plot(med_num_nodes,tf,'g:')
plot(med_num_nodes,tt,'c:')
xlabel('Number of gates')
ylabel('Number of terminals')
set(gca,'yscale','log')
set(gca,'xscale','log')
ylim([0.2*min(num_terminals_long) 5*max(num_terminals_long)])
%xlim([1e0 2*max(med_num_nodes)]);
fixfigs(1,3,14,12)


% Just the median terminal/node relation
figure(2)
clf
hold on
scatter(med_num_nodes,med_num_terminals_vec,'r','MarkerFaceColor','r')
xlabel('Number of gates')
ylabel('Number of terminals')
set(gca,'yscale','log')
set(gca,'xscale','log')
%ylim([0.2*min(med_num_terminals_vec) 5*max(med_num_terminals_vec)])
%xlim([0.2*min(med_num_nodes) 5*max(med_num_nodes)]);
fixfigs(2,3,14,12)



%% Additional plots

% TSV/MIV Prediction vs actual cutsize at each step
figure(3)
clf
hold on
cind = 5;
plot(tsvs_actual{cind},'k-')
plot(tsvs_to{cind},'r-')
plot(tsvs_estimated{cind},'b-')

xlabel('Layer number')
ylabel('Number of TSVs')
%ylim([0 400])
fixfigs(3,3,14,12)

%% Comparison of Via estimation methods (tot, vs to-only) and actual cutsize
figure(5)
clf
plot(num_tiers_vec,tsvs_estimated_tot,'b')
hold on
%plot(num_tiers_vec,tsvs_estimated_tot_2,'g')
plot(num_tiers_vec,tsvs_estimated_to,'r')
plot(num_tiers_vec,tsvs_actual_tot,'k')
xlim([2 16])
xlabel('Number of tiers')
ylabel('Total TSVs in design')
%set(gca,'yscale','log')
set(gca,'xscale','log')
fixfigs(5,3,14,12)

%% Rent param fit info
figure(6)
clf
scatter(kvec,pvec)
% set(gca,'yscale','log')
% set(gca,'xscale','log')
xlabel('k (fit)')
ylabel('p (fit)')
fixfigs(6,3,14,12)

figure(7)
clf
plot(pvec)
xlabel('Fit case')
ylabel('p (fit)')
fixfigs(7,3,14,12)

figure(8)
clf
plot(kvec)
set(gca,'yscale','log')
fixfigs(8,3,14,12)
xlabel('Fit case')
ylabel('k (fit)')