close all
clear all
%%
filename = 'netlists/industry2.hgr';
area_constraint = 0.5;
max_partition_level = 9;

%%
disp('Recursively bipartitioning netlist...')
[terminals blocks cuts] = recursive_bipartition_eig(filename,area_constraint,max_partition_level);
%%
disp('Bipartitioning complete!')

figure(1)
clf
hold on

med_num_nodes = zeros(1,length(blocks));
med_num_terminals_vec = zeros(1,length(blocks));
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
    
    med_num_nodes(bind) = median(num_nodes);
    med_num_terminals_vec(bind) = median(num_terminals);
end

xlabel('Number of gates')
ylabel('Number of terminals')
       
set(gca,'yscale','log')
set(gca,'xscale','log')

%%
figure(2)
clf
hold on
scatter(med_num_nodes,med_num_terminals_vec,'r','MarkerFaceColor','r')

xlabel('Number of gates')
ylabel('Number of terminals')
set(gca,'yscale','log')
set(gca,'xscale','log')
ylim([1e0 2*max(med_num_terminals_vec)])
xlim([1e0 2*max(med_num_nodes)]);
fixfigs(1:2,3,14,12)
%%

back_ind = 2;
[k p] = calc_rent_params(med_num_terminals_vec(end),med_num_nodes(end),med_num_terminals_vec(end-back_ind),med_num_nodes(end-back_ind));
k
p

%% estimate number of tsvs used for each partition level

tsvs_actual{1} = 0;
tsvs_estimated{1} = 0;
total_nodes = length(blocks{1}{1});
for cind = 2:length(cuts)
    num_tiers = length(cuts{cind});
    nt_tot = zeros(1,num_tiers);
    for i = 1:num_tiers
        nt_tot(i) = sum(sum( cuts{cind}(1:i,i:end) ));
    end
    tsvs_actual{cind} = nt_tot;
    
    [nt_max nt_tot nt_to nt_through Tacmat] = estimate_tsvs_required(total_nodes,num_tiers,k,p);
    tsvs_estimated{cind} = nt_tot;
end

%%
figure(3)
clf
set(gca,'ColorOrder',[0 0 1; 1 0 0]);
set(gca,'LineStyleOrder','-|--|-.|:');
hold all
for cind = 2:length(cuts)
    plot(tsvs_estimated{cind})
    plot(tsvs_actual{cind})
end
fixfigs(3,3,14,12)
