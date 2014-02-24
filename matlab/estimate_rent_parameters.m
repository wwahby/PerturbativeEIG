close all
clear all
%%
filename = 'netlists/p2.hgr';
area_constraint = 0.5;
max_partition_level = 5;

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

%%
fixfigs(1:2,3,14,12)
