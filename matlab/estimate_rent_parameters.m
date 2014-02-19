filename = 'netlists/structP.hgr';

[terminals blocks] = recursive_bipartition_eig(filename,0.4);

figure(1)
clf
hold on
% [FIX] median plot doesn't seem to work
med_num_nodes = zeros(1,length(blocks));
med_num_terminals = zeros(1,length(blocks));
for bind = 2:length(blocks)
    for nind = 1:length(blocks{bind})
        num_nodes = zeros(1,length(blocks{bind}) );
        num_terminals = zeros(1,length(blocks{bind}) );
        
        num_nodes(nind) = length(blocks{bind}{nind});
        num_terminals(nind) = terminals{bind}(nind);
        
        
        scatter(num_nodes,num_terminals,'b')
    end
    
    med_num_nodes(bind) = median(num_nodes);
    med_num_terminals(bind) = median(num_terminals);
end
        
set(gca,'yscale','log')
set(gca,'xscale','log')

figure(2)
clf
scatter(med_num_nodes,med_num_terminals)
set(gca,'yscale','log')
set(gca,'xscale','log')
