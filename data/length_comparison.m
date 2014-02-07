
num_nodes_vec = [ 149         833        1952        3014        6514       12637       12752       15406 ];
lsq_vec = [ 0.2652    0.0472    0.0562    0.0860    0.0693    0.0031    0.0300    0.0108];
lsq_pe_vec = [ 0.2526    0.0779    0.1084    0.1380    0.1018    0.0357    0.0711    0.0932];
lsq_p_vec = [ 11.0539    1.3850    0.6925    0.2000    0.3517   46.0708    0.1645    0.1588];

figure(1)
clf
loglog(num_nodes_vec,lsq_vec,'k')
hold on
loglog(num_nodes_vec,lsq_pe_vec,'b')
loglog(num_nodes_vec,lsq_p_vec,'r')
xlabel('Number of nodes')
ylabel('squared wirelength (au)')
xlim([1e2 3e4])

fixfigs(1,3,14,12)