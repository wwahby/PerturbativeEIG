function sim_ratio = compare_partitions(p1a,p2a,p2b)

length_a = length(p2a);
length_b = length(p2b);

matchaa = sum(ismember(p1a,p2a));
matchab = sum(ismember(p1a,p2b));

sim_ratio_a = matchaa/length_a;
sim_ratio_b = matchab/length_b;

sim_ratio = max(sim_ratio_a,sim_ratio_b);