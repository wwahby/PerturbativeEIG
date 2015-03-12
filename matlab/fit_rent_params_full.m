function [k, p, kvec, pvec] = fit_rent_params_full(nodes,terminals)
% fits data in nodes (x) and terminals (y)  to y = k*x^p
% assumes all node/terminal data is available.
% removes the top "half" of the data at each iteration and calculates fits
%   half here means data points corresponding to nodes > med(nodes)
% chooses median fit from all the options
% This truncation is done to avoid being influenced by the io-limited
% region in the nodes/terminals relation

% Clear out any zero values so we don't blow up the fit
zeromask = ((nodes == 0) | (terminals == 0));
nodes = nodes(~zeromask);
terminals = terminals(~zeromask);

% Find fits for smaller subsets of the data and take median values, to eliminate impact of io-limited region
kvec = zeros(1,length(nodes)-1);
pvec = zeros(1,length(kvec));
% for f_ind = 1:length(nodes)-1
%     xv = nodes(f_ind:end);
%     yv = terminals(f_ind:end);
%     
%     c = fit(xv',yv','power1');
%     coeffs = coeffvalues(c);
%     kvec(f_ind) = coeffs(1);
%     pvec(f_ind) = coeffs(2);
% end

gen_ind = 0;
max_trunc = max(nodes);
trunc_cond = nodes < max_trunc;
while (sum(trunc_cond) > 1)
    gen_ind = gen_ind + 1;
    
    xv = nodes(trunc_cond);
    yv = terminals(trunc_cond);
    
    c = fit(xv',yv','power1');
    coeffs = coeffvalues(c);
    kvec(gen_ind) = coeffs(1);
    pvec(gen_ind) = coeffs(2);
    
    max_trunc = 0.5*max_trunc;
    trunc_cond = nodes < max_trunc;
end

kvec = kvec(kvec > 0);
pvec = pvec(kvec > 0);
    

% this doesn't necessarily pick same (k,p) pair. Should still be OK
p = median(pvec);
k = median(kvec);