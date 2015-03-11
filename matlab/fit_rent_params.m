function [k, p, kvec, pvec] = fit_rent_params(nodes,terminals)
% fits data in nodes (x) and terminals (y)  to y = k*x^p

% Clear out any zero values so we don't blow up the fit
zeromask = ((nodes == 0) | (terminals == 0));
nodes = nodes(~zeromask);
terminals = terminals(~zeromask);

% Find fits for smaller subsets of the data and take median values, to eliminate impact of io-limited region
kvec = zeros(1,length(nodes)-1);
pvec = zeros(1,length(kvec));
for f_ind = 1:length(nodes)-1
    xv = nodes(f_ind:end);
    yv = terminals(f_ind:end);
    
    c = fit(xv',yv','power1');
    coeffs = coeffvalues(c);
    kvec(f_ind) = coeffs(1);
    pvec(f_ind) = coeffs(2);
end


% this doesn't necessarily pick same (k,p) pair. Should still be OK
p = median(pvec);
k = median(kvec);