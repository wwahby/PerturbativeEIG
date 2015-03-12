function [k p] = calc_rent_params_simple(T1,N1,T2,N2)

p = log(T1/T2)/log(N1/N2);

k = sqrt(T1*T2/(N1*N2)^p);