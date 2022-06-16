function mems2nuc=findingmem2nuc(kbase,basedist)
%% Scales Membrane-Nucleus Connection way up if they get too close. 
syms k x
a=-5;         %-2.1

f(x)=abs(k)*(exp(a*x))+k;


mems2nuc=double(solve(f(basedist)==kbase,k));

 f(x)=subs(f(x),k,mems2nuc);

% fplot(f(x),[0.00000001 15])
%  
%  double(f(0.00000001))
%  double(f(15))
%  f;
end