close all; clc; clear;
% The matrix is 4to4 lattice with m0=0 anc csw=1
% its trace with m0=0.0 is 7.6514e+02 + 8.3632e-17i
% To find it, simply Use in= inv(A1); true_tr = trace(in)

global A1;
load('LQCD_A1.mat')
%A1=A2; clear A2;
dimensions= size(A1);
n=dimensions(1);
m0 =-0.1;
A1=A1+m0*speye(n,n);
global A1;

product = @(A,b) A*b;
invert = @(A,b)  gmresnomsg(A, b, 1E-10, 200)*b;  %Use alternatively M\b;
second_term = @(A,b)  (A*b - gmresnomsg(A1, b,1E-10, 200) ); 
%----------------Hutchinson Parameters----------
trace_tol= 1E-1;
N= 100000;
counter=1;
counters_1=[];
counters_2=[];
global min_iters;
min_iters= 5;
N_rough =5;
%----------------------------------------------


%----------------Hopping matrices----------
k=1/(8+2*m0); %Hopping Parameter
B=(4+m0)*speye(n); %Diagonal
H = 2*(B-A1); %Off-Diag
temp= speye(n);
A=temp;
%----------------------------------------------

t_value=1.96; % 95% confidence
aux=1;
for N_neu=50%10:10:1000
    
    %------- consider 10 more Hopping terms every iter    
    for j=1:10
        temp = temp*H*k;
        A = A + temp;
    end
    A=2*k*A;
    %%
    %load('A_50_terms.mat');
    
    
    
    %-----------------Hutchinson method for FIRST term-----------------------------
    rough_tr_1 = Rough_hutchinson(A, N_rough, product);
    [First_term_traces(aux), counters_1(aux)] = Hutchinson(A, N, trace_tol, rough_tr_1, product);    
    %-----------------------------------------------------------
    
     %-----------------Hutchinson method for Second term-----------------------------
    rough_tr_2 = Rough_hutchinson(A, N_rough, second_term);
    [Second_term_traces(aux), counters_2(aux)] = Hutchinson(A, N, trace_tol, rough_tr_2, second_term);    
    %-----------------------------------------------------------
    
    
    
    
    %-------------------------------------------------------------------
    fprintf("N_neu= %d \t counter_1= %d \t trace_1 %f \n", N_neu, counters_1(aux), real(First_term_traces(aux)));
    fprintf("N_neu= %d \t counter_2= %d \t trace_2 %f \n", N_neu, counters_2(aux), real(Second_term_traces(aux)));
    fprintf("N_neu= %d \t trace %f \n \n", N_neu, real(Second_term_traces(aux) - real(First_term_traces(aux))  ));
    
    aux=aux+1;
    A=A/(2*k);
    
end



save("First_term_traces.mat","First_term_traces");
save("n_ests_1.mat","counters_1");

save("Second_term_traces.mat","Second_term_traces");
save("n_ests_2.mat","counters_2");


Traces = First_term_traces - Second_term_traces;
save("Traces.mat", "Traces");
