% This script test how good is the hopping expansion approximation as 
% more terms are included


close all; clc; clear;
% The matrix is 4to4 lattice with m0=0 anc csw=0

load('LQCD_A1.mat')
%A1=A2; clear A2;
dimensions= size(A1);
n=dimensions(1);
m0 =0.0;%-0.1;
A1=A1+m0*speye(n,n);


k=1/(8+2*m0); %Hopping Parameter
B=(4+m0)*speye(n); %Diagonal 
H = 2*(B-A1); %Off-Diag

% Check if re-writting of D is correct
% D=B*( speye(n)-inv(B)*0.5*H);
% 
% dif=D-A1;
% norm(dif,"fro")




max=zeros(1,10);
counter=1;
for N_neu=10:10:100
    
    temp= speye(n);
    A=temp;
    for i=1:N_neu
        
        temp = temp*H*k;
        A = A + temp;
    end
    A=2*k*A;
    
    for j=1:5
        x=rand(n,1);
        temp = norm(A*A1*x-x)
        if temp>max(counter)
            max(counter)=temp;
        end
    end
    counter=counter+1;
end

save("max_norms.mat", "max")


semilogy( 10:10:100, max(:), "*r")
hold on
saveas(gcf,'1E2.png')
hold off