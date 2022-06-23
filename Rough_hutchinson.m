function rough_tr = Rough_hutchinson (A, N_rough, func)
n=size(A,1);
measurements=zeros(1,N_rough);
for i=1:N_rough
    x = round(rand(n,1));
    x(x==0)=-1;
    measurements(i)= x'* func(A,x);     %inv_A*BIG_X; %%inv_A*X;
    
end
rough_tr = sum(measurements(:))/N_rough;

end