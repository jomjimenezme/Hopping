function [tr, counter ] = Hutchinson(A, N, trace_tol, rough_tr, func)
global min_iters;
n=size(A,1);
measurements=zeros(1,N);
counter=1;
stop = norm(rough_tr)*trace_tol;
for i=1:N
    x = round(rand(n,1));
    x(x==0) =-1;
    measurements(i) = x'* func(A,x);     %inv_A*BIG_X; %%inv_A*X;
    
    
    tr=sum(measurements)/i;
    v= sum((measurements- tr).*conj(measurements- tr))/i;
    
    RMSD = sqrt(v/i);
    
    
    counter=counter+1;
    if(i > min_iters && norm(RMSD) < stop)
        break
    end
    
    %RMSD
end


end