function x = gmresnomsg(A,b,tol, max_iters)
    [x,~] = gmres(A, b, [], tol, max_iters); 
end