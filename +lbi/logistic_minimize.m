function beta = logistic_minimize(X,y)
[~,n] = size(X);
C = -repmat(y,1,size(X,2)) .* X;
cvx_precision best
           cvx_begin quiet
               variable s(n)
               minimize ( sum(log(1 + exp(C*s))));
               
           cvx_end
beta = s;
end