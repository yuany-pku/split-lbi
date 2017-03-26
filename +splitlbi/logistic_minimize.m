function beta = logistic_minimize(X,y,D,nu)
[N,n] = size(X);
C = -repmat(y,1,size(X,2)) .* X;
DtopD = D' * D;
cvx_precision best
           cvx_begin quiet
               variable s(n)
               minimize ( sum(log(1 + exp(C*s))) / N + s' * DtopD * s/ 2/ nu );
               
           cvx_end
beta = s;
end