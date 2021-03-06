function beta = linear_minimize(X,y,D,nu)
[N,n] = size(X);
ytopX = y' * X;
DtopD = D' * D;
XtopX = X' * X;
beta = pinv(full(XtopX / N + DtopD / nu)) * (ytopX' / N);
end