function beta = linear_minimize(X,y)
n = size(X,1);
ytopX = y' * X;
XtopX = X' * X;
beta = pinv(full(XtopX)) * (ytopX' / n);
end