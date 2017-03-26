function IC = Check_IC_1(X, D, I, s)
% DESCRIPTION
% Computes IC_1 of genlasso
if nargin < 4
    error('The input is incomplete!');
else
    [~,p] = size(X);
    m = size(D,1);
    if size(D,2) ~= p, error('The number of columns of D and X should be equal!'); end
    if isempty(I) || max(max(I),length(I)) > m, error('Wrong size of true signal set!'); end
    if size(s,1) < size(s,2), s = s';end
    if size(s,2) > 1 || size(s,1) ~= m, error('D\beta^{\star} should be a vector with length m!'); end
    if sum(s(I)==0) >= 1, error('The value of true signal set of D\beta^{\star} should be non-zeros!');end
end
XTX = X' * X;
DI = D(I,:);
DJ = D;
DJ(I,:) = [];
[V,S,U] = svd(DJ);
ind = find(diag(S)>1e-12);
U(:,ind) = [];
V(:,ind) = [];
if (length(ind)==p)
    error('Full rank D_J');
end
A = U*pinv(U'*XTX*U)*U';
Omega_s = pinv(DJ')*((XTX*A-eye(p))*(DI'*s(I)));

cvx_begin quiet
    variables u(size(V,2)) t
    minimize  (t)
    subject to
        abs(V*u - Omega_s) <= t
cvx_end
IC = t;
end