function irr = Check_irr(X, D, I, nu)
% DESCRIPTION
% This file compute IRR(\nu).
if nargin < 4
    error('The input is incomplete!');
else
    [~,p] = size(X);
    m = size(D,1);
    if size(D,2) ~= p, error('The number of columns of D and X should be equal!'); end
    if length(I) == 0 || max(max(I),length(I)) > m, error('Wrong size of true signal set!'); end
    if nu <= 0, error('nu should be greater than 0!'); end
end
[m,~] = size(D);
n = size(X,1);
J = 1:m;
J(I) = [];
Sigma = eye(m) - D * pinv(nu * X'*X/n + D'*D) * D';
irr = norm(Sigma(J,I) * inv(Sigma(I,I)),'inf');
end
