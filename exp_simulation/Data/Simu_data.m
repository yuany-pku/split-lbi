% DESCRIPTION
% Generate 100 times X and y for simulation experiments %
run ../Linear_setting_sim;
for i = 1:100
    rng(i)
    X = randn(n,p);
    y = X * beta_true + randn(n,1) * sigma;
    save(strcat('X',num2str(i),'_lbi_linear.txt'),'X','-ascii');
    save(strcat('y',num2str(i),'_lbi_linear.txt'),'y','-ascii');
end

