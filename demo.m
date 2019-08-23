clc;
warning off
ds = {'YALE_165n_1024d_15c_zscore_uni'};

for iData = 1
    dataset = ds{iData}
    data_file = fullfile([dataset, '.mat']);
    kernel_file = fullfile([dataset, '_allkernel.mat']);
    
    load(data_file)
    load(kernel_file)
    
    alpha=1;
    beta=20;
    mu=0.5;
    rho = 1.1;
    [result,Z,E,obj,iter]= ADMM(K,y,alpha,beta, mu,rho);
    %fprintf('#%.3f %.3f %.3f %.3f',alpha,beta,gamma);
    fprintf('# rank=%.2f> MEA-ACC:%.4f   MEA-NMI:%.4f   MEA-Purity:%.4f\n', sum(svd(Z)), result(1,1), result(2,1),result(3,1));
    fprintf('# %d> STD-ACC:%.4f   STD-NMI:%.4f   STD-Purity:%.4f\n', rank(Z), result(1,2), result(2,2),result(3,2));
end
