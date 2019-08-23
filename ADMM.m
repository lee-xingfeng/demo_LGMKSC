%\|Z\|_*+alpha*\|E\|_{21}+beta*\sum w_i \|K-K_i\|_F^2  s.t. K = KZ+E
%\|J\|_*+alpha*\|E\|_{21}+beta*\sum w_i \|K-K_i\|_F^2  s.t. K = KZ+E, Z = J

function [result,Z,E,obj,iter]= ADMM(H,y,alpha,beta, mu,rho)
nCluster=length(unique(y));
nn=length(y);
m = size(H,3);

K = zeros(nn);
g = ones(1,m)/m;
for j = 1:m
    K = K + g(j)*H(:,:,j);
end
I = eye(nn);
E = zeros(nn);
Z = eye(nn);
J = Z;
%args
epsilon = 1e-7; maxIter = 3000;  max_mu = 1e5;

Y1=zeros(nn);
Y2=zeros(nn);
obj=[];
tic
for iter = 1:maxIter
    Zold=Z;
    
    %Update Z
    % norm2(K-K*Z-E+Y1/mu)^2+norm2(Z-J+Y2/mu)^2
    Z= (I+K'*K)\(K'*(K-E+Y1/mu)+J-Y2/mu);
    %Z(find(Z<0))=0;
    
    % Update J
    M=Z+Y2/mu;
    [U, S, V] = svd(M, 'econ');
    diagS = diag(S);
    svp = length(find(diagS > 1/mu));
    diagS = max(0,diagS - 1/mu);
    if svp < 0.5 %svp = 0
        svp = 1;
    end
    J= U(:,1:svp)*diag(diagS(1:svp))*V(:,1:svp)';
    J(find(J<0))=0;%J=(J+J')/2;
    
    % Update E
    G = K - K*Z+Y1/mu;
    E = updateE21(G,mu,alpha);
    %E = updateE(G,mu,alpha);
    
    %Update K
    sumH=zeros(nn);
    sum1=0;
    for j=1:m
        temp(j)=norm(H(:,:,j)-K,'fro');
        omega(j)=1/(2*temp(j));
        sumH=sumH+omega(j)*H(:,:,j);
        sum1=sum1+omega(j);
    end
    K = (2*beta*sumH+mu*(E-Y1/mu)*(I-Z'))/(2*beta*sum1*eye(nn)+mu*(I-Z)*(I-Z'));
    K(find(K<0))=0;K=(K+K')/2;
    
    leq1 = K - K*Z - E;
    leq2 = Z - J;
    Y1 = Y1 + mu*leq1;
    Y2 = Y2 + mu*leq2;
    mu  = min(rho*mu,max_mu);
    %mu = 1.125*mu;
    obj(iter)=norm(Z-Zold,'fro');
    if((iter>5)&(norm(Z-Zold,'fro') < norm(Zold,'fro') * epsilon))
        break
    end
end
toc
Z = Z - diag(diag(Z));
Z = abs(Z);
Z=(Z+Z')/2;
for ij=1:10
    CKSym = BuildAdjacency(thrC(Z,0.7));
    grps = SpectralClustering(CKSym,nCluster);
    grps = bestMap(y,grps);
    res(ij,:) =  ClusteringMeasure(grps,y);
end
result(1,1)=mean(res(:,1));
result(2,1)=mean(res(:,2));
result(3,1)=mean(res(:,3));
result(1,2)=std(res(:,1));
result(2,2)=std(res(:,2));
result(3,2)=std(res(:,3));
end
