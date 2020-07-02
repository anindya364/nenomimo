function [miu, phi, al, p, z, churn] = ...
    gibbsMSNBURRMM(Y, K, alpha, Nsamp,z)

% Use Markov chain Monte Carlo simulation to cluster the data Y into a
% mixture of K univariate MSNBurr.
% Outputs of function are samples from the posterior distributions, so
% that theta(i) = [mu(i,:) phi(i,:) al(i,:) z ], i = 1..Nsamp
% Copyright (C) 2020 Pravitasari, et al

N = length(Y);
%z = drawMultinom(ones(k,N));
for j=1:K
  yj = Y(find(z == j));
  miu(1,j) = mean(yj);
  phi(1,j) = std(yj);
  al(1,j)=1;

end
p(1,:) = full(sparse(1, z, 1, 1, K));

% Go! draw the parameters
for i=2:Nsamp
  % Mu
  for j=1:K
    n = sum(z == j);
    yj = Y(find(z == j));
        
    phi_temp=phi(i-1,j);
    al_temp=al(i-1,j);
    omega_temp=mean(yj);
    sigma=std(yj);
    miu(i,j)=drawmu(yj,phi_temp,al_temp,omega_temp,sigma);
  end

  % Sigma
  for j=1:K
    n = sum(z == j);
    yj = Y(find(z == j));
    phi_temp=phi(i-1,j);
    miu_temp=miu(i-1,j);
    al_temp=al(i-1,j);
    a=2;
    b=std(yj);
    phi(i,j)=drawphi(yj,phi_temp,miu_temp,al_temp,a,b);
  end
  
  % Alpha
  for j=1:K
    n = sum(z == j);
    yj = Y(find(z == j));
    
    miu_temp=miu(i-1,j);
    phi_temp=phi(i-1,j);
    to=2;
    l=1;
    u=2;
    a=1;
    b=std(yj)/sqrt(a);
    al(i,j)=drawal(yj,miu_temp,phi_temp,to,l,u,a,b);
    
  end
  
  for j=1:K
    tmp_pr(j,:) = normalLike(Y, miu(i-1,j), phi(i-1,j));
  end
  n = tabulate(z);
  n = n(:,2)';

  % Scale likelihoods by class memberships times prior
  pri = repmat((n'+alpha/K)/(sum(n)-1+alpha), 1, N);
  idxs = sub2ind(size(pri), z, [1:N]);
  pri(idxs) = pri(idxs) - 1/(sum(n)-1+alpha);
  
  tz = drawMultinom(pri .* tmp_pr);
  churn(i) = sum(tz ~= z);
  z = tz;
  p(i,:) = n;
  
end


function x = drawNormal(mu, sigSq)
% Draw one sample from a Gaussian with mean mu and variance sigSq
x = randn(1)*sqrt(sigSq) + mu;


function pr = normalLike(y, mu, sigSq)
% Evaluate the likelihood of the points y under the Gaussian with mean
% mu and variance sigSq
pr = 1/sqrt(2*pi*sigSq) .* exp(-(y-mu).^2/(2*sigSq));


