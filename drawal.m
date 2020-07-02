function X= drawal(data,miu,phi,to,l,u,a,b)
%Draw the posterior sample of skewness parameter alpha
% Copyright (C) 2020 Pravitasari, et al

    f = @(x)exp((to-1)*(log(x-l)+log(u-x))-sum((exp(1/sqrt(2*pi*phi))*((data-miu)/phi))-((x+1)*sum(log(1+((exp(-(exp(1/sqrt(2*pi*phi))*((data-miu)/phi)))/x)))))));
    g = @(x)gsbetapdf(x,to,l,u);
    grnd = @()gamrnd(a,b);
    rng('default') % For reproducibility
    X = accrejrnd(f,g,grnd,1e-1000,100,1);
    X = mean(X);
end