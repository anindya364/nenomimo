function X= drawphi(data,phi,miu,al,a,b)
%Draw the posterior sample of dispersion parameter phi
% Copyright (C) 2020 Pravitasari, et al

    f = @(x)exp(-((a+1)*log(x))-sum((exp(1/sqrt(2*pi*phi))*((data-miu)/x))-((al+1)*sum(log(1+((exp(-(exp(1/sqrt(2*pi*phi))*((data-miu)/x)))/al)))))));
    g = @(x)invgampdf(x,a,b);
    grnd = @()gamrnd(a,b);
    rng('default') % For reproducibility
    X = accrejrnd(f,g,grnd,1e-2000,100,1);
    X = mean(X);
end