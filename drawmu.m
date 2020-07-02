function X= drawmu(data,phi,al,omega,sigma)
%Draw the posterior sample of location parameter mu
% Copyright (C) 2020 Pravitasari, et al

    f = @(x)exp(-(((x-omega)^2)/2*sigma)-sum((exp(1/sqrt(2*pi*phi))*((data-x)/phi))-((al+1)*sum(log(1+((exp(-(exp(1/sqrt(2*pi*phi)))*((data-x)/phi)))/al))))));
    g = @(x)normpdf(x,omega,sigma);
    grnd = @()normrnd(omega,sigma);
    rng('default') % For reproducibility
    X = accrejrnd(f,g,grnd,1e-1000,100,1);
    X = mean(X);
end