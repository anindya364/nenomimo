function y = invgampdf(x,r,s)
% p.d.f of invers gamma(r,s)
y= (s^r)/(gamma(r))*(x^(-r-1))*exp(-s/x);
end