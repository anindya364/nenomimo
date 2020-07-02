function y = gsbetapdf(x,a,l,u)
% p.d.f of generalized symmetrical beta(a,a,l,u)
y= (gamma(2*a)/gamma(a^2))*(((x-l)*(u-l))^(a-1))/((u-l)^(2*a-1));
end