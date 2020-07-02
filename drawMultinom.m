function X = drawMultinom(k)

% Draw size(k,2) samples from a Multinomial distribution 
% Copyright (C) 2020 Pravitasari, et al

k = cumsum(k);
kmax = max(max(k))+1;
u = repmat(rand(1,size(k,2)).*k(end,:), size(k,1), 1);
m = (u < k) .* (kmax-k);
X = argmax(m);
