function X = accrejrnd(f,g,grnd,c,m,n)
%Function for generate samples with Acceptance Rejection

    X = zeros(m,n); % Preallocate memory
    for i = 1:m*n
        accept = false;
        while accept == false
            u = rand();
            v = grnd();
            if c*u <= f(v)/g(v)
               X(i) = v;
               accept = true;
            end
        end
    end
end
