function v = Gamma(n,a,b)
    mu = a/b;
    sig = sqrt(a/(b*b));
    v = zeros(n,1);
    x = a/b;
    v(1) = x;
        for i = 2:n
             can = normrnd(mu,sig,1);
             accept = min(1,(gampdf(can,a,b)/gampdf(x,a,b))/(normpdf(can,mu,sig)/normpdf(x, mu, sig)));
             u = rand;
             if(u<accept)
                 x = can;
             end
             v(i) = x;
        end
        v = v.^-1;
end