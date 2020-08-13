function [ trproduct ] = product(k,f,g)
%Differential Transform of a product
%
S=0;
for l=0:k
    S=S+f(l).*g(k-l);
end
trproduct=S;

end

