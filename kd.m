function Y = kd(X)
%Kronecker Delta function
%   
Y = zeros(size(X));
Y(X == 0) = 1;

end

