function [dyp] = delyp(t,y) 
    dyp=zeros(2,1);
    %dyp = t-pi;
    dyp(1)= t-2;
    dyp(2)= t./2;
    
end
