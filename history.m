function y = history(t)
    %y = cos(t);
    y=zeros(6,1);
    y(1)=exp(t);
    y(2)=exp(t);
    y(3)=exp(t);
    y(4)=t^2;
    y(5)=2*t;
    y(6)=2;
end