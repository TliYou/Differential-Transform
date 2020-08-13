function yp = ddefun(t,y,ydel,ypdel) 
    yp = zeros(6,1); 
    %yp = 1 + y - 2*ydel^2 - ypdel;
    yp(1)=y(2);
    yp(2)=y(3);
    yp(3)=ypdel(3,1)*ydel(1,1)+(y(1))^(2/3)+ydel(5,3);
    yp(4)=y(5);
    yp(5)=y(6);
    yp(6)=1/2*ypdel(6,2)+ydel(5,2)*ydel(1,1);
end

