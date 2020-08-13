function dy = dely(t,y) 
    dy=zeros(3,1);
    %dy = t/2;
    dy(1)=t./3;
    dy(2)=t-1;
    dy(3)=t-1/2.*exp(-t);
end

