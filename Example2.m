syms x y;
syms k;
r=10; %order of polynomial approximating the solution, r>2
%interval [0,t_1]
U1=zeros(1,r+1); %initiation of a vector U, ATTENTION, indexes shifted by 1 compared to the paper (the same by B,H,F)
U1(1:3)=[1 1 1/2]; %loading of the initial conditions
U2=zeros(1,r+1);
U2(1:3)=[0 0 1];
F=@(k) gamma(2/3+1)./(gamma(k+1).*gamma(2/3-k+1)); %DT of the outec compound of F
FF=F(0:r-1); %vector of values for F
eqn1='(k+1)*(k+2)*(k+3)*x==exp(-2)*product(k,@(k) 1./factorial(k), @(k) U1(k+1).*(1/3)^k)+H+2*kd(k)-(-1)^k/factorial(k)'; % recursive equation, x is U1(k+4)
eqn2='(k+1)*(k+2)*(k+3)*(1-1/(2^(k+1)))*y==2*product(k,@(k) kd(k-1), @(k) U1(k+1).*(1/3)^k)-2/(3^k)*U1(k+1)';
B=zeros(r,r); %initiation of a matrix with Bells coefficients
B(1,1)=1; % loading of the first value (and the first row)
k=0; %first step (no need to calculate B)
H=1;
U1(k+4)=solve(eval(eqn1),x); %U1(4)
U2(k+4)=solve(eval(eqn2),y); %U2(4)
for k=1:r-3 %calculation of values U(5) to U(r+1)
  for l=2:k+1 % filling the (k+1)th row of Bells matrix in the columns 2 to the diagonal
  for i= 1:k+2-l
  B(k+1,l)=B(k+1,l)+((i*(l-1))/(k))*U1(i+1)*B(k+1-i,l-1);
  end
  end
  H=dot(FF,B(k+1,:)); %calculation of H(k) as the scalar product of LL and (k+1)th row of B
  U1(k+4)=solve(eval(eqn1),x); %next coefficient of U1
  U2(k+4)=solve(eval(eqn2),y); %next coefficient of U2
end
reseniu1=dot(U1,x.^(0:r));
reseniu2=dot(U2,x.^(0:r));
display(U1);
display(U2);
%
%comparison - first interval
tspan1=linspace(0,0.35,8);
hresu1=subs(reseniu1,x,tspan1) %u1 using DT evaluated at 8 points of the first interval
hresu2=subs(reseniu2,x,tspan1) %u2 using DT evaluated at 8 points of the first interval
tspan = [1e-20 0.621];
sol = ddensd(@ddefun,@dely,@delyp,@history,tspan);
tn = linspace(1e-20,0.35,8);
yn = deval(sol,tn);
hu1=yn(1,:); %u1 using ddensd evaluated at 8 points of the first interval
hu2=yn(4,:); %u2 using ddensd evaluated at 8 points of the first interval
%display(r)
%er1=hresu1-hu1
%er2=hresu2-hu2
%
%second interval
t1=eval(solve(x-1/2*exp(-x)));
V1=zeros(1,r+1); 
t1k=t1.^[0:r];
V1(1)=sum(t1k.*U1); %loading initial values for V1
coef1=[1:r];
V1(2)=sum(t1.^[0:r-1].*coef1.*U1(2:r+1));
coef2=factorial([0:r-2]+2)./(2.*factorial([0:r-2]));
V1(3)=sum(t1.^[0:r-2].*coef2.*U1(3:r+1));
V2=zeros(1,r+1);
V2(1)=sum(t1k.*U2); %loading initial values for V2
V2(2)=sum(t1.^[0:1:r-1].*coef1.*U2(2:r+1));
V2(3)=sum(t1.^[0:1:r-2].*coef2.*U2(3:r+1));
%
E2=@(k) t1.*kd(k)+kd(k-1)-((-1).^k).*exp(-t1)./(2.*factorial(k));
EE2=E2(0:r-1); 
G2=@(k) gamma(2/3+1)./(gamma(k+1).*gamma(2/3-k+1)).*(V1(1)).^(2/3-k); 
GG2=G2(0:r-1); 
LL=U2(2:r+1).*[1:r];
%
%eqn3='(k+1)*(k+2)*(k+3)*x==exp(-2)*product(k,@(k) exp(t1)./factorial(k), @(k) V1(k+1).*(1/3)^k)+H2+F2'; % recursive equation, x is U1(k+4)
eqn3='(k+1)*(k+2)*(k+3)*x==N1+H2+F2';
%eqn4='(k+1)*(k+2)*(k+3)*(1-1/(2^(k+1)))*y==2*product(k,@(k) ((t1-1).*kd(k)+kd(k-1)), @(k) V1(k+1).*(1/3)^k)';
eqn4='(k+1)*(k+2)*(k+3)*y==((1/2)^(k+1))*(k+1)*(k+2)*(k+3)*N2+N3';
B1=zeros(r,r); 
B2=zeros(r,r); 
B1(1,1)=1; 
B2(1,1)=1;
k=0; %first step (no need to calculate B)
H2=V1(1);
F2=0;
N1=exp(-2)*exp(t1)*sum((t1/3).^[0:r].*U1);
N2=sum(gamma([k+4:r+1])./(gamma([1:r-k-2]).*gamma(k+4)).*U2(k+4:r+1).*(t1/2).^[0:r-k-3]);
N3=2*(t1-1)*sum(U1.*(t1/3).^[0:r]);
V1(k+4)=solve(eval(eqn3),x); %V1(4)
V2(k+4)=solve(eval(eqn4),y); %V2(4)
for k=1:r-3 %calculation of values U(5) to U(r+1)
  for l=2:k+1 %filling (k+1)th row of Bells matrix in the columns 2 to the diagonal
  for i= 1:k+2-l
  B1(k+1,l)=B1(k+1,l)+((i*(l-1))/(k))*V1(i+1)*B1(k+1-i,l-1);
  B2(k+1,l)=B2(k+1,l)+((i*(l-1))/(k))*E2(i+1)*B2(k+1-i,l-1);
  end
  end
  H2=dot(GG2,B1(k+1,:)); 
  F2=dot(LL,B2(k+1,:));
  Z=0;
  W=0;
  for ll=0:k
      S=sum(gamma([k-ll+1:r+1])./(gamma([1:r+ll-k+1]).*gamma(k-ll+1)).*U1(k-ll+1:r+1).*(t1/3).^[0:r+ll-k]);
      Z=Z+exp(t1)/factorial(ll)*(1/3)^(k-ll)*S;
      W=W+((t1-1)*kd(ll)+kd(ll-1))*(1/3)^(k-ll)*S;
  end
  N1=exp(-2)*Z;
  N2=sum(gamma([k+4:r+1])./(gamma([1:r-k-2]).*gamma(k+4)).*U2(k+4:r+1).*(t1/2).^[0:r-k-3]);
  N3=2*W;
  V1(k+4)=solve(eval(eqn3),x); %next coefficient of V1
  V2(k+4)=solve(eval(eqn4),y); %next coefficient of V2
end
reseniv1=dot(V1,(x-t1).^(0:r));
reseniv2=dot(V2,(x-t1).^(0:r));
%comparison - second interval
tspan2=linspace(0.352,0.621,8);
hresv1=subs(reseniv1,x,tspan2) %u1 using DT evaluated at 8 points of the second interval
hresv2=subs(reseniv2,x,tspan2) %u2 using DT evaluated at 8 points of the second interval
%tspan3 = [0.352 0.621];
%sol3 = ddensd(@ddefun,@dely,@delyp,@history,tspan3);
%ttn = linspace(0.352,0.621,8);
yyn = deval(sol,tspan2);
hv1=yyn(1,:); %u1/v1 using ddensd evaluated at 8 points of the second interval
hv2=yyn(4,:); %u2/v2 using ddensd evaluated at 8 points of the second interval
%display(r)
er3=hresv1-hv1;
er4=hresv2-hv2;
%t1=linspace(1e-20,0.3517,100);
%t2=linspace(0.352,0.621,100);
%y1=subs(reseniu1,x,t1);
%y2=subs(reseniv1,x,t2);
%yy1=subs(reseniu2,x,t1);
%yy2=subs(reseniv2,x,t2);
%tvalues=[t1 t2];
%yvalues=[y1 y2];
%yyvalues=[yy1 yy2];
%zn = deval(sol,t1);
%zzn=deval(sol3,t2);
%plot(tvalues,yvalues,tvalues,yyvalues);

