syms x;
syms k;
r=25; %order of polynomial aproximating the solution, r>1
%first interval [1,e], center of U at 1
U=zeros(1,r+1); %initiation of a vector U, ATTENTION, indexes shifted by 1 compared to the paper (the same by B,H,F)
U(1)=exp(1); %loading of the initial condition
F=@(k)((-1).^(k+1))./k; %DT of logarithm at 1
FF=[0 F(1:r-1)]; %vector of values for F, the first place according to the value at t0
E=@(k) 1./factorial(k); %DT of exponential function at 0
EE=[1 E(1:r-1)]; %vector of values for E, the first place according to the value at 0
eqn1='(k+1)*x==U(k+1)-kd(k-1)-kd(k)+H'; % recursive equation, x at place of tranformed the highest derivative - here U(k+2)
B=zeros(r,r); %initiation of a matrix with Bells coefficients
B(1,1)=1; % loading of the first value (and the first row)
k=0; %first step (no need to calculate B)
H=1; %H[1]=E[0]
U(k+2)=solve(eval(eqn1),x); %U(2)
for k=1:r-1 %calculation of values U(3) to U(r+1)
  for l=2:k+1 % filling the (k+1)th row of Bells matrix in the columns 2 to the diagonal
  for i= 1:k+2-l
  B(k+1,l)=B(k+1,l)+((i*(l-1))/(k))*FF(i+1)*B(k+1-i,l-1);
  end
  end
  H=dot(EE,B(k+1,:)); %calculation of H(k) as the scalar product of vector of values EE and (k+1)th row of matrix B
  U(k+2)=solve(eval(eqn1),x); %next coefficient of the vector U
end
reseni1=dot(U,(x-1).^(0:r)); %approx. solution in the form of a polynomial of variable x
%
%
%second interval [e,e^e], center UU at e
UU=zeros(1,r+1); %initiation of a vector UU, ATTENTION, indexes shifted by 1 compared to the paper 
UU(1)=(exp(1))^(exp(1)); %loading the initial condition - or the correct polyval(reseni1,exp(1))
L=@(k)((-1).^(k+1))./(k.*(exp(k))); %DT pf logarithm at e
LL=[1 L(1:r-1)]; %vector of values for L, the first place according to value at e
eqn2='(k+1)*x==UU(k+1)-kd(k-1)-exp(1)*kd(k)+H'; % recursive equation, x is UU(k+2)
BB=zeros(r,r); %initiation of a matrix with Bells coefficients
BB(1,1)=1; %loading of the first value (and the first row)
k=0; %first step (no need to calculate B)
H=exp(1); %H[1]=U(0)[1]
UU(k+2)=solve(eval(eqn2),x); %U(2)
for k=1:r-1 %calculation of values U(3) to U(r+1)
  for l=2:k+1 % filling the (k+1)th row of Bells matrix in the columns 2 to the diagonal
  for i= 1:k+2-l
  BB(k+1,l)=BB(k+1,l)+((i*(l-1))/(k))*LL(i+1)*BB(k+1-i,l-1);
  end
  end
  H=dot(U(1:r),BB(k+1,:)); %calculation of H(k) as the scalar product of vector of values U and (k+1)th row of matrix B
  UU(k+2)=solve(eval(eqn2),x); %next coefficient of the vector UU
end
reseni2=dot(UU,(x-exp(1)).^(0:r)); %polynomial form of the approx. solution (highest power to lowest)
%
%
%comparison with exact solution
presne=exp(x);
display(r)
%first interval
tspan1=linspace(1,exp(1),5);
presnehodnoty1=subs(presne,x,tspan1);
hodnotyreseni1=subs(reseni1,x,tspan1);
error1=presnehodnoty1-hodnotyreseni1
%second interval
tspan2=linspace(exp(1),(exp(1))^(exp(1)),5)
presnehodnoty2=subs(presne,x,tspan2);
hodnotyreseni2=subs(reseni2,x,tspan2);
error2=presnehodnoty2-hodnotyreseni2;

tspan = [1 exp(1)];
sol1 = ddensd(@ddefun2,@dely2,@delyp2,@history2,tspan);
yn1 = deval(sol1,tspan1);
e1=presnehodnoty1-yn1;
tspan = [exp(1) (exp(1))^(exp(1))];
sol2 = ddensd(@ddefun2,@dely2,@delyp2,@history2,tspan);
yn2 = deval(sol2,tspan2);
e2=presnehodnoty2-yn2;