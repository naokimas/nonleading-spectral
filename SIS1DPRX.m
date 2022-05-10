function dy=SIS1DPRX(t,x,betaprx,alphaprx)
%global N r1 r2 r3 D alpha
global K
D=1;
%K=2; %0.2
dy=zeros(1,1);
%L=A-diag(sum(A));
%r1

dy=-D*x+K*(1-betaprx*x).*alphaprx.*x;