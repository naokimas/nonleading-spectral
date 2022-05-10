function dy=SIS1(t,x)
%global N r1 r2 r3 D alpha
global A K
D=1;
dy=zeros(length(A),1);
%L=A-diag(sum(A));
%r1
coup=A*x;
dy=-D*x+K*(1-x).*coup;