function dy=GLV1(t,x)
%function dy=GLV1(t,x)
%global N r1 r2 r3 D alpha
global A c K
alp1=0.5;
%K=0.2;
dy=zeros(length(A),1);
%L=A-diag(sum(A));
%coup=A*x;
coup=A*x;
dy=alp1.*x-c*x.*x+K*x.*coup;