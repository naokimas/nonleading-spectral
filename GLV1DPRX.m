function dy=GLV1DPRX(t,x,betaprx,alphaprx)
%global N r1 r2 r3 D alpha
global c K
alp1=0.5;
%c=14.9963+1;

dy=alp1.*x+(K*alphaprx-c)*betaprx*x.^2;