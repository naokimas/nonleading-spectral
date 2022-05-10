function dy=double_well1DPRX(t,x,betaprx,alphaprx)
%global N r1 r2 r3 D alpha
global K 
r1=1;
r2=2;
r3=5;
%D=0.1;
alpha=1;

dy=-(x-r1).*(x-r2).*(x-r3)+K*alphaprx*x.^alpha;