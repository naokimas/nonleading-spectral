function dy=double_well_work2(t,x)
%global N r1 r2 r3 D alpha
global A K
r1=1;
r2=2;
r3=5;
%D=0.1;
alpha=1;
%c=0;
dy=zeros(length(A),1);
%L=A-diag(sum(A));
%r1
dy=-(x-r1).*(x-r2).*(x-r3)+K.*A*x.^alpha;