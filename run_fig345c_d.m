
 global K A

t0 = 0; % the start time 
tf = 500; % the end time 
x_low = .01;%initial condition
%x_hi=10;
%generate networks save edgelist in Ei.tex file, i is the number of network 
%%save the eigen vectors corresponding to largest eigen value and the
%%optimal eigen value as the coulumns in ai.tex

network_data=[];
for i=1:100   %%number of networks to be considered ## 1 for single network result. 
    network_edge=sprintf('E%d.txt',i); %%%%%loading the edgelist of network
    edgelist=load(network_edge);
    
    network_ev=sprintf('a%d.txt',i); %%%loading the eigen vectors
    ev=load(network_ev);
 
    metad = sprintf('metadata%d.txt',i); %%%reading the storred eigen values
    fileID = fopen(metad);

format long
cell_data= textscan(fileID,'%s%f','Delimiter','=','headerLines',0);
parameters=cell_data{2};
beta_maxev=parameters(1);
betaprxN1=parameters(2);
e_val_max=parameters(7);
alphaprxN1=parameters(8);
alphaprxN2=parameters(9);

    
    A=edgeL2adj(edgelist); %%%% Converting edgelist to a adjacency matrix
    
   


n = length(A); % number of nodes in the network
x0l = ones(n,1)*x_low;  % Initial conditions
%x0h = ones(n,1)*x_hi;
options = [];
%data_gao=[];
data_laurence=[];
%[e_vec, e_val]=eig(A);
deg1=sum(A);
betaprx=beta_maxev; %% beta for maximum eigen value

alphaprx=e_val_max; %% Maximum eigen value



for K=0.:0.2:4
    [t,x] = ode45(@SIS1,[t0,tf],x0l,options); %%solving the system for equilibrium
    xl_ss = x(end,:); %%%equilibribum 
R_low_prx=sum(ev(:,1).*xl_ss'); %%%%Calculating  R using eigen vector corresponding to largest eigen value
R_low_n2=sum(ev(:,3).*xl_ss'); %%Calculating  R using eigen vector corresponding to optimal eigen value
%R_hi_naoki2=sum(a(:,3).*xh_ss');
[t,xll1]=ode45(@SIS1DPRX,[0,200],x_low,options,betaprx,alphaprx); %% solving One dimensionl model corresponding to largest eigen value and corresponding beta
[t,xlN1]=ode45(@SIS1DPRX,[0,200],x_low,options,1,alphaprx);%% solving One dimensionl model corresponding to largest eigen value and  beta=1
[t,xlN2]=ode45(@SIS1DPRX,[0,200],x_low,options,1,alphaprxN2);  %% solving One dimensionl model corresponding to optimal eigen value and beta=1
%data_laurence=[data_laurence;D R_low_prx R_hi_prx xll1(end) xhl1(end) R_low_naoki1 R_hi_naoki1 R_low_naoki2 R_hi_naoki2];
data_laurence=[data_laurence;K xl_nn alphaprx R_low_prx xll1(end) alphaprxN1 R_low_prx xlN1(end) alphaprxN2 R_low_n2 xlN2(end) ];

end
    %save datadouble_wellA3E6.mat data_gao data_laurence
    network_data=[network_data;i*ones(length(data_laurence(:,1)),1),data_laurence];
end
%save network_Double_well162datak4.mat network_data

l1=5;l1D=6;
n1=8;n1D=9;
n2=11;n2D=12;
mid_val=20; %%
Df=41; %%row number corresponding to final couplinbg strength considered for each network
m=41;n=size(network_data,1)/m;
k=network_data(1:m,2);

laurence_error=abs(network_data(:,l1)-network_data(:,l1D));
laurence_error_mat=reshape(laurence_error,m,n);
laurence_num=reshape(network_data(:,7),m,n);
%laurence_error_mat_scalled=laurence_error_mat./laurence_num(20,:);
laurence_error_mat_scalled=bsxfun(@rdivide, laurence_error_mat, laurence_num(Df,:)); %% scaling by dividing with final numerical value
laurence_error_mat_scalled(isnan(laurence_error_mat_scalled)) = 0;


%optimizer1_error=network_data(:,20);
optimizer1_error=abs(network_data(:,l1)-network_data(:,n1D));
optimizer1_error_mat=reshape(optimizer1_error,m,n);
optimizer1_num=reshape(network_data(:,l1),m,n);


index1=find(optimizer1_error_mat(mid_val,:)>0.6); %%discarding high error which does not change the state with high coupling strength
% index1=find(abs(optimizer1_error_mat)>1);
% optimizer1_error_mat(index1)=10;
optimizer1_error_mat(:,index1)=[];
optimizer1_num(:,index1)=[];
% index3=find(mean(optimizer1_error_mat)>.8);
% optimizer1_error_mat(:,index3)=[];
%optimizer1_error_mat_scalled=optimizer1_error_mat./optimizer1_num(20,:);
optimizer1_error_mat_scalled=bsxfun(@rdivide, optimizer1_error_mat, optimizer1_num(Df,:));
optimizer1_error_mat_scalled(isnan(optimizer1_error_mat_scalled)) = 0;



optimizer2_error=abs(network_data(:,n2)-network_data(:,n2D));
optimizer2_error_mat=reshape(optimizer2_error,m,n);
optimizer2_num=reshape(network_data(:,13),m,n);

index2=find(optimizer2_error_mat(mid_val,:)>0.6);
optimizer2_error_mat(:,index2)=[];
optimizer2_num(:,index2)=[];
%optimizer2_error_mat_scalled=optimizer2_error_mat./optimizer2_num(20,:);
optimizer2_error_mat_scalled=bsxfun(@rdivide, optimizer2_error_mat, optimizer2_num(Df,:));
optimizer2_error_mat_scalled(isnan(optimizer2_error_mat_scalled)) = 0;

l_e1m=mean(laurence_error_mat_scalled');
l_e1s=std(laurence_error_mat_scalled');
l_e2m=mean(optimizer1_error_mat_scalled');
l_e2s=std(optimizer1_error_mat_scalled');
l_e3m=mean(optimizer2_error_mat_scalled');
l_e3s=std(optimizer2_error_mat_scalled');
errorbar(k(1:end),l_e1m(1:end),l_e1s(1:end))
hold on
errorbar(k(1:end),l_e2m(1:end),l_e2s(1:end))
errorbar(k(1:end),l_e3m(1:end),l_e3s(1:end))