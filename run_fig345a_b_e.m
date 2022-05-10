
 global K A
%%K is the coupling strength
% A is the network adjacency matrix
t0 = 0; % the start time 
tf = 200; % the end time 
x_low = .01;%0.0001;
load E109.txt -ascii %%edhe list of the network
load a109.txt -ascii % eigen vector of the network

  edgelist=E109;

 A=edgeL2adj(edgelist); %converting the edgelist to a adjacency matrix

 %%%reading other data like eigen values
    metad = sprintf('metadata109.txt');
    fileID = fopen(metad);

format long
cell_data= textscan(fileID,'%s%f','Delimiter','=','headerLines',0);
parameters=cell_data{2};
beta_maxev=parameters(1);
betaprxN1=parameters(2);
e_val_max=parameters(7);
alphaprxN1=parameters(8);
alphaprxN2=parameters(9);

%ev=eig_vecs_corrected_netscience_lcc;
   ev=a109;


n = length(A); % number of nodes in the network
x0l = ones(n,1)*x_low;  % Initial conditions
options = [];

deg1=sum(A);

betaprx=beta_maxev;

alphaprx=e_val_max;

network_data=[];


for K=0:.01:1.5
    tic
    [t,x] = ode45(@double_well_work2,[t0,tf],x0l,options); % change the function for different dynamical system
    xl_ss = x(end,:);
 
R_low_prx=sum(ev(:,1).*xl_ss'); %R from spectral method
R_low_naoki2=sum(ev(:,3).*xl_ss'); %R from optimal eigen vector
[t,xll1]=ode45(@double_well1DPRX,[0:.001:100],x_low,options,betaprx,alphaprx); %% approximation of spectral method
[t,xlN1]=ode45(@double_well1DPRX,[0:.001:100],x_low,options,1,alphaprx); %% approximation of spectral method with beta 1
[t,xlN2]=ode45(@double_well1DPRX,[0:.001:100],x_low,options,1,alphaprxN2); %%approximation with optimal eigenvector/eigenvalue
%data_laurence=[data_laurence;D R_low_prx R_hi_prx xll1(end) xhl1(end) R_low_naoki1 R_hi_naoki1 R_low_naoki2 R_hi_naoki2];
network_data=[network_data;K R_low_prx xll1(end) xlN1(end) R_low_naoki2 xlN2(end) ];
toc
save network_data_Double_well.mat data_laurence
end

figure(1);
%%plotting the bifurcation diagram corresponding to spectral method and
%%spectral method with beta=1
subplot(1,3,1)
plot(network_data(:,1),network_data(:,2))
hold on
plot(network_data(:,1),network_data(:,3))
plot(network_data(:,1),network_data(:,4))
%%%%plotting the bifurcation diagram corresponding to optimal eigenvector
%%%%and eigenvalue
subplot(1,3,2)
plot(network_data(:,1),network_data(:,5))
hold on
plot(network_data(:,1),network_data(:,6))

subplot(1,3,3)
plot(data_laurence(:,1),abs(data_laurence(:,2)-data_laurence(:,3))/data_laurence(end,2))
hold on
plot(data_laurence(:,1),abs(data_laurence(:,2)-data_laurence(:,4))/data_laurence(end,2))
plot(data_laurence(:,1),abs(data_laurence(:,5)-data_laurence(:,6))/data_laurence(end,5))