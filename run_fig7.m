


global A

deg=sum(A);
t0=0;T=50;y0=.01*ones(1,998);N=10000; y01=0.01;
 h = ( T - t0 )/( N -1); % Calulate and store the step - size
 t = linspace( t0 ,T , N ); % A vector to store the time values .
 y = zeros (numel(y0) , N); % Initialize the Y matrix .
 load E109.txt -ascii
 load a109.txt -ascii

  edgelist=E109;
  
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

ev=a109;

betaprx=beta_maxev;
%beta_prx=e_vec_cen'*(diag(deg1))*e_vec_cen/e_val_max/(e_vec_cen'*e_vec_cen);
alphaprx=e_val_max;
%     A=edgeL2adj(edgelist);
%load netscience-lcc.mat -ascii
%edgelist=netscience_lcc;
 A=edgeL2adj(edgelist);
 nn=length(A);
 data=[];
 D=0.2;
 for K=0:.01:4
 y = zeros (numel(y0) , N);
  y11=zeros(1,N);
   yopt1=zeros(1,N);
yopt2=zeros(1,N);
 
 y(:, 1) = y0(:); 
 
 y11(:,1)=y01;
   yopt1(:,1)=y01;
yopt2(:,1)=y01;


 for i = 1:( N -1)
  
y(:, i +1) = y(:, i) + h * SIS1(t(i), y(:, i),K) +D*sqrt(h)*randn(nn,1); % Update approximation y at t+h
 y11(:, i +1) = y11(:, i) + h * SIS1DPRX(t(i), y11(:, i),betaprx,alphaprx,K);% +D*sqrt(h)*randn(); % Update approximation y at t+h
 yopt1(:, i +1) = yopt1(:, i) + h * SIS1DPRX(t(i), yopt1(:, i),1,alphaprx,K);% +D*sqrt(h)*randn(); % Update approximation y at t+h
yopt2(:, i +1) = yopt2(:, i) + h * SIS1DPRX(t(i), yopt2(:, i),1,alphaprxN2,K);% +D*sqrt(h)*randn(); % Update approxim
 y(:,i+1)=mod(y(:,i+1),1);
 end
 
 xh_ss = y(:,end);
 
%xnn=mean(deg'.*y(:,end))/mean(deg);
R_low_prx=sum(ev(:,1).*xh_ss);

R_low_naoki2=sum(ev(:,3).*xh_ss);

data=[data;K,R_low_prx,y11(end),yopt1(end),R_low_naoki2,yopt2(end)];

y=[];
 
 end
 %save data_SIS_noise1.mat data
 figure(2)
 subplot(1,2,1)
plot(data(:,1),data(:,2))
hold on
plot(data(:,1),data(:,3))
plot(data(:,1),data(:,4))

 subplot(1,2,2)
plot(data(:,1),data(:,5))
hold on
plot(data(:,1),data(:,6))


