global c A r1 r2 r3 D

r1=1;r2=2;r3=5;D=0.001; 
%c=0;
Nn=100; %% number of nodes
N=Nn;
%load A_BAN100ad4.mat
t0 = 0; % the start time 
tf = 500; % the end time 
x_low = .001;%0.0001;
x0 = ones(Nn,1)*x_low;
x1=0:.05:10;
 y1=-(x1-r1).*(x1-r2).*(x1-r3);
   
TFmin = islocalmin(y1);
TFmax = islocalmax(y1);
%plot(x1,y1,x1(TF),y1(TF),'r*')

loc_minx=x1(TFmin); loc_miny=y1(TFmin);
loc_maxx=x1(TFmax); loc_maxy=y1(TFmax);

options=[];
data_network=[];



     %% generating networks
  %% For BA scalefree network uncomment next two lines
      el=preferential_attachment(Nn,2);
      A=edgeL2adj(el);
%%%%%%%%%%%For Holme-Kim network uncomment next three lines
%       network_edge=sprintf('HolmeKim_networks/E%d.txt',ni);
% edgelist=load(network_edge);
% A=edgeL2adj0(edgelist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%for real networks load the network adjacency matrix

cc=largestcomponent(A);
     A=A(cc,cc);
    N = length(A);
x0 = ones(N,1)*x_low;
deg=sum(A);

data=[];
unique_degree=unique(deg);
for i=1:length(unique_degree)
    [~,indices]=find(deg==unique_degree(i));
    freq(i)=length(indices);
end

     
%for D=0:.001:.2 %%for parameter D use this loop instead
for c=0:.001:1   %%c represents the parameter u in the paper 
     theta_x=[];
    %c=c1(i);
    [t1,x] = ode45(@double_well_work2,[t0,tf],x0,options); 
    xl_ss = x(end,:);
    xl_ss=xl_ss';
    %xlss=xl_ss';
    [x_eff,beta_eff]=betaspace(A,xl_ss);

  [t,xlg1]=ode45(@double_well1D,[0,200],x_low,options,beta_eff);
 %   data=[data;c beta_eff x_eff xlg1(end) R_low_prx xll1(end)];

for theta1=.01:.01:10
         k_o=-(loc_miny+c)/theta1/D;
         state_hi=find(unique_degree>k_o);
         state_low=find(unique_degree<k_o);
         theta2=0;
         xana=[];
         for k=1:length(unique_degree)
             
             [x_k,y_k]=intersections(x1,y1,x1,(-c-D*theta1*k)*ones(1,length(y1)));
             if unique_degree(k)<k_o
                 xl_index=find(x_k<loc_minx);
                
                 xk=x_k(xl_index); %lower state
             elseif unique_degree(k)>=k_o    
                 xh_index=find(x_k>=loc_maxx);
                 
                     xk=x_k(xh_index); %higher state
             end
         
             theta2=theta2+k*xk*freq(k)/N/mean(deg);
          xana=[xana, unique_degree(k), xk];
         end
         theta_x=[theta_x;theta1,theta2,xana];
                 
       
end

 index11=find(abs(theta_x(:,1)-theta_x(:,2))==min(abs(theta_x(:,1)-theta_x(:,2))));
         theta_xnew=theta_x(index11,:);
         
         unique_x=theta_xnew(4:2:end);
         x_anallytic=zeros(N,1);
         
         for i=1:length(unique_degree)  
            x_anallytic(find(deg==unique_degree(i)))=unique_x(i);
         end
         %data 1-c, 2-beta_eff, 3-x_eff 4-xlg1(end) 5-meanfield gao,
         data=[data;c,beta_eff x_eff xlg1(end) mean(theta_xnew(3:2:end).*theta_xnew(4:2:end))/mean(theta_xnew(3:2:end))];
         %%comment the above line and uncomment line below for figure
         %%4a,b,c,d
%    data=[data;D,beta_eff x_eff xlg1(end) mean(theta_xnew(3:2:end).*theta_xnew(4:2:end))/mean(theta_xnew(3:2:end))];
end
plot(data(:,1),data(:,3))
hold on
plot(data(:,1),data(:,4))
plot(data(:,1),data(:,5))





   
    