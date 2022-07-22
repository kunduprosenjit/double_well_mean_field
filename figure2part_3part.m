global c A r1 r2 r3 D
tic

r1=1;r2=2;r3=5;D=0.001; 
%c=0; %% c represents the parameter u
Nn=100; %% number of nodes change it to 300,1000,3000,10000 for generating data for fig 2c,d,e,f
N=Nn;
%load A_BAN100ad4.mat
t0 = 0; % the start time 
tf = 200; % the end time 
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


 for ni=1:100 %% number of network simmulation
     %% generating networks
   tic 
  %% For BA scalefree network uncomment next two lines
      el=preferential_attachment(Nn,2);
      A=edgeL2adj(el);
%%%%%%%%%%%For Holme-Kim network uncomment next three lines
%       network_edge=sprintf('HolmeKim_networks/E%d.txt',ni);
% edgelist=load(network_edge);
% A=edgeL2adj0(edgelist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
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
%%finding eigen vaulues and vectors for Laurence
 %fileID = fopen('textfile.txt', 'a');
% [e_vec,e_val]=eig(A);
% e_val_max=e_val(end,end);
% e_vec_cen=e_vec(:,end)/sum(e_vec(:,end));
% betaprx=e_vec_cen'*(diag(deg))*e_vec_cen/e_val_max/(e_vec_cen'*e_vec_cen);
% alphaprx=e_val_max;
     
for c=0:.01:2 %% for parameter u use loop of c instead
%for c=0:.01:4    %% c represents u in the paper
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
         %data 1-c, 2-beta_eff, 3-x_eff 4-xlg1(end) 5-R_low_prx 6xll1(end) 7-meanfield gao, 8-meanfield laurence
    data=[data;c,beta_eff x_eff xlg1(end) mean(theta_xnew(3:2:end).*theta_xnew(4:2:end))/mean(theta_xnew(3:2:end))];
%     d1=[c beta_eff x_eff xl1(end)];
%     fprintf(fileID, '%f %f %f %f \n', d1);
end
data_network=[data_network;ni*ones(length(data(:,1)),1) data];
%save network10_100holmekim_10000N_test_new_varyingD_6.mat data_network
 toc
 end
 
 %data 2-c, 3-beta_eff, 4-x_eff 5-xlg1(end) 6-R_low_prx 7-xll1(end) 8-meanfield gao, 9-meanfield laurence
 %fclose(fileID);
    %save resilience2work_numerical.mat data
    %load network100BAdata_200N.mat
m=201;n=size(data_network,1)/m; at=101;
c_param=data_network(1:m,2); %%parameter
%c=data_network(1:m,2)
gao_numerical=data_network(:,4);
gao_numerical=reshape(gao_numerical,m,n);
gao_1D=data_network(:,5);
gao_1D=reshape(gao_1D,m,n);
gao_mean_field=data_network(:,6);
gao_mean_field=reshape(gao_mean_field,m,n);
threshold=1.5*ones(length(c_param),1);     %%%choose thereshold 
%threshold=2.5*ones(length(D_param),1); 
intersection_data=[];
for i=1:n
    [x1_gaon,y1_gaon]=intersections(c_param,threshold,c_param,gao_numerical(:,i));
    [x1_gao1D,y1_gao1D]=intersections(c_param,threshold,c_param,gao_1D(:,i));
    [x1_gaoMF,y1_gaoMF]=intersections(c_param,threshold,c_param,gao_mean_field(:,i));
    intersection_data=[intersection_data;x1_gaon,y1_gaon,x1_gao1D,y1_gao1D,x1_gaoMF(1),y1_gaoMF(1)];
end

gao_tran_diff=abs(intersection_data(:,1)-intersection_data(:,3));
gao_MF_diff=abs(intersection_data(:,1)-intersection_data(:,5));




diff_gao_mean=mean(abs(intersection_data(:,1)-intersection_data(:,3)));
diff_gao_std=std(abs(intersection_data(:,1)-intersection_data(:,3)));
diff_mf_mean=mean(abs(intersection_data(:,1)-intersection_data(:,5)));
diff_mf_std=std(abs(intersection_data(:,1)-intersection_data(:,5)));

%%Note these values for the network size and save them in a matrix for
%%error bar calculation 
%store [number of node, diff_gao_mean, diff_gao_std,diff_mf_mean,
%diff_mf_std] as the columns of matrix data
%%the eror bars can be plotted using the command 
%%errorbar(data(:,1),data(:,2),data(:,3))
%%hold on
%%errorbar(data(:,1),data(:,4),data(:,5))

% FOR Figure 3ab
 figure(3)
relative_error_gao=abs(data_network(:,4)-data_network(:,5))./data_network(:,4);
relative_error_mean_field=abs(data_network(:,4)-data_network(:,6))./data_network(:,4);
relative_error_mean_field_matrix=reshape(relative_error_mean_field,m,n);
relative_error_gao_matrix=reshape(relative_error_gao,m,n);

relative_error_mean_field_matrix=relative_error_mean_field_matrix';
relative_error_gao_matrix=relative_error_gao_matrix';
error_gao_mean=mean(relative_error_gao_matrix(:,at:end));
error_mean_field_mean=mean(relative_error_mean_field_matrix(:,at:end));
error_gao_std=std(relative_error_gao_matrix(:,at:end));
error_mean_field_std=std(relative_error_mean_field_matrix(:,at:end));

errorbar(c_param(at:end),error_gao_mean,error_gao_std )
hold on
errorbar(c_param(at:end),error_mean_field_mean,error_mean_field_mean )

%%for 3c collect the mean and std of the relative error at c_param=1.5 and plot the errorbar 
% 
%    
%     