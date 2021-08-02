%%����ʾ�������ø��ֺ���������WCF����
clc,clear
% ��������������Ҫ�������ݺ�ʱ��㣬 ���������[b,v,w,r]
% ����WCFչ����Ҫ����B,V,W,R,nablaR,laplaceR,n,alpha
% B ��������һ��
% V �������Ķ���
% W ������������
% Rϵͳ�ɿ��ȵ����
% nablaR ϵͳ�ɿ��ȹ��ڲ����ĵ���
% laplaceR ϵͳ�ɿ��ȹ��ڲ������׵���
% n ������
% alpha ���Ŷ�
 
%%
%�������ã�
alpha=0.1;
n=50;%������
N=2000;
%�����������
m=[2,1.4,1.8,2,2.5];
eta=[5,6,7,9,9];
c=5;%��������
model=4*ones(1,c);
%% system structure
R_para=sym('Rc',[1,c]);
%ϵͳ�ṹͳһ�޸ģ������Ǵ���ϵͳ�������������͵�ϵͳֱ������ṹ������
%����R_para(1)*(1-(1-R_para(2))*(1-R_para(3)))��2��3�Ų������������һ������������ϵͳ��
R_struct=prod(R_para);
R_structure=matlabFunction(R_struct);
%%



%%
z_alpha=norminv(1-alpha,0,1);%��λ��
zz=prepare(); %necessary for Weibull Distribution
syms d1 d2 d3
R_weib=exp(-(d3/d1)^d2);
d1R_weib=diff(R_weib,d1);
d1R_weib=matlabFunction(d1R_weib);
%d1R_weib=@(d1,d2,d3) d2*d3^d2/d1^(d2+1)*exp(-(d3/d1)^d2);
d2R_weib=diff(R_weib,d2);
d2R_weib=matlabFunction(d2R_weib);
%d2R_weib=@(d1,d2,d3) -log(d3/d1)*(d3/d1)^d2*exp(-(d3/d1)^d2);
z=-[0.5772,-1.9781,5.4444,-23.5528,117.6895,-712.4501];%��
opt.d1R_weib=d1R_weib;
opt.d2R_weib=d2R_weib;
opt.zz=zz;
opt.z=z;
opt.R_structure=R_structure;
%%



%%
%simulated data
para=[eta',m'];
Data=simulationdata(N,n,c,model,para);
%%








%%
%����
data=log(Data);%��ֵ�ֲ�
for j=1:10
   j
   t_0=j/20;
   result=zeros(0,N);
   result_d=zeros(0,N);
   True=truevalue(c,model,para,t_0);%using function to calculate the true value of components reliability
   True_in = num2cell(True);
   R_true=R_structure(True_in{:});  %����ֱ�����
   X(j)=R_true;
   
for i=1:N
    data_d=data(:,(n*(i-1)+1):(n*i));
   %%
    %RWCF
    result(i)= RWCF_system(data_d, c, model, t_0,  opt,alpha);
   %%
   %delta
   result_d(i)=Delta_method(alpha,data_d,t_0,model,R_structure,c,opt);
end
 Coverage_w(j)=sum(result<R_true)/N;
 Coverage_d(j)=sum(result_d<R_true)/N;
 result=sort(result);
 Quantile_w(j)=result(N*(1-alpha));
 result_d=sort(result_d);
 Quantile_d(j)=result_d(N*(1-alpha));
end

plot(X,Coverage_w,'--'); hold on
plot(X,Coverage_d,'-.');hold on;
%plot(x,coverage_o(ss,:),'g');hold on;
plot(X,(1-alpha)*ones(1,10),'-')
xlim([X(10),X(1)]);
xlabel('reliability')
ylabel('coverage rate')
ylim([0.5,1]) 
legend('R-WCF','delta','criterion');




plot(X,Quantile_w,'--'); hold on
plot(X,Quantile_d,'-.');hold on;
%plot(x,coverage_o(ss,:),'g');hold on;
plot(X,X,'-')
xlim([X(10),X(1)]);
xlabel('reliability')
ylabel('quantile')
ylim([min(X),1]) 
legend('R-WCF','delta','criterion');