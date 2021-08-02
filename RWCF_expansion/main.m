%%程序示例，调用各种函数来计算WCF下限
clc,clear
% 计算特征量还需要输入数据和时间点， 输出列向量[b,v,w,r]
% 计算WCF展开需要输入B,V,W,R,nablaR,laplaceR,n,alpha
% B 特征量的一阶
% V 特征量的二阶
% W 特征量的三阶
% R系统可靠度点估计
% nablaR 系统可靠度关于部件的导数
% laplaceR 系统可靠度关于部件二阶导数
% n 样本量
% alpha 置信度
 
%%
%参数设置：
alpha=0.1;
nn=[5,10,20,50];
c=3;%部件个数
type=1;%series system，目前不能修改
N=2000;
%仿真参数设置
m=[2,2,2];
eta=[5,5,5];

%%
R_para=sym('Rc',[1,c]);
R_structure=prod(R_para);%系统结构统一修改。
R_structure=matlabFunction(R_structure);

%%
z_alpha=norminv(1-alpha,0,1);%分位点
zz=prepare();
syms d1 d2 d3
R_weib=exp(-(d3/d1)^d2);
d1R_weib=diff(R_weib,d1);
d1R_weib=matlabFunction(d1R_weib);
%d1R_weib=@(d1,d2,d3) d2*d3^d2/d1^(d2+1)*exp(-(d3/d1)^d2);
d2R_weib=diff(R_weib,d2);
d2R_weib=matlabFunction(d2R_weib);
%d2R_weib=@(d1,d2,d3) -log(d3/d1)*(d3/d1)^d2*exp(-(d3/d1)^d2);
z=-[0.5772,-1.9781,5.4444,-23.5528,117.6895,-712.4501];%矩

%%
%simulated data
 for i=1:c
    Data(i,:)=(-log(rand(1,N*nn(4)))).^(1/m(i))*eta(i);%Weibull分布
 end


%%
%计算
data=log(Data);%极值分布
for ss=1:4
    n=0.8*nn(ss)
for j=1:50
  j
   t_0=j/10;
   result=zeros(0,N);
   True=exp(-(t_0./eta).^m);
   R_true=prod(True);  %串联直接相乘
   X(j)=R_true;
   
for i=1:N
    %%
    %WCF
    t=log(t_0);
    
    temp=[];
    for cc=1:c
        datal=data(cc,(n*(i-1)+1):(n*i));
    temp=[temp,characteristic_quantity_Extremevalue(datal,t,zz)];
    end
    R_c=temp(4,:);
    B=temp(1,:);
    V=temp(2,:);
    W=temp(3,:);
    dv=temp(5,:);
    
    nablaR=zeros(1,c);
    for cc=1:c
        temp=R_c;
        temp(cc)=[];
        nablaR(cc)=prod(temp);
    end
    laplaceR=zeros(c,c);
    for cc=1:c
        for ccc=1:c
            temp=R_c;
            if cc==ccc
                laplaceR(cc,ccc)=0;
            else
                temp([cc,ccc])=[];
                laplaceR(cc,ccc)=prod(temp);
            end
        end
    end
    R=prod(R_c);
    result(i)=RWCF(B,V,W,R,dv,nablaR,laplaceR,n,alpha);
   %%
   %delta
   %mle_est=[moment_est(Data(j,1:n));moment_est(Data(j,n+1:2*n));moment_est(Data(j,2*n+1:3*n))];
   mle_est=[];
   for cc=1:c
       datal=data(cc,(n*(i-1)+1):(n*i));
       Datal=exp(datal);
       mle_est=[mle_est;mle_for_weibull(Datal)];
   end
    Inform=[];
    for sss=1:c
        eta_mle=mle_est(sss,1);
        m_mle=mle_est(sss,2);
        temp=inv([m_mle^2/eta_mle^2,-1/eta_mle-z(1)/eta_mle;-1/eta_mle-z(1)/eta_mle,(pi^2/6+(1+z(1))^2)/(m_mle^2)]);%(1+z(1))^2/(m_mle^2)+pi^2/(6*m_mle^2)
        Inform=blkdiag(Inform,temp);
    end
   diff_r_mle=[];
    for sss=1:c
        diff_r_mle=[diff_r_mle,d1R_weib(mle_est(sss,1),(mle_est(sss,2)),t_0),d2R_weib(mle_est(sss,1),(mle_est(sss,2)),t_0)];
    end
    Est=(exp(-(t_0./mle_est(:,1)).^mle_est(:,2)))';

    diff_est=[];
    for sss=1:c
        temp=Est;
        temp(sss)=[];
        diff_est([2*sss-1,2*sss])=prod(temp);
    end
    diff_r_sys_over_para_mle=diff_r_mle.*diff_est; 

    temp=(diff_r_sys_over_para_mle*Inform*(diff_r_sys_over_para_mle'))^0.5;
    result_d(i)=prod(Est)-z_alpha*temp/n^0.5;
end
 Coverage_w(ss,j)=sum(result<R_true)/N;
 Coverage_d(ss,j)=sum(result_d<R_true)/N;

end
end
%%
nn=[5,10,20,50];
for jj=1:4
    subplot(2,2,jj);
    set(gcf,'Position',[50/0.277 50/0.277 200/0.277 100/0.277]);
    plot(X,Coverage_w(jj,:),'--'); hold on;
    plot(X,Coverage_d(jj,:),'-.'); hold on;
    plot((0.2:0.01:1),(1-alpha)*ones(1,81),'-');
    title(['n=',num2str(nn(jj))]);
    xlabel('reliability')
    ylabel('coverage rate')
    xlim([0.75,0.95]);
    ylim([0.5,1]);
end
legend('R-WCF','delta','criterion','location','southeast');




