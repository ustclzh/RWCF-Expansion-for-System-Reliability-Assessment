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
N=2000;
%仿真参数设置
m=[3,3,3];
eta=[8,8,8];
%%

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

%%
%simulated data
 for i=1:c
    Data(i,:)=(-log(rand(1,N*nn(4)))).^(1/m(i))*eta(i);%Weibull分布
 end
%%

%%
%calculation
for ss=1:4
    n=0.8*nn(ss);
data=log(Data);%极值分布
for j=1:20
  j
   t_0=j/2;
   result=zeros(0,N);
   True=exp(-(t_0./eta).^m);
   R_true=1-prod(1-True);  %并联
   X(j)=R_true;
   
for i=1:N
    %%
    %WCF
    t=log(t_0);
    
    temp=[];
    for cc=1:c
    temp=[temp,characteristic_quantity_Extremevalue(data(cc,(n*(i-1)+1):(n*i)),t,zz)];
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
        nablaR(cc)=prod(1-temp);
    end
    laplaceR=zeros(c,c);
    for cc=1:c
        for ccc=1:c
            temp=R_c;
            if cc==ccc
                laplaceR(cc,ccc)=0;
            else
                temp([cc,ccc])=[];
                laplaceR(cc,ccc)=-prod(1-temp);
            end
        end
    end
    R=1-prod(1-R_c);
    result(i)=RWCF(B,V,W,R,dv,nablaR,laplaceR,n,alpha);
   %%
   %delta
   %mle_est=[moment_est(Data(j,1:n));moment_est(Data(j,n+1:2*n));moment_est(Data(j,2*n+1:3*n))];
   mle_est=[];
   for cc=1:c
       mle_est=[mle_est;mle_for_weibull(Data(cc,(n*(i-1)+1):(n*i)))];
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
        diff_est([2*sss-1,2*sss])=-prod(1-temp);
    end
    diff_r_sys_over_para_mle=diff_r_mle.*diff_est; 

    temp=(diff_r_sys_over_para_mle*Inform*(diff_r_sys_over_para_mle'))^0.5;
    result_d(i)=1-prod(1-Est)-z_alpha*temp/n^0.5;
end
 Coverage_w(ss,j)=sum(result<R_true)/N;
 Coverage_d(ss,j)=sum(result_d<R_true)/N;
end
end
%%

for jj=1:4
    for i=1:20
        
        if Coverage_w(jj,i)-1+alpha>0
        Coverage_w_p(jj,i)=abs(Coverage_w(jj,i)-1+alpha);
        Coverage_w_n(jj,i)=NaN;
        end
        if Coverage_w(jj,i)-1+alpha<=0
        Coverage_w_n(jj,i)=abs(Coverage_w(jj,i)-1+alpha);
        Coverage_w_p(jj,i)=NaN;
        end
         if Coverage_d(jj,i)-1+alpha>0
        Coverage_d_p(jj,i)=abs(Coverage_d(jj,i)-1+alpha);
        Coverage_d_n(jj,i)=NaN;
        end
        if Coverage_d(jj,i)-1+alpha<=0
        Coverage_d_n(jj,i)=abs(Coverage_d(jj,i)-1+alpha);
        Coverage_d_p(jj,i)=NaN;
        end
    end
end



%%
%ploting
nn=[5,10,20,50];
for jj=1:4
    subplot(2,2,jj);
    set(gcf,'Position',[50/0.277 50/0.277 200/0.277 100/0.277]);
    plot(X,Coverage_w(jj,:),'-- '); hold on;
    plot(X,Coverage_d(jj,:),'-.');
    plot((0.2:0.01:1),(1-alpha)*ones(1,81),'-');
    %plot((0.2:0.01:1),(1-alpha)*ones(1,81),'-');
    hold off
    
    
    title(['n=',num2str(nn(jj))]);
    xlabel('reliability')
    ylabel('coverage rate')
    xlim([0.75,0.95]);
    ylim([0.5,1]);
end
legend('R-WCF','delta','criterion','location','southeast');




