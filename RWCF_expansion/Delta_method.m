function y=Delta_method(alpha,data,t_0,model,R_structure,c,opt)
%data,原始数据，未经过对数变换
z=-[0.5772,-1.9781,5.4444,-23.5528,117.6895,-712.4501];%矩
z_alpha=norminv(1-alpha,0,1);%分位点
n=length(data(1,:));
switch model(1)
    case 1%Bino
        
        
    case 2%Exp
        
        
    case 3%Lognormal
        data=log(data);%transform to normal
        
        mle_est=zeros(1,2*c);
        Est=zeros(1,c);
        diff_r_sys_over_para_mle=zeros(1,2*c);
        Sigma=zeros(1,2*c);
        
        for cc=1:c
            temp=[mean(data(cc,:)),(n-1)*(var(data(cc,:)))^(1/2)/n];
            mle_est((2*cc-1):2*cc)=temp;
            Est(cc)=normcdf((temp(1)-log(t_0))/temp(2));
            diff_r_sys_over_para_mle((2*cc-1):2*cc)=[normpdf((temp(1)-log(t_0))/temp(2)^0.5)/temp(2)^0.5,-normpdf((temp(1)-log(t_0))/temp(2)^0.5)*((temp(1)-log(t_0))/temp(2)^1.5)/2];
            Sigma((2*cc-1):2*cc)=[temp(2)/n,2*(temp(2))^2*(n-1)/(n^2)];
        end
        
        %mle=mle_est;%MLE
        Sigma=diag(Sigma);%渐近方差
        temp=system_diff(Est,R_structure);
        for sss=1:c
            diff_est([2*sss-1,2*sss])=temp(sss);
        end
        partial=diff_est.*diff_r_sys_over_para_mle;
        temp=(partial*Sigma*(partial'))^0.5;
        Est_in=num2cell(Est);
        y=R_structure(Est_in{:})-z_alpha*temp;
        
    case 4%Weibull
        d1R_weib=opt.d1R_weib;
        %d1R_weib=@(d1,d2,d3) d2*d3^d2/d1^(d2+1)*exp(-(d3/d1)^d2);
        d2R_weib=opt.d2R_weib;
        data=log(data);
        mle_est=zeros(c,2);
        for cc=1:c
            mle_est(cc,:)=mle_for_weibull(data(cc,:));
        end
        Inform=[];
        for sss=1:c
            eta_mle=mle_est(sss,1);
            m_mle=mle_est(sss,2);
            temp=inv([m_mle^2/eta_mle^2,-1/eta_mle-z(1)/eta_mle;-1/eta_mle-z(1)/eta_mle,(pi^2/6+(1+z(1))^2)/(m_mle^2)]);%(1+z(1))^2/(m_mle^2)+pi^2/(6*m_mle^2)
            Inform=blkdiag(Inform,temp);
        end
        diff_r_mle=zeros(1,2*c);
        for sss=1:c
            diff_r_mle((2*sss-1):2*sss)=[d1R_weib(mle_est(sss,1),(mle_est(sss,2)),t_0),d2R_weib(mle_est(sss,1),(mle_est(sss,2)),t_0)];
        end
        Est=(exp(-(t_0./mle_est(:,1)).^mle_est(:,2)))';
        diff_est=[];
        temp=system_diff(Est,R_structure);
        for sss=1:c
            diff_est([2*sss-1,2*sss])=temp(sss);
        end
        diff_r_sys_over_para_mle=diff_r_mle.*diff_est;
        temp=(diff_r_sys_over_para_mle*Inform*(diff_r_sys_over_para_mle'))^0.5;
        Est_in=num2cell(Est);
        y=R_structure(Est_in{:})-z_alpha*temp/n^0.5;
        
end
end