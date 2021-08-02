function y=DELTA_weibull(data,t_0,alpha)
n=length(data)/3;
z_alpha=norminv(1-alpha,0,1);%·ÖÎ»µã
syms d1 d2 d3
R_weib=exp(-(d3/d1)^d2);
d1R_weib=diff(R_weib,d1);
d1R_weib=matlabFunction(d1R_weib);
%d1R_weib=@(d1,d2,d3) d2*d3^d2/d1^(d2+1)*exp(-(d3/d1)^d2);
d2R_weib=diff(R_weib,d2);
d2R_weib=matlabFunction(d2R_weib);
%d2R_weib=@(d1,d2,d3) -log(d3/d1)*(d3/d1)^d2*exp(-(d3/d1)^d2);
z=-[0.5772,-1.9781,5.4444,-23.5528,117.6895,-712.4501];%¾Ø
mle_est=[mle_for_weibull(data(1:n));mle_for_weibull(data(n+1:2*n));mle_for_weibull(data(2*n+1:3*n))];
Inform=[];
for sss=1:3
eta_mle=mle_est(sss,1);
m_mle=mle_est(sss,2);
temp=inv([m_mle^2/eta_mle^2,-1/eta_mle-z(1)/eta_mle;-1/eta_mle-z(1)/eta_mle,(pi^2/6+(1+z(1))^2)/(m_mle^2)]);%(1+z(1))^2/(m_mle^2)+pi^2/(6*m_mle^2)
Inform=blkdiag(Inform,temp);
end
diff_r_mle=[d1R_weib(mle_est(1,1),(mle_est(1,2)),t_0),d2R_weib(mle_est(1,1),(mle_est(1,2)),t_0),d1R_weib(mle_est(2,1),(mle_est(2,2)),t_0),d2R_weib(mle_est(2,1),(mle_est(2,2)),t_0),d1R_weib(mle_est(3,1),(mle_est(3,2)),t_0),d2R_weib(mle_est(3,1),(mle_est(3,2)),t_0)];
Est=(exp(-(t_0./mle_est(:,1)).^mle_est(:,2)))';
diff_r_sys_over_para_mle=diff_r_mle.*[Est(2)*Est(3),Est(2)*Est(3),Est(1)*Est(3),Est(1)*Est(3),Est(1)*Est(2),Est(1)*Est(2)]; 
temp=(diff_r_sys_over_para_mle*Inform*(diff_r_sys_over_para_mle'))^0.5;
y=prod(Est)-z_alpha*temp/n^0.5;
 



