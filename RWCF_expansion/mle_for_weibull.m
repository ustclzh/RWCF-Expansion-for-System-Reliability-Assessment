function y=mle_for_weibull(data)
%global n;
%n=length(data);
e=1;
u_l=0;
u_h=100;
u=u_h;
while e>0.0001
   u=(u_l+u_h)/2;
   x=u;
   f=1/x-sum(log(data).*data.^x)/(sum(data.^x))+mean(log(data));
   if f<0
       u_h=u;
   elseif f>0
       u_l=u;
   elseif f==0
       break;
   end
   e=u_h-u_l;
end
m_mle=u;
eta_mle=(mean(data.^m_mle))^(1/m_mle);
y=[eta_mle,m_mle];
end
