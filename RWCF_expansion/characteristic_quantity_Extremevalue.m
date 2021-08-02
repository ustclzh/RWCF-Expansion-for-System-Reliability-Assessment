function f= characteristic_quantity_Extremevalue(data,t,z)%输入数据应该是极值分布
n=length(data);
gamma=0.57721;
sigma=6^(1/2)*((n-1)*var(data)/n)^(1/2)/pi;
mu=mean(data)+gamma*sigma;
R=exp(-exp((t-mu)/sigma));
r1=-gamma;
r2=pi^2/6+gamma^2;
a1=log(log(R^(-1)));
%global c1 x y a f_x f_y f_xr f_yr f_xy f_xx f_yy
%c1=pi/(6^0.5);
% syms x y a;
% f=exp(-exp(psi(1)+c1*(a-x)/((y-x^2)^(1/2))));
% fx=diff(f,x);
% fy=diff(f,y);
% fxr=diff(fx,a)/(-exp(a)*exp(-exp(a)));
% fyr=diff(fy,a)/(-exp(a)*exp(-exp(a)));
% fxx=diff(fx,x);
% fxy=diff(fx,y);
% fyy=diff(fy,y);
f_x=z{1};
f_y=z{2};
f_xr=z{3};
f_yr=z{4};
f_xy=z{5};
f_xx=z{6};
f_yy=z{7};
% f_x=matlabFunction(fx);
% f_y=matlabFunction(fy);
% f_xr=matlabFunction(fxr);
% f_yr=matlabFunction(fyr);
% f_xy=matlabFunction(fxy);
% f_xx=matlabFunction(fxx);
% f_yy=matlabFunction(fyy);
fx=f_x(a1,r1,r2);
fy=f_y(a1,r1,r2);
fxr=f_xr(a1,r1,r2);
fyr=f_yr(a1,r1,r2);
fxx=f_xx(a1,r1,r2);
fxy=f_xy(a1,r1,r2);
fyy=f_yy(a1,r1,r2);

% fx=(R*log(R))*((r2-r1^2)^(-1/2)-r1*(a-r1)*(r2-r1^2)^(-3/2));
% fy=(R*log(R))*(a-r1)*(r2-r1^2)^(-3/2)/2;
% fxx=R^(-1)*(1+(log(R))^(-1))*fx^2-(R*log(R))*((a-r1)*(r2-r1^2)^(-3/2)+3*r1*(a-r1)*(r2-r1^2)^(-5/2));
% fyy=R^(-1)*(1+(log(R))^(-1))*fy^2-3*fy*(r2-r1^2)^(-1)/4;
% fxy=R^(-1)*(1+(log(R))^(-1))*fx*fy-(R*log(R))*((r2-r1^2)^(-3/2)/2-3*r1*(a-r1)*(r2-r1^2)^(-5/2)/2);
z=-[0.5772,-1.9781,5.4444,-23.5528,117.6895,-712.4501];%矩
a20=(z(2)-z(1)^2)/n;
a02=(z(4)-z(2)^2)/n;
a11=(z(3)-z(1)*z(2))/n;
a30=(z(3)-3*z(1)*z(2)+2*z(1)^3)/n^2;
a03=(z(6)-3*z(2)*z(4)+2*z(2)^3)/n^2;
a12=(z(5)-z(1)*z(4)-2*z(2)*z(3)+2*z(1)*z(2)^2)/n^2;
a21=(z(4)-z(2)^2-2*z(1)*z(3)+2*z(2)*z(1)^2)/n^2;
a22=(z(2)*z(4)+2*z(3)^2-z(2)^3-z(1)^2*z(4)-4*z(1)*z(2)*z(3)+3*z(1)^2*z(2)^2)/n^2;
a31=(3*z(2)*z(3)-3*z(1)*z(2)^2-3*z(1)^2*z(3)+2*z(1)^3*z(2))/n^2;
a13=(3*z(3)*z(4)-3*z(1)*z(2)*z(4)-3*z(2)^2*z(3)+2*z(1)*z(2)^3)/n^2;
a40=(3*z(2)^2-6*z(1)^2*z(2)+3*z(1)^4)/n^2;
a04=(3*z(4)^2-6*z(2)^2*z(4)+3*z(2)^4)/n^2;
b=n*(fxx*a20/2+fyy*a02/2+fxy*a11);
v=n*(fx^2*a20+fy^2*a02+2*fx*fy*a11);
dv=n*(2*fx*fxr*a20+2*fy*fyr*a02+2*fx*fyr*a11+2*fy*fxr*a11);
w=n^2*(fx^3*a30+fy^3*a03+3*fx^2*(fy*a21+fxx*a40/2+fyy*a22/2+fxy*a31/2)+3*fy^2*(fx*a12+fyy*a04/2+fxx*a22/2+fxy*a13/2)+3*fx*fy*(fxx*a31/2+fyy*a13/2+fxy*a22));
f=[b,v,w,R,dv]';

  








