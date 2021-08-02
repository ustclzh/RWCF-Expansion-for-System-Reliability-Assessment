function y= characteristic_quantity_Normal(data,t)%this is actually normal data, not log-normal
n=length(data);
mu=mean(data);
sigma=((n-1)*var(data)/n)^(1/2);
R=normcdf((mu-t)/sigma);
fx=normpdf(norminv(R));
fy=-norminv(R)*normpdf(norminv(R))/2;
fxx=-norminv(R)*normpdf(norminv(R));
fyy=(3-norminv(R)^2)*norminv(R)*normpdf(norminv(R))/4;
fxy=-(norminv(R)^2+1)*normpdf(norminv(R))/2;
fxr=-norminv(R);
fyr=-1/2+(norminv(R))^2/2;
%a11=0;
a20=1/n;
a02=2/n;
%a12=0;
%a21=0;
%a30=0;
a03=8/n^2;
a22=2/n^2;
a40=3/n^2;
a04=6/n^2;
%a13=0;
%a31=0;
b=n*(fxx*a20/2+fyy*a02/2);
v=n*(fx^2*a20+fy^2*a02);
dv=2*n*(fx*fxr*a20+fy*fyr*a02 );
w=n^2*(fy^3*a03+3*fx^2*fxx*a40/2+3*fx^2*fyy*a22/2+3*fy^2*fxx*a22/2+3*fx*fy*fxy*a22+3*fy^2*fyy*a04/2);
y=[b,v,w,R,dv]';







