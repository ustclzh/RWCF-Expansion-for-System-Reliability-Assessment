function y= characteristic_quantity_Exponential(data,t)
n=length(data);
R=exp(t/mean(data));
b=R*log(R)*(2+log(R))/(2*n);
v=(R*log(R))^2/n;
w=2*(R*log(R))^3/(n^2)+9*(R*log(R))^3*(2+log(R))/(2*n^2);
y=[b,v,w,R]';
