function z=prepare()
c1=pi/(6^0.5);
syms x y a;
f=exp(-exp(psi(1)+c1*(a-x)/((y-x^2)^(1/2))));
fx=diff(f,x);
fy=diff(f,y);
fxr=diff(fx,a)/(-exp(a)*exp(-exp(a)));
fyr=diff(fy,a)/(-exp(a)*exp(-exp(a)));
fxx=diff(fx,x);
fxy=diff(fx,y);
fyy=diff(fy,y);
f_x=matlabFunction(fx);
f_y=matlabFunction(fy);
f_xr=matlabFunction(fxr);
f_yr=matlabFunction(fyr);
f_xy=matlabFunction(fxy);
f_xx=matlabFunction(fxx);
f_yy=matlabFunction(fyy);
z={f_x,f_y,f_xr,f_yr,f_xy,f_xx,f_yy};
end









