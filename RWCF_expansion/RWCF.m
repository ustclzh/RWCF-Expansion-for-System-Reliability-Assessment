%RWCF展开的基本框架
function y=RWCF(B,V,W,R,dv,nablaR,laplaceR,n,alpha)%均为行向量
h=(sum(nablaR.*nablaR.*V))^(-1/2);
L1=h*sum(nablaR.*B)+h*sum(nablaR.*V)/2-h^3*((V.*nablaR)*laplaceR*((V.*nablaR)'))-h^3*sum(nablaR.^3.*dv.*V)/2;
L2=h^2*sum(nablaR.*nablaR.*W)+3*h*sum(nablaR.*B)-3*h^2*sum(nablaR.^3.*V.*B)+3*h*sum(nablaR.*V)/2+3*h*((V.*nablaR)*laplaceR*((V.*nablaR)'))-9*h^3*((V.*nablaR)*laplaceR*((V.*nablaR)'))-h^3*sum(nablaR.^3.*dv.*V)/2;
g1=-3*L1/2+L2/6;
g3=h^2*L1/2-h^2*L2/6;
z_alpha=norminv(1-alpha);
y=R-z_alpha/(h*n^(1/2))+(g3*z_alpha^2+g1*h^2)*h^(-3)/n;




















