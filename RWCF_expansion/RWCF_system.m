%RWCF展开的基本框架
function y=RWCF_system(data, c, model, t_0,  opt,alpha)%均为行向量
n=length(data(1,:));
zz=opt.zz;
R_structure=opt.R_structure;
t=log(t_0);
temp=zeros(5,c);
for cc=1:c
    switch model(cc)
        case 1

        case 2

        case 3
            temp(:,cc)=characteristic_quantity_Normal(data(cc,:),t);
        case 4
            temp(:,cc)=characteristic_quantity_Extremevalue(data(cc,:),t,zz);
    end
end
R_c=temp(4,:);
B=temp(1,:);
V=temp(2,:);
W=temp(3,:);
dv=temp(5,:);
nablaR=system_diff(R_c,R_structure);
laplaceR=system_2_diff(R_c,R_structure);
R_c_in=num2cell(R_c);
R=R_structure(R_c_in{:});
y=RWCF(B,V,W,R,dv,nablaR,laplaceR,n,alpha);
end

%RWCF展开的基本框架
function y=RWCF(B,V,W,R,dv,nablaR,laplaceR,n,alpha)%均为行向量
h=(sum(nablaR.*nablaR.*V))^(-1/2);
L1=h*sum(nablaR.*B)+h*sum(nablaR.*V)/2-h^3*((V.*nablaR)*laplaceR*((V.*nablaR)'))-h^3*sum(nablaR.^3.*dv.*V)/2;
L2=h^2*sum(nablaR.*nablaR.*W)+3*h*sum(nablaR.*B)-3*h^2*sum(nablaR.^3.*V.*B)+3*h*sum(nablaR.*V)/2+3*h*((V.*nablaR)*laplaceR*((V.*nablaR)'))-9*h^3*((V.*nablaR)*laplaceR*((V.*nablaR)'))-h^3*sum(nablaR.^3.*dv.*V)/2;
g1=-3*L1/2+L2/6;
g3=h^2*L1/2-h^2*L2/6;
z_alpha=norminv(1-alpha);
y=R-z_alpha/(h*n^(1/2))+(g3*z_alpha^2+g1*h^2)*h^(-3)/n;
end
































