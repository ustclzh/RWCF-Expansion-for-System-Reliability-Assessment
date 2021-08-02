function y=simulationdata(N,n,c,m,para)
y=zeros(c,N*n);
for i=1:c
    switch m(i)
        case 1
            y(i,:)=(rand(1,N*n)<para(i,1));
        case 2
            y(i,:)=(-log(rand(1,N*n)))*para(i,1);
        case 3
            y(i,:)=exp(randn(1,N*n)*para(i,2)+para(i,1));
        case 4
            y(i,:)=(-log(rand(1,N*n))).^(1/para(i,2))*para(i,1);
    end
    
end






