function y=truevalue(c,m,para,t)
y=zeros(1,c);
for i=1:c
    switch m(i)
        case 1
            y(i)=para(i,1);
        case 2
            y(i)=exp(t/para(i,1));
        case 3
            y(i)=normodf((-log(t)+para(i,1))/para(i,2));
        case 4
            y(i)=exp((t/para(i,1))^para(i,2));
    end
    
end