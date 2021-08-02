function y=filling_Normal(data_complete,data_censored)%进来的是正态数据
censored=unique(data_censored);
freq=zeros(1,length(censored));
for i=1:length(censored)
    freq(i)=sum(data_censored==(censored(i)));
end
miu=mean(data_complete);
sigma=(var(data_complete))^0.5;

for i=1:100
    X=[];
    for j=1:length(censored)
        k=1:freq(j);
       Y= norminv((1-normcdf((censored(j)-miu)/sigma))*k/(freq(j)+1)+normcdf((censored(j)-miu)/sigma))*sigma+miu;
       X=[X,Y];
    end
    miu=mean([X,data_complete]);
    sigma=(var([X,data_complete]))^0.5;
end
y=[data_complete,X];
end
