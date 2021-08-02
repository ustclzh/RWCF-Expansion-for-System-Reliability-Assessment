function y=filling_Exponential(data_complete,data_censored)%进来的是正态数据
censored=unique(data_censored);
freq=zeros(1,length(censored));
for i=1:length(censored)
    freq(i)=sum(data_censored==(censored(i)));
end
lambda=mean(data_complete);
for i=1:100
    X=[];
    for j=1:length(censored)
        k=1:freq(j);
       Y=-lambda*log((1-(k/(1+freq(j))))*exp(-censored(j)/lambda));
       X=[X,Y];
    end
    lambda=mean([X,data_complete]);
end
y=[data_complete,X];
end
