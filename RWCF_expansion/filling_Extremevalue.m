function y=filling_Extremevalue(data_complete,data_censored)%�������Ǽ�ֵ�ֲ���Ҳ������Ҫ�ȶ��������ֲ�ȥ����
censored=unique(data_censored);
freq=zeros(1,length(censored));
for i=1:length(censored)
    freq(i)=sum(data_censored==(censored(i)));
end
sigma=6^0.5*(var(data_complete))^0.5/pi;
miu=mean(data_complete)+0.57721*sigma;
for i=1:100
    X=[];
    for j=1:length(censored)
       k=1:freq(j);
       temp=(1-k/(freq(j)+1))*exp(-exp((censored(j)-miu)/sigma));
       Y= log(-log(temp))*sigma+miu;
       X=[X,Y];
    end
sigma=6^0.5*(var([data_complete,X]))^0.5/pi;
miu=mean([data_complete,X])+0.57721*sigma;
end
y=[data_complete,X];
end