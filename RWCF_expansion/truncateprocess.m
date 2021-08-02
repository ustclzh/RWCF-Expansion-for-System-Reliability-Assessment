function y=truncateprocess(Data,model)%model stands for the components lifetime model 1 for exponential, 2 for weibull(extremevalue), 3 for lognormal(normal).
m=length(Data(:,1));
n=length(Data(1,:));
y=[];
alpha=0.2;
c=n*alpha;
datacomplete=Data(:,1:(n-c));
datacensored=Data(:,(n-c+1):n);
datacensored=mycensored(datacensored);
for i=1:m
    switch model
        case 1
            temp=filling_Expoential(datacomplete(i,:),datacensored(i,:)); 
        case 2
            temp=filling_Extremevalue(datacomplete(i,:),datacensored(i,:));
        case 3
            temp=filling_Normal(datacomplete(i,:),datacensored(i,:));
    end
y=[y;temp];
end
end


function y=mycensored(datacensored)%ÕâÀïÎÒÃÇ½«Êý¾Ý½øÐÐÉ¾Ê§´¦Àí£»
L=min(min(datacensored))-10;
n=length(datacensored(1,:));
m=length(datacensored(:,1));
C=(rand(m,n)+1)/2;
y=(datacensored-L).*C+L;
end



function y=filling_Extremevalue(data_complete,data_censored)%½øÀ´µÄÊÇ¼«Öµ·Ö²¼£¬Ò²¾ÍÊÇÐèÒªÏÈ¶ÔÍþ²¼¶û·Ö²¼È¥¶ÔÊý
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




function y=filling_Normal(data_complete,data_censored)%½øÀ´µÄÊÇÕýÌ¬Êý¾Ý
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




function y=filling_Exponential(data_complete,data_censored)%½øÀ´µÄÊÇÕýÌ¬Êý¾Ý
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


