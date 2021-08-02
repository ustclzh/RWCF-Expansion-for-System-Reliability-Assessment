function y=Weibull_siimulation(eta,m,n)
data=rand(1,n);
y=(-log(data))^(1/m)*eta;
