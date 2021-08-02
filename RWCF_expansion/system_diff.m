function y=system_diff(est,structure)
for i=1:length(est)
   est_1=est;
   est_0=est;
   est_1(i)=1;
   est_0(i)=0;
   est_1=num2cell(est_1);
   est_0=num2cell(est_0);
   diff(i)=structure(est_1{:})-structure(est_0{:});
end
y=diff;
end








