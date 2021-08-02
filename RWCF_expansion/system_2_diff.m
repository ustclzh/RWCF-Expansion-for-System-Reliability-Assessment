function y=system_2_diff(est,structure)
for i=1:length(est)
    for j=length(est)
        if i==j
           diff(i,j)=0;
        else
           est_11=est;
           est_10=est;
           est_01=est;
           est_00=est;
           est_11(i)=1;
           est_11(j)=1;
           est_01(i)=0;
           est_01(j)=1;
           est_10(i)=1;
           est_10(j)=0;
           est_00(i)=0;
           est_00(j)=0;
           est_11=num2cell(est_11);
           est_10=num2cell(est_10);
           est_01=num2cell(est_01);
           est_00=num2cell(est_00);
           diff(i,j)=structure(est_11{:})-structure(est_10{:})+structure(est_00{:})-structure(est_01{:});
        end
   
    end
end
y=diff;
end

