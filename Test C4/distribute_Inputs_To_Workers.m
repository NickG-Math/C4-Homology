function [rangeN,rangeM,distributedInputs]= distribute_Inputs_To_Workers(rangeN,rangeM)
total=rangeN*rangeM;
for i=1:12 %Slightly increase the range of n so that we can assign equal load to workers
    if mod(total,12)==0
        break
    else
        rangeN=rangeN+1;
        total=rangeN*rangeM;
    end
end
[n,m]=ind2sub([rangeN,rangeM],1:total);
weight=n+2.^m;
[~,I]=sort(weight);
n=n(I);
m=m(I);
inputs=[n;m];
counter=1;
while ~isempty(inputs)
    distributedInputs(:,counter)=inputs(:,1);
    distributedInputs(:,counter+1)=inputs(:,end);
    inputs(:,1)=[];
    inputs(:,end)=[];
    counter=counter+2;
end
end