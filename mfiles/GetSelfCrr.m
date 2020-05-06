function [ Crrs, IDX ] = GetSelfCrr(Ct,N,t)

%%
Ct(Ct(:)<0)=0;

[IDX, Crrs] = hierclust2nmf(Ct,N);

for i = 1:N
    Crrs(:,i) = mean(Ct(:,IDX==i),2);
end

sT = length(Crrs);
if nargin < 3
    aucs = trapz(Crrs(1:round(sT/3),:));
else
    mean_crr = mean(Crrs, 2);
    start_frame = find(mean_crr > 0.1*mean(mean_crr), 1, 'first');
    end_frame = find(t > t(start_frame)+1, 1, 'first');
    aucs = trapz(Crrs(start_frame:end_frame,:));
end
[~,idx]=sort(aucs);
Crrs = Crrs(:,idx);


end

