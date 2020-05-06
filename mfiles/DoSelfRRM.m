function [ chosen, all_fits ] = DoSelfRRM(RRM,Ct,Cp,t,N,fTail)
% N = number of CRRs to generate

%%
if ~exist('fTail', 'var')
    fTail = find(t>3, 1, 'first');
end

if ~exist('N', 'var')
    N = 3;
end
%%
CrrSelf = GetSelfCrr(Ct,N,t);
median_rss = zeros(N,1);

for i = 1:N
    Crr = CrrSelf(:,i);
    [ all_fits(i) ] = DoMuscleRRM(RRM,Ct,Cp,t,Crr,fTail);
    median_rss(i) = iqrMean(all_fits(i).rss);
end

[~, min_idx] = min(median_rss);
chosen = all_fits(min_idx);
