function [ est ] = DoMuscleRRM(RRM,Ct,Cp,t,Crr,fTail,estKepRR)
% N = number of CRRs to generate

%%
if ~exist('fTail', 'var') || isempty(fTail)
    fTail = find(t>3, 1, 'first');
end

%%
if exist('estKepRR', 'var')
    [pkRRM, fittedCt] = RRM(Ct,Crr,t, estKepRR);
else
    [pkRRM, fittedCt, estKepRR] = RRM(Ct,Crr,t);
end

estKtRR = RRIFT(Cp(fTail:end), Crr(fTail:end), t(fTail:end), estKepRR);
estVeRR = estKtRR/estKepRR;

estKt = pkRRM(:,1) .* estKtRR;
estVe = pkRRM(:,2) .* estVeRR;
estKep = pkRRM(:,3);
if size(pkRRM,2)>3
    estVp = pkRRM(:,4) .* estKtRR;
else
    estVp = zeros(size(estKt));
end

%%
rss = sum((Ct-fittedCt).^2);

est = struct;

est.params.Kt = estKt;
est.params.kep = estKep;
est.params.ve = estVe;
est.params.vp = estVp;

est.paramsRR.Kt = estKtRR;
est.paramsRR.kep = estKepRR;
est.paramsRR.ve = estVeRR;

est.fittedCt = fittedCt;
est.Crr = Crr;
est.rss = rss;
end

