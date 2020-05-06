%% Description
% - Loads fitting results from b01_process.m
% - Makes figures, and shows CCC values

%% Initialize
clearvars
fclose('all');
addpath('./mfiles')

%% Configuration

inDir = './data/TCGA-GBM-Results/postprocessed'; % Input directory (output from c02_doRRIFT.m)

%% Main

matFiles = dir([inDir '/*.mat']);
curFile = matFiles(1).name;
    
load(fullfile(inDir,curFile));

%% Concatenate all the fitted parameters for tumour

estVals.ETM = estETM.params;
estVals.Musc = estMusc.params;
estVals.Self = estSelf.params;
estVals.MuscPop = estMuscPop.params;
estVals.SelfPop = estSelfPop.params;

%% Parameter maps
fieldsMethod = fieldnames(estVals);
fieldsParams = fieldnames(estVals.(fieldsMethod{1}));
maps = struct;
for i=1:length(fieldsMethod)
    for j=1:length(fieldsParams)
        maps.(fieldsParams{j}).(fieldsMethod{i}) = zeros(size(maskCt));
        maps.(fieldsParams{j}).(fieldsMethod{i})(maskCt) = ...
            estVals.(fieldsMethod{i}).(fieldsParams{j});
    end
end

% Crop the image
fieldsA = fieldnames(maps);
fieldsB = fieldnames(maps.(fieldsA{1}));
for i=1:length(fieldsA)
    for j=1:length(fieldsB)
        maps.(fieldsA{i}).(fieldsB{j}) = AutoCrop(maps.(fieldsA{i}).(fieldsB{j}));
    end
end

slice = 10;
majula = jet(100);
% Ktrans
figure('Position',[50,50,500,700]);
clims = [0 0.2];
subplot(3,3,1)
imagesc(maps.Kt.ETM(:,:,slice));
caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
title('Ktrans-ETM')

subplot(3,3,2)
imagesc(maps.Kt.MuscPop(:,:,slice));
caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
title('Ktrans-MuscRRM')

subplot(3,3,3)
imagesc(maps.Kt.SelfPop(:,:,slice));
caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
title('Ktrans-SelfRRM')

% ve
clims = [0 0.4];
subplot(3,3,4)
imagesc(maps.ve.ETM(:,:,slice));
caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
title('ve-ETM')

subplot(3,3,5)
imagesc(maps.ve.MuscPop(:,:,slice));
caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
title('ve-MuscRRM')

subplot(3,3,6)
imagesc(maps.ve.SelfPop(:,:,slice));
caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
title('ve-SelfRRM')

% vp
clims = [0 0.1];
subplot(3,3,7)
imagesc(maps.vp.ETM(:,:,slice));
caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
title('vp-ETM')

subplot(3,3,8)
imagesc(maps.vp.MuscPop(:,:,slice));
caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
title('vp-MuscRRM')

subplot(3,3,9)
imagesc(maps.vp.SelfPop(:,:,slice));
caxis(clims); colormap(majula); set(gca,'XColor','none','YColor','none');
title('vp-SelfRRM')

%% 2D Histograms - Self-RRM vs Tofts model
overlayInfo = [];

doLog = 0;
chosenFieldA = 'ETM';
chosenFieldB = 'SelfPop';
chosenValsA = estVals.(chosenFieldA);
chosenValsB = estVals.(chosenFieldB);

figure('Position',[100,300,1800,400]);

% 2D hist kt

valsA = chosenValsA.Kt;
valsB = chosenValsB.Kt;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.1) = NaN;
    valsB(valsB > 0.1) = NaN;
end

subplot(1,3,1)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,1) = a1;
overlayInfo(2,1) = a2;
sum(isfinite(h.Data(:)))/2./length(h.Data);
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('Self-RRM')
title('kt')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.05:0.2)
customizeFig(gca);

% 2D hist ve
valsA = chosenValsA.ve;
valsB = chosenValsB.ve;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.5) = NaN;
    valsB(valsB > 0.5) = NaN;
end

subplot(1,3,2)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,2) = a1;
overlayInfo(2,2) = a2;
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('Self-RRM')
title('ve')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.1:0.5)
customizeFig(gca);

% 2D hist vp

valsA = chosenValsA.vp;
valsB = chosenValsB.vp;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.05) = NaN;
    valsB(valsB > 0.05) = NaN;
end

subplot(1,3,3)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,3) = a1;
overlayInfo(2,3) = a2;
sum(isfinite(h.Data(:)))/2./length(h.Data);
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('Self-RRM')
title('vp')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.01:0.05)
customizeFig(gca);

disp('------------')
disp([chosenFieldB ' - ' chosenFieldA]) 
disp('CCC and Pearson correlation coefficient for Ktrans, ve, and vp')
disp(overlayInfo)
disp('')
%% 2D Histograms - Muscle-RRM vs Tofts model
overlayInfo = [];

doLog = 0;
chosenFieldA = 'ETM';
chosenFieldB = 'MuscPop';
chosenValsA = estVals.(chosenFieldA);
chosenValsB = estVals.(chosenFieldB);

figure('Position',[100,300,1800,400]);

% 2D hist kt

valsA = chosenValsA.Kt;
valsB = chosenValsB.Kt;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.1) = NaN;
    valsB(valsB > 0.1) = NaN;
end

subplot(1,3,1)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,1) = a1;
overlayInfo(2,1) = a2;
sum(isfinite(h.Data(:)))/2./length(h.Data);
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('Musc-RRM')
title('kt')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.05:0.2)
customizeFig(gca);

% 2D hist ve
valsA = chosenValsA.ve;
valsB = chosenValsB.ve;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.5) = NaN;
    valsB(valsB > 0.5) = NaN;
end

subplot(1,3,2)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,2) = a1;
overlayInfo(2,2) = a2;
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('Musc-RRM')
title('ve')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.1:0.5)
customizeFig(gca);

% 2D hist vp

valsA = chosenValsA.vp;
valsB = chosenValsB.vp;

valsA(valsA<0) = NaN;
valsB(valsB<0) = NaN;

if doLog
    valsA = log10(valsA);
    valsB = log10(valsB);
    valsA(valsA < -5) = NaN;
    valsB(valsB < -5) = NaN;
else
    valsA(valsA > 0.05) = NaN;
    valsB(valsB > 0.05) = NaN;
end

subplot(1,3,3)
h=histogram2(valsA,valsB,100,'DisplayStyle','tile','ShowEmptyBins','on');
[a1,a2] = CCC(valsA,valsB);
overlayInfo(1,3) = a1;
overlayInfo(2,3) = a2;
sum(isfinite(h.Data(:)))/2./length(h.Data);
imagesc(h.XBinEdges,h.YBinEdges,imgaussfilt(log10(h.Values'+1),0.5))
set(gca,'YDir','normal')
colormap('jet')
caxis([0 3])
hold on;
plot([min(valsA) max(valsA)],[min(valsA) max(valsA)],'w')
xlabel('ExtToftsModel')
ylabel('Musc-RRM')
title('vp')
colorbar
pbaspect([1 1 1])
%set(gca,'YTick',0:0.01:0.05)
customizeFig(gca);

disp('------------')
disp([chosenFieldB ' - ' chosenFieldA]) 
disp('CCC and Pearson correlation coefficient for Ktrans, ve, and vp')
disp(overlayInfo)
disp('')
%% Plot of reference tissue curves
figure('Position',[300,300,500,400])
plot(t, CrrSelf, t, Crr, 'LineWidth', 2)
ylabel('Concentration [mM]')
xlabel('Time [min]')
legend('SelfRef', 'Muscle','location','east')
legend boxoff
customizeFig(gca);
ylim([-0.01, 0.4])
%%
