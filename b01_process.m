%% Description
% - Loads pre-processed data
% - Fit Muscle-RRM, Self-RRM, and Extended Tofts model
% - Saves results to file

% Estimated runtime: ~20 seconds

%% Initialize 
clearvars
fclose('all');
addpath(genpath('./mfiles'))

%% Configuration

inDir = './data/TCGA-GBM-Results/preprocessed'; % Input data (from previous script)
outDir = './data/TCGA-GBM-Results/postprocessed'; % Output directory

k=3; % Number of clusters
doOverwrite = true; % Overwrite existing files

%% Main

if ~exist(outDir,'dir')
    mkdir(outDir)
end

matFiles = dir([inDir '/*.mat']);

tic;
for i=1:length(matFiles)
    curFile = matFiles(i).name;
    outFile = fullfile(outDir,curFile);
    if exist(outFile,'file') && ~doOverwrite
        continue
    end
    
    % Display patient ID in console as a way of tracking progress
    curFile
    
    %% Load data
    load(fullfile(inDir,curFile)); 
    % Provides: 'Ct','Cp','Crr','t','maskCt','maskCrr','maskCp' ... and more
    %% Basic pre-processing
    
    % Population AIF from literature
    CpPop = GeorgiouAif(t,t(7));

    % Identify voxels with negligible enhancement.
    % In this case, 'negligible' is when the maximum concentration is below 0.01 mM.
    enhancementMask = max(Ct) > 0.01 & max(Ct)<10;
    % Keep track of how many voxels there were
    numGoodVox = sum(enhancementMask(:));
    numVox = sum(maskCt(:));
    % Remove voxels with negligible enhancement
    Ct = Ct(:,enhancementMask);
    maskCt(maskCt) = enhancementMask;
    %% Fit MuscleRRM, SelfRRM, & ETM
    
    % Muscle-RRM using measured AIF tail (not in manuscript)
    [estMusc] = DoMuscleRRM(@CERRM,Ct,Cp,t,Crr);
    
    % Self-RRM using measured AIF tail (not in manuscript)
    [estSelf] = DoSelfRRM(@CERRM,Ct,Cp,t,k);
    CrrSelf = estSelf.Crr;
    
    % Muscle-RRM using population AIF from literature 
    [estMuscPop] = DoMuscleRRM(@CERRM,Ct,CpPop,t,Crr);
    
    % Self-RRM using population AIF from literature 
    [estSelfPop, estSelfPop_all] = DoSelfRRM(@CERRM,Ct,CpPop,t,k);
    
    % Extended Tofts model 
    estETM = struct;
    [estTmp, fittedTmp] = Tofts_LLSQ(Ct, Cp, t, 1);
    estETM.params.Kt = estTmp(:,1);
    estETM.params.ve = estTmp(:,1)./estTmp(:,2);
    estETM.params.kep = estTmp(:,2);
    estETM.params.vp = estTmp(:,3);
    estETM.rss = sum((Ct-fittedTmp).^2);
    
    % Remove fitted curves. This reduces filesize of output mat.
    estMusc = rmfield(estMusc, 'fittedCt');
    estSelf = rmfield(estSelf, 'fittedCt');
    estMuscPop = rmfield(estMuscPop, 'fittedCt');
    estSelfPop = rmfield(estSelfPop, 'fittedCt');
    estSelfPop_all = rmfield(estSelfPop_all, 'fittedCt');
    %% Save results
    
    save(outFile,...
        'estETM','estMusc','estSelf','estMuscPop','estSelfPop','estSelfPop_all',...
        'Crr','CrrSelf','Cp','t',...
        'maskCt','maskCrr','numVox','numGoodVox');
end
toc
%%

%%