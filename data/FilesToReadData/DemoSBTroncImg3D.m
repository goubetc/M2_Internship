cd('~/Documents/M2_Internship/data/FilesToReadData');
addpath('~/Documents/M2_Internship/matLabImplem/ART/Efficient-ART-Split-Bregman-Reconstruction');
addpath('~/Documents/M2_Internship/data/FilesToReadData/demoSBtronc');
addpath('~/Documents/M2_Internship/matLabImplem/sift-0.9.19/sift');

nameMain    = '../multips_166_2012_femR_1L_120nm_tomo1_/multips_166_2012_femR_1L_120nm_tomo1_';
pathSave    = '~/Documents/M2_Internship/data//multips_166_2012_femR_1L_120nm_tomo1_';
nameMainSave    = 'multips_166_2012_femR_1L_120nm_stack';

pathResSave= '~/Documents/M2_Internship/data/res/SB_Reconstruction/';
nameResSave= 'Tronc_Images';


%Load troncatured image
load(fullfile(pathSave,[nameMainSave '_target_3D_img']));
load(fullfile(pathSave,[nameMainSave '_target2_3D_img']));
load(fullfile(pathSave,[nameMainSave '_target3_3D_img']));


%%% Normalize gray scale
targets1 = abs(targets1);
targetsNorm = zeros(size(targets1));
for i = 1:size(targets1,1);
    minint = min(min(squeeze(targets1(i,:,:))));
    maxint = max(max(squeeze(targets1(i,:,:))));
    targetsNorm(i,:,:)=((squeeze(targets1(i,:,:))-minint).*0.05)./(maxint-minint);
end

for i = 1:size(targets1,2);
    minint = min(min(squeeze(targetsNorm(:,i,:))));
    maxint = max(max(squeeze(targetsNorm(:,i,:))));
    targetsNorm(:,i,:)=((squeeze(targetsNorm(:,i,:))-minint).*0.05)./(maxint-minint);
end

for i = 1:size(targets1,2);
    minint = min(min(squeeze(targetsNorm(:,:,i))));
    maxint = max(max(squeeze(targetsNorm(:,:,i))));
    targetsNorm(:,:,i)=((squeeze(targetsNorm(:,:,i))-minint).*0.05)./(maxint-minint);
end

targets1 = targetsNorm;

%Put images in dar circle
[X Y]   = meshgrid(1:156,1:156);
X       = X - 156/2;
Y       = Y - 156/2;
ind     = ((X.^2+Y.^2)<(156/2-5)^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Split Bregman senarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setting parameters :
%   target cell array, nbProj table (full dose, half dose, 1/4 dose, 1/10
%   dose), Number of iterations
nbImages = size(targets1,1);


ImgSize         = size(targets1, 2);
targets         = {targets1};%, targets2, targets3};
numProj         = [ceil(ImgSize(1)*1.6)];%, ceil((ImgSize(1)*1.6)/2), ceil((ImgSize(1)*1.6)/4), ceil((ImgSize(1)*1.6)/10)];
nBreg           = [200]; % ,5000, 10000]; 



% Reserving memory space
nbTargets       = length(targets);
nbProj          = length(numProj);
nbNBreg         = length(nBreg);

recImg2          = zeros(nbTargets, nbProj,nbNBreg, nbImages, ImgSize(1), ImgSize(1));
recImgnoisy2     = zeros(nbTargets, nbProj,nbNBreg, nbImages, ImgSize(1), ImgSize(1));

exTime          = zeros(nbTargets, nbProj, nbNBreg);
err             = cell(nbTargets, nbProj, nbNBreg);
errnoisy        = cell(nbTargets, nbProj, nbNBreg);
ctrsts          = cell(nbTargets, nbProj, nbNBreg);
ctrstsnoisy     = cell(nbTargets, nbProj, nbNBreg);

%reconstruction call
for i = 1:nbTargets
    for p = 1:nbProj
        for it = 1:nbNBreg
            [recImg2(i,p,it,:,:,:), errStruct, ctrst, exTime(i,p,it)] = SB_3D_Call_Low_Dose(cell2mat(targets(1,i)), numProj(p), nBreg(it));
            err(i,p,it) = {struct2cell(errStruct)};
            ctrsts(i,p,it) = {ctrst};
%             [recImgnoisy2(i,p,it,:,:,:), errStruct, ctrst, exTime(i,p,it)] = SB_3D_Call_Low_Dose(cell2mat(targets(1,i)), numProj(p), nBreg(it), 1e4);
%             errnoisy(i,p,it) = {struct2cell(errStruct)};
%             ctrstsnoisy(i,p,it) = {ctrst};
        end
    end
end


%save results in matLab File
save([fullfile(pathResSave, nameResSave) '_3D_targets1_Fuuly_50-1000it'],'recImg2', 'exTime', 'err', 'ctrsts', 'recImgnoisy2', 'errnoisy', 'ctrstsnoisy','-mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Results Analysis:

load([fullfile(pathResSave, nameResSave) '_3D_targets1_Fuuly_50-1000it'], '-mat');

% plot targets
figure; imagesc(squeeze(targets1(2,:,:))); colormap gray;colormap(flipud(colormap)); colorbar; 
    caxis('auto'); 
    title('Target image'); drawnow;
    
figure; imagesc(squeeze(targets2(2,:,:))); colormap gray;colormap(flipud(colormap)); colorbar; 
    caxis('auto'); 
    title('Target image'); drawnow;
    
figure; imagesc(squeeze(targets3(2,:,:))); colormap gray;colormap(flipud(colormap)); colorbar; 
    caxis('auto'); 
    title('Target image'); drawnow;



 
targetIm = 1;
nBit = 4;    %max4
nBproj = 1;

% plot target image
    figure; imagesc(targets{targetIm}); colormap gray;colormap(flipud(colormap)); colorbar; 
    caxis('auto'); 
    title(['Target image ' num2str(targetIm)]); drawnow;

% plot FBP 
    [errorFBP, pNoisesFBP, tv, u] = plotFBP(targets{1,targetIm}, numProj(nBproj),1e4);
    
    errorsFBP = zeros(1, length(numProj));
    pNoisesFBP = zeros(1, length(numProj));
    tvs = zeros(1, length(numProj));

    for i = 1:length(numProj)
        [errorsFBP(i), pNoisesFBP(i), tvs(i), imgFBP] = plotFBP(targets{1,targetIm}, numProj(i));
    end
    
% plot image results
    imgRec = printRecImg3D(targetIm, nBproj, nBit, recImg2, ImgSize, numProj, nBreg, 2);

% plot all reconstructed image for 1 target images
    % for a specific iteration number
    for i = 1:length(numProj)
        printRecImg3D(targetIm, i, nBit, recImg2, ImgSize, numProj, nBreg, 2);
    end
    % for a specific projection number
    for i = 1:length(nBreg)
        printRecImg3D(targetIm, nBproj, i, recImg, ImgSize, numProj, nBreg, 2);
    end

% compare execution time
    % between nb proj
    timeProj = plotTimeProj(targetIm, nBit, exTime, numProj, nBreg);   
    % between nb iter
    timeIter = plotTimeIter(targetIm, nBproj, exTime, numProj, nBreg);

% compare errors
    %between nb proj
    [errorProj, ePSNR, tvs] = plotErrorProj3D(targetIm, nBit, err, numProj, nBreg, recImg2, targets,2);
    %between nb iterations
    errorIter = plotErrorIter(targetIm, nBproj, err, numProj, nBreg);
    
% evaluate edges
targetIm = 1;
nBit = 1;
nBproj = 1;

    imgRec = printRecImg(targetIm, nBproj, nBit, recImg, ImgSize, numProj, nBreg);
    %[errorFBP, pnoise, tv, imgRec] = plotFBP(targets{1,targetIm}, numProj(nBproj));
    %imgRec = targets{1,targetIm};
    BW = edge(imgRec, 'canny');%09435);
    figure; imagesc(BW); colormap gray;colormap(flipud(colormap)); 
        caxis('auto'); 
    
    BWtarget = edge(targets{1, targetIm}, 'canny');%09435);
    figure; imagesc(BWtarget); colormap gray;colormap(flipud(colormap)); 
        caxis('auto');
    
    tvEdges = zeros(1, nbProj);
    for i = 1:nbProj
        imgRec = printRecImg(targetIm, i, nBit, recImg, ImgSize, numProj, nBreg);
        %[errorFBP, pnoise, tv, imgRec] = plotFBP(targets{1,targetIm}, numProj(i));
        BW = edge(imgRec, 'canny');%09435);
        figure; imagesc(BW); colormap gray;colormap(flipud(colormap)); 
        caxis('auto'); 
        tvEdges(i) = totalV(BWtarget, BW);
    end
        
% plot lines
    %SB
    % lines
    linesSB = zeros(nbProj, ImgSize);
    % legend for plot
    legend2 = {};
    for i = 1:nbProj
       legend2{i} = [num2str(numProj(i)) ' projections']; %legend2(i) = numProj(i);
       %imgRec = printRecImg(targetIm, i, nBit, recImg, ImgSize, numProj, nBreg);
       %linesSB(i,:) =  imgRec(ceil(ImgSize/2),:);
       linesSB(i,:) = ctrsts{1, 1, i}(1, :);
    end
    % plot
    figure; plot(linesSB(:,:)'); 
    legend(legend2);
    drawnow;
    %FBP
    linesFBP = zeros(nbProj, ImgSize);
    legend2 = {};
    for i = 1:nbProj
       legend2{i} = [num2str(numProj(i)) ' projections']; %legend2(i) = numProj(i);
       [errorFBP, imgRec] = plotFBP(targets{1,targetIm}, numProj(i));
       linesFBP(i,:) =  imgRec(ceil(ImgSize/2), :);
    end
    % plot
    
    figure; plot(linesFBP(:,:)'); 
    legend(legend2);
    drawnow;
    
    
%%% Compute SIFT Descriptor
    %%% Target image
    t = abs(targets{1, targetIm});
    %%% SB Image    
    imgRec = printRecImg3D(targetIm, 4, nBit, recImg2, ImgSize, numProj, nBreg, 2);
    
    %%% FBP Image
    %[errorFBP, pnoise, tv, imgRec] = plotFBP(targets{1,targetIm}, numProj(1));
    
    %%% Edges
    %    imgRec = edge(imgRec, 'canny');
    %    t = edge(t, 'canny');
    
    %%% Display image
    figure; imagesc(abs(imgRec)); colormap gray;colormap(flipud(colormap)); 
        caxis('auto');
    
    %%% Normalization
    maxt = max(max(t));
    lambda = (1/maxt);
    t = t * lambda;
        
    t2 = abs(imgRec);
    t2 = t2 * (1/(max(max(t2))));
    
    
    %%% SIFT
        % target
    [FRAMESTarget,DESCRTarget]=sift(t);
    figure;
    plotsiftdescriptor(DESCRTarget);
    
        % Rec image
    [FRAMESimRec,DESCRimRec]=sift(abs(t2));
    figure;
    plotsiftdescriptor(DESCRimRec);
    
        % matches
    matches = siftmatch(DESCRTarget, DESCRimRec);
    figure;
    plotmatches(t, t2, FRAMESTarget, FRAMESimRec, matches);
    
    
% % AUC
% AUC = sum((Y(1:end-1)+Y(2:end))/2.*...
%   (X(2:end)-X(1:end-1)));