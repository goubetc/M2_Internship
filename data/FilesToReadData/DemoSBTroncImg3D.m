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
% targets1 = abs(targets1);
% targetsNorm = zeros(size(targets1));
% for i = 1:size(targets1,1);
%     minint = min(min(squeeze(targets1(i,:,:))));
%     maxint = max(max(squeeze(targets1(i,:,:))));
%     targetsNorm(i,:,:)=((squeeze(targets1(i,:,:))-minint).*255)./(maxint-minint);
% end
% imhist(squeeze(targetsNorm(1,:,:)));
% 
% for i = 1:size(targets1,2);
%     minint = min(min(squeeze(targetsNorm(:,i,:))));
%     maxint = max(max(squeeze(targetsNorm(:,i,:))));
%     targetsNorm(:,i,:)=((squeeze(targetsNorm(:,i,:))-minint).*0.05)./(maxint-minint);
% end
% 
% for i = 1:size(targets1,2);
%     minint = min(min(min(targets2)));
%     maxint = max(max(max(targets2)));
%     targetsNorm=((targets2-minint).*255)./(maxint-minint);
% end

%targets1 = targetsNorm;
nbImages = size(targets1,1);

%Put images in dar circle
[X Y]   = meshgrid(1:156,1:156);
X       = X - 156/2;
Y       = Y - 156/2;
ind     = ((X.^2+Y.^2)<(156/2-5)^2);
for i= 1:nbImages
    targets1(i,:,:)  = squeeze(targets1(i,:,:)).*double(ind);
    %targets2(i,:,:)  = squeeze(targets2(i,:,:)).*double(ind);
    %targets3(i,:,:)  = squeeze(targets3(i,:,:)).*double(ind);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Split Bregman senarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setting parameters :
%   target cell array, nbProj table (full dose, half dose, 1/4 dose, 1/10
%   dose), Number of iterations

targets1 = squeeze(targets1(1:5,:,:));
%targets1 = squeeze(targets1(1,:,:));
%targets1 = permute(repmat(targets1,[1 1 3]),[3 1 2]);

ImgSize         = size(targets1, 2);
targets         = {targets1};%, targets2, targets3};
numProj         = [ceil((ImgSize(1)*1.6)/10), ceil((ImgSize(1)*1.6)/7), ceil((ImgSize(1)*1.6)/6), ceil((ImgSize(1)*1.6)/4), ceil((ImgSize(1)*1.6)/2)];
nBreg           = 300;


% Reserving memory space
nbTargets        = length(targets);
nbProj           = length(numProj);
%nbNBreg         = length(nBreg);

%recImg          = zeros(nbTargets, nbProj, ceil(nBreg/500)+3, nbImages, ImgSize(1), ImgSize(1));
recImgnoisy2D    = zeros(nbTargets, nbProj, ceil(nBreg/50)+2, nbImages, ImgSize(1), ImgSize(1));
recImgnoisy3D    = zeros(nbTargets, nbProj, ceil(nBreg/50)+2, nbImages, ImgSize(1), ImgSize(1));
uBests2D         = zeros(nbTargets, nbProj, nbImages, ImgSize(1), ImgSize(1));
uBests3D         = zeros(nbTargets, nbProj, nbImages, ImgSize(1), ImgSize(1));

exTime3D         = zeros(nbTargets, nbProj);
exTime2D         = zeros(nbTargets, nbProj);
errnoisy2D       = cell(nbTargets, nbProj);
errnoisy3D       = cell(nbTargets, nbProj);
%ctrsts          = cell(nbTargets, nbProj);
ctrstsnoisy2D    = cell(nbTargets, nbProj);
ctrstsnoisy3D    = cell(nbTargets, nbProj);

mu          = 5;
lambda      = 5;
gamma       = 0;
alpha       = 0.05; %lambda/100;%0.05;
% noise   = 0.001;
noise   = 0.001;

%reconstruction call
for i = 1:nbTargets
    for p = 1:nbProj
%             [recImg2(i,p,:,:,:,:), errStruct, ctrst, exTime(i,p)] = SB_3D_Call_Low_Dose(cell2mat(targets(1,i)), numProj(p), nBreg);
%             err(i,p) = {struct2cell(errStruct)};
%             ctrsts(i,p) = {ctrst};
            %3D
            [recImgnoisy3D(i,p,:,:,:,:), errStruct, ctrst, exTime3D(i,p), uBests3D(i,p,:,:,:)] = SB_3D_Call_Low_Dose(cell2mat(targets(1,i)), numProj(p), nBreg, noise, gamma, alpha, lambda, mu);
            errnoisy3D(i,p) = {struct2cell(errStruct)};
            ctrstsnoisy3D(i,p) = {ctrst};
            %2D
            [recImgnoisy2D(i,p,:,:,:,:), errStruct, ctrst, exTime2D(i,p), uBests2D(i,p,:,:,:)] = SB_3D_2D_Call_Low_Dose(cell2mat(targets(1,i)), numProj(p), nBreg, noise, gamma, alpha, lambda, mu);
            errnoisy2D(i,p) = {struct2cell(errStruct)};
            ctrstsnoisy2D(i,p) = {ctrst};
    end
end


figure; plot(errnoisy{1,1}{1});
%[recImg(i,p,:,:), errStruct, ctrst, exTime(i,p)] = SB_Call_Low_Dose(cell2mat(targets(1,i)), numProj(p), nBreg(it));

%save results in matLab File
save([fullfile(pathResSave, nameResSave) '_3D_and_2D_res'],'recImgnoisy3D','recImgnoisy2D', 'errnoisy3D','errnoisy2D','exTime3D','exTime2D','uBests3D','uBests2D','-mat');


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
    imgRec = printRecImg3D(targetIm, nBproj, recImg2, ImgSize, numProj, nBreg, 2);

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