cd('/home/goubet/Documents/data/FilesToReadData');
addpath('/home/goubet/Documents/matLabImplem/ART/Efficient-ART-Split-Bregman-Reconstruction');
addpath('/home/goubet/Documents/data/FilesToReadData/demoSBtronc');
addpath('/home/goubet/Documents/matLabImplem/sift-0.9.19/sift');

nameMain    = '../multips_166_2012_femR_1L_120nm_tomo1_/multips_166_2012_femR_1L_120nm_tomo1_';
pathSave    = '/home/goubet/Documents/data//multips_166_2012_femR_1L_120nm_tomo1_';
nameMainSave    = 'multips_166_2012_femR_1L_120nm_stack';

pathResSave= '/home/goubet/Documents/data/res/SB_Reconstruction/';
nameResSave= 'Tronc_Images';


%Load troncatured image
load(fullfile(pathSave,[nameMainSave '_target_img']));
load(fullfile(pathSave,[nameMainSave '_target2_img']));
load(fullfile(pathSave,[nameMainSave '_target3_img']));


%Put images in dar circle
[X Y]   = meshgrid(1:156,1:156);
X       = X - 156/2;
Y       = Y - 156/2;
ind     = ((X.^2+Y.^2)<(156/2-5)^2);

target  = target.*double(ind);
target2  = target2.*double(ind);
target3  = target3.*double(ind);

% %plot targets
% figure; imagesc(target); colormap gray;colormap(flipud(colormap)); colorbar; 
%     caxis('auto'); 
%     title('Target image'); drawnow; %selected(i)]);
%     
% figure; imagesc(target2); colormap gray;colormap(flipud(colormap)); colorbar; 
%     caxis('auto'); 
%     title('Target image'); drawnow; %selected(i)]);
% 
% figure; imagesc(target3); colormap gray;colormap(flipud(colormap)); colorbar; 
%     caxis('auto'); 
%     title('Target image'); drawnow; %selected(i)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Split Bregman senarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setting parameters :
%   target cell array, nbProj table (full dose, half dose, 1/4 dose, 1/10
%   dose), Number of iterations
ImgSize         = size(target, 1);
targets         = {target}%, target2, target3};
numProj         = [ImgSize(1), ceil(ImgSize(1)/2), ceil(ImgSize(1)/4), ceil(ImgSize(1)/10)];
nBreg           = [20, 100, 200, 1000]; % ,5000, 10000]; 

% Reserving memory space
nbTargets       = length(targets);
nbProj          = length(numProj);
nbNBreg         = length(nBreg);

recImg          = zeros(nbTargets, nbProj,nbNBreg, ImgSize(1), ImgSize(1));
recImgnoisy     = zeros(nbTargets, nbProj,nbNBreg, ImgSize(1), ImgSize(1));

exTime          = zeros(nbTargets, nbProj, nbNBreg);
err             = cell(nbTargets, nbProj, nbNBreg);
errnoisy        = cell(nbTargets, nbProj, nbNBreg);
ctrsts          = cell(nbTargets, nbProj, nbNBreg);
ctrstsnoisy     = cell(nbTargets, nbProj, nbNBreg);

%reconstruction call
for i = 1:nbTargets
    for p = 1:nbProj
        for it = 1:nbNBreg
            [recImg(i,p,it,:,:), errStruct, ctrst, exTime(i,p,it)] = SB_Call_Low_Dose(cell2mat(targets(1,i)), numProj(p), nBreg(it));
            err(i,p,it) = {struct2cell(errStruct)};
            ctrsts(i,p,it) = {ctrst};
            [recImgnoisy(i,p,it,:,:), errStruct, ctrst, exTime(i,p,it)] = SB_Call_Low_Dose(cell2mat(targets(1,i)), numProj(p), nBreg(it), 1e4);
            errnoisy(i,p,it) = {struct2cell(errStruct)};
            ctrstsnoisy(i,p,it) = {ctrst};
        end
    end
end


%save results in matLab File
save([fullfile(pathResSave, nameResSave) '_t1-3_5000-10000it'],'recImg', 'exTime', 'err', 'ctrsts', 'recImgnoisy', 'errnoisy', 'ctrstsnoisy','-mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Results Analysis:

load([fullfile(pathResSave, nameResSave) '_t1_20-1000it'], '-mat');

% plot targets
figure; imagesc(target); colormap gray;colormap(flipud(colormap)); colorbar; 
    caxis('auto'); 
    title('Target image'); drawnow;
    
figure; imagesc(target2); colormap gray;colormap; colorbar; 
    caxis('auto'); 
    title('Target image'); drawnow;
    
figure; imagesc(target3); colormap gray;colormap; colorbar; 
    caxis('auto'); 
    title('Target image'); drawnow;



 
targetIm = 1;
nBit = 1;    %max4
nBproj = 4;

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
    imgRec = printRecImg(targetIm, nBproj, nBit, recImg, ImgSize, numProj, nBreg);

% plot all reconstructed image for 1 target images
    % for a specific iteration number
    for i = 1:length(numProj)
        printRecImg(targetIm, i, nBit, recImg, ImgSize, numProj, nBreg);
    end
    % for a specific projection number
    for i = 1:length(nBreg)
        printRecImg(targetIm, nBproj, i, recImg, ImgSize, numProj, nBreg);
    end

% compare execution time
    % between nb proj
    timeProj = plotTimeProj(targetIm, nBit, exTime, numProj, nBreg);   
    % between nb iter
    timeIter = plotTimeIter(targetIm, nBproj, exTime, numProj, nBreg);

% compare errors
    %between nb proj
    [errorProj, ePSNR, tvs] = plotErrorProj(targetIm, nBit, err, numProj, nBreg, recImg, targets);
    %between nb iterations
    errorIter = plotErrorIter(targetIm, nBproj, err, numProj, nBreg);
    
% evaluate edges
targetIm = 1;
nBit = 4;
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
       imgRec = printRecImg(targetIm, i, nBit, recImg, ImgSize, numProj, nBreg);
       linesSB(i,:) =  imgRec(ceil(ImgSize/2),:);
       %linesSB(i,:) = ctrsts{1, 4, i}(1, :);
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
    imgRec = printRecImg(targetIm, 4, nBit, recImg, ImgSize, numProj, nBreg);
    
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