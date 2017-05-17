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
load(fullfile(pathSave,[nameMainSave '_target_3D_img_25']));
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

targets1 = squeeze(targets1(1:2,:,:));
ImgSize         = size(targets1, 2);
targets         = {targets1};%, targets2, targets3};
numProj         = [ceil(ImgSize(1)*1.6)/4];%ceil(ImgSize(1)*1.6)/8, ceil(ImgSize(1)*1.6)/4, ceil(ImgSize(1)*1.6)/2, ceil(ImgSize(1)*1.6)];%ceil(ImgSize(1)*1.6), ceil((ImgSize(1)*1.6)/2), ceil((ImgSize(1)*1.6)/4), ceil((ImgSize(1)*1.6)/10)];
nBreg           = 2500; % ,5000, 10000]; 



% Reserving memory space
nbTargets       = length(targets);
nbProj          = length(numProj);

recImg2          = zeros(nbTargets, nbProj, ceil(nBreg/500)+3, nbImages, ImgSize(1), ImgSize(1));
recImgnoisy2     = zeros(nbTargets, nbProj, ceil(nBreg/500)+3, nbImages, ImgSize(1), ImgSize(1));

exTime          = zeros(nbTargets, nbProj);
err             = cell(nbTargets, nbProj);
errnoisy        = cell(nbTargets, nbProj);
ctrsts          = cell(nbTargets, nbProj);
ctrstsnoisy     = cell(nbTargets, nbProj);

%reconstruction call
for i = 1:nbTargets
    for p = 1:nbProj
%             [recImg2(i,p,:,:,:,:), errStruct, ctrst, exTime(i,p)] = SB_3D_2D_Call_Low_Dose(cell2mat(targets(1,i)), numProj(p), nBreg);
%             err(i,p) = {struct2cell(errStruct)};
%             ctrsts(i,p) = {ctrst};
            [recImgnoisy2(i,p,:,:,:,:), errStruct, ctrst, exTime(i,p)] = SB_3D_2D_Call_Low_Dose(cell2mat(targets(1,i)), numProj(p), nBreg, 0.05);
            errnoisy(i,p) = {struct2cell(errStruct)};
            ctrstsnoisy(i,p) = {ctrst};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%ù
%Tune parameters

targets1 = squeeze(targets1(1:2,:,:));
ImgSize         = size(targets1, 2);
targets         = {targets1};%, targets2, targets3};
numProj         = [ceil(ImgSize(1)*1.6)/4, ceil(ImgSize(1)*1.6)/2];%ceil(ImgSize(1)*1.6), ceil((ImgSize(1)*1.6)/2), ceil((ImgSize(1)*1.6)/4), ceil((ImgSize(1)*1.6)/10)];
nBreg           = 1500; % ,5000, 10000]; 

%test 1: gamma, alpha, lambda, mu
        % ref: gamma = 10-2, alpha = lambda = mu = 1
test = [10e-1, 1, 1, 1;
        10e-2, 1, 1, 1;
        10e-4, 1, 1, 1;
        10e-6, 1, 1, 1;
        10e-2, 1, 0.1, 0.1;
        10e-2, 1, 0.01, 0.01;
        10e-2, 1, 1, 2;
        10e-2, 1, 1, 5;
        10e-2, 0.1, 1, 1;
        10e-2, 10, 1, 1];
    
%test 2 ref: gamma = 10-2, alpha = lambda = mu = 1
test = [10e-2, 1, 1, 1; % ref
        0, 1, 1, 1;         %change gamma
        10e-2, 1, 2, 2;     % change lambda and mu
        10e-2, 1, 5, 5;
        10e-2, 1, 1, 6;     % change mu
        10e-2, 1, 1, 6;
        10e-2, 1, 1, 8;
        10e-2, 1, 1, 20;
        10e-2, .1, 1, 1;    % change alpha
        10e-2, .05, 1, 1;
        10e-2, .01, 1, 1;
        10e-2, .005, 1, 1;
        10e-2, .001, 1, 1;
    ];

% test 3 try two best parameters
test = [0, 0.05, 5, 5;
        0, 0.05, 5, 20];

nbTargets       = length(targets);
nbProj          = length(numProj);

recImgnoisyTest3     = zeros(nbTargets, nbProj, size(test,1), ceil(nBreg/500)+3, 2, ImgSize(1), ImgSize(1));
errnoisyTest3        = cell(nbTargets,size(test,1), nbProj);
ctrstsnoisyTest3          = cell(nbTargets,size(test,1), nbProj);
exTimeTest3          = zeros(nbTargets,size(test,1), nbProj);

%reconstruction call
for i = 1:nbTargets
    for p = 1:nbProj
%             [recImg2(i,p,:,:,:,:), errStruct, ctrst, exTime(i,p)] = SB_3D_2D_Call_Low_Dose(cell2mat(targets(1,i)), numProj(p), nBreg);
%             err(i,p) = {struct2cell(errStruct)};
%             ctrsts(i,p) = {ctrst};
        for numtest = 1:size(test,1)
            gamma =     test(numtest,1);
            alpha =     test(numtest,2);
            lambda =    test(numtest,3);
            mu =        test(numtest,4);
            nInner =    1;
            [recImgnoisyTest3(i,p,numtest,:,:,:,:), errStruct, ctrst, exTimeTest3(i,p)] = SB_3D_2D_Call_Low_Dose(cell2mat(targets(1,i)), numProj(p), nBreg, 1e8, gamma, alpha, lambda, mu);
            errnoisyTest3(i,p,numtest) = {struct2cell(errStruct)};
            ctrstsnoisyTest3(i,p,numtest) = {ctrst};
        end
        
    end
end

%%%%%%%
% plot convergence test 1
figure, plot(errnoisy{1,1,2}{1},'b');
hold on;
plot(errnoisy{1,1,1}{1},'y');
plot(errnoisy{1,1,3}{1},'m');
plot(errnoisy{1,1,4}{1},'c');
legend('\gamma=10e-2','\gamma=10e-1','\gamma=10e-4', '\gamma=10e-6');

figure, 
plot(errnoisy{1,1,2}{1},'y');hold on;
plot(errnoisy{1,1,5}{1},'r');
plot(errnoisy{1,1,6}{1},'g');
legend('\mu=1, \lambda=1','\mu=0.1, \lambda=0.1','\mu=0.01, \lambda=0.01');

figure,
plot(errnoisy{1,1,2}{1},'y');hold on;
plot(errnoisy{1,1,7}{1},'b--');
plot(errnoisy{1,1,8}{1},'g--');
legend('\mu=1','\mu=2','\mu=5');

figure,
plot(errnoisy{1,1,2}{1},'y');hold on;
plot(errnoisy{1,1,9}{1},'m--');
plot(errnoisy{1,1,10}{1},'c--');
hold off;
legend('\alpha=1','\alpha=0.1','\alpha=10');
%%%%%%%%%%%%%%%%%%%%%%%%
%ref: gamma = 10-2, alpha = lambda = mu = 1
%        1  10e-2, 1, 1, 1; % ref
%        2  0, 1, 1, 1;         %change gamma
%        3  10e-2, 1, 2, 2;     % change lambda and mu
%        4  10e-2, 1, 5, 5;
%        5  10e-2, 1, 1, 6;     % change mu
%        6  10e-2, 1, 1, 6;
%        7  10e-2, 1, 1, 8;
%        8  10e-2, 1, 1, 20;
%        9  10e-2, .1, 1, 1;    % change alpha
%        10 10e-2, .05, 1, 1;
%        11 10e-2, .01, 1, 1;
%        12 10e-2, .005, 1, 1;
%        13 10e-2, .001, 1, 1;

%%%%%%%
% plot convergence test 2
figure, plot(errnoisyTest2{1,1,1}{1},'b');
hold on;
plot(errnoisyTest2{1,1,2}{1},'c');
legend('gamma=10e-2','gamma=0');

figure, 
plot(errnoisyTest2{1,1,1}{1},'b');hold on;
plot(errnoisyTest2{1,1,3}{1},'c');
plot(errnoisyTest2{1,1,4}{1},'g');
legend('lambda = mu = 1','lambda = mu = 2','labda = mu = 5');

figure,
plot(errnoisyTest2{1,1,1}{1},'b');hold on;
%plot(errnoisyTest2{1,1,5}{1},'y--');
plot(errnoisyTest2{1,1,6}{1},'c--');
plot(errnoisyTest2{1,1,7}{1},'g--');
plot(errnoisyTest2{1,1,8}{1},'m--');
legend('\mu=1','\mu=6','\mu=8', '\mu = 20');

figure,
plot(errnoisyTest2{1,1,1}{1},'b');hold on;
plot(errnoisyTest2{1,1,9}{1},'c--');
plot(errnoisyTest2{1,1,10}{1},'g--');
plot(errnoisyTest2{1,1,11}{1},'m--');
plot(errnoisyTest2{1,1,12}{1},'y--');
plot(errnoisyTest2{1,1,13}{1},'r--');
hold off;
legend('\alpha=.1','\alpha=0.05','\alpha=.005', '\alpha = .001');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%
% compare visually reconstructions
% 
%%% alpha
im1 = squeeze(recImgnoisyTest2(1,1,9,3,1,:,:));
figure; imagesc(im1); colormap gray;colormap(flipud(colormap)); colorbar; 
    %caxis('auto'); 
    title('alpha = 0.1'); drawnow;
im2 = squeeze(recImgnoisyTest2(1,1,10,3,1,:,:));
figure; imagesc(im2); colormap gray;colormap(flipud(colormap)); colorbar; 
    %caxis('auto'); 
    title('alpha = 0.05'); drawnow;
im3 = squeeze(recImgnoisyTest2(1,1,11,3,1,:,:));
figure; imagesc(im3); colormap gray;colormap(flipud(colormap)); colorbar; 
    %caxis('auto'); 
    title('alpha = 0.01'); drawnow;
im4 = squeeze(recImgnoisyTest2(1,1,12,3,1,:,:));
figure; imagesc(im4); colormap gray;colormap(flipud(colormap)); colorbar; 
    %caxis('auto'); 
    title('alpha = 0.005'); drawnow;    
im5 = squeeze(recImgnoisyTest2(1,1,13,3,1,:,:));
figure; imagesc(im5); colormap gray;colormap(flipud(colormap)); colorbar; 
    %caxis('auto'); 
    title('alpha = 0.001'); drawnow;

%%% mu
im1 = squeeze(recImgnoisyTest2(1,1,5,3,1,:,:));
figure; imagesc(im1); colormap gray;colormap(flipud(colormap)); colorbar; 
    %caxis('auto'); 
    title('\mu = 6'); drawnow;
im3 = squeeze(recImgnoisyTest2(1,1,7,3,1,:,:));
figure; imagesc(im3); colormap gray;colormap(flipud(colormap)); colorbar; 
    %caxis('auto'); 
    title('\mu = 8'); drawnow;
im4 = squeeze(recImgnoisyTest2(1,1,8,3,1,:,:));
figure; imagesc(im4); colormap gray;colormap(flipud(colormap)); colorbar;
    %caxis('auto'); 
    title('\mu = 20'); drawnow;
    
    
%%% Lambda and mu
im1 = squeeze(recImgnoisyTest2(1,1,3,3,1,:,:));
figure; imagesc(im1); colormap gray;colormap(flipud(colormap)); colorbar; 
    %caxis('auto'); 
    title('\mu = \lambda = 6'); drawnow;
im3 = squeeze(recImgnoisyTest2(1,1,7,3,1,:,:));
figure; imagesc(im3); colormap gray;colormap(flipud(colormap)); colorbar; 
    %caxis('auto'); 
    title('\mu = \lambda = 8'); drawnow;
    
    
%find minimum in each test
minimums = zeros(size(test,1),1);
for i = 1:size(test,1)
   minimums(i) = min(errnoisy{1,1,i}{1});
end

% [recImg(i,p,:,:), errStruct, ctrst, exTime(i,p)] = SB_Call_Low_Dose(cell2mat(targets(1,i)), numProj(p), nBreg);

%save results in matLab File
save([fullfile(pathResSave, nameResSave) 'parameters_tunig'], 'exTime','ctrstsnoisyTest2', 'ctrsts', 'recImgnoisy2','recImgnoisyTest2', 'errnoisy', 'errnoisyTest2', 'ctrstsnoisy','-mat');%,'err','recImg2','-mat');

%_3D_2D_targets1_low_dose_to_Fuly_2500it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Results Analysis:

load([fullfile(pathResSave, nameResSave) 'parameters_tunig'], '-mat');

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

recimgnoisy3 = zeros(size(recImgnoisy2,5), size(recImgnoisy2,6), size(recImgnoisy2,4));
for i = 1:size(recImgnoisy2,4);
    recimgnoisy3(:,:,i) = squeeze(recImgnoisy2(1,3,8,i,:,:));
end


recimgnoisy3 = abs(recimgnoisy3*1e4);
fileID = fopen([fullfile(pathResSave, nameResSave) '25_projection.raw'],'w');
fwrite(fileID,recimgnoisy3, 'uint32');
fclose(fileID);

  
targetIm = 1;
nBit = 1;    %max4
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