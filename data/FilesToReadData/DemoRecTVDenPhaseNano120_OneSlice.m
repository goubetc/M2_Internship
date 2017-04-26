% DemoRecTVDenPhaseNano120_OneSlice

cd('C:\CreatisMain\DATA\PhaseNano120');

nameMain    = 'multips_166_2012_femR_1L_120nm_tomo1_';
pathSave    = 'F:\PhaseNano120';
nameMainSave    = 'multips_166_2012_femR_1L_120nm_stack';

numProj     = 2000;
Np          = [2048 2048 numProj];
numPStack   = 32;
indLoop     = 1:numPStack:Np(1);
indStack    = 1;

% ANGLE_BETWEEN_PROJECTIONS = 0.090000
theta       = (1:numProj)*180/numProj;

% Correction
rotAxis     = 1043.5;

% -------------------------------------------------------------------------
% Load data
load(fullfile(pathSave,[nameMainSave '_' num2str(indStack)]),'sinoThis');
sinoThis            = permute(sinoThis,[1 3 2]);
Np                  = size(sinoThis);
size(sinoThis)
figure; imagesc(squeeze(sinoThis(:,:,1))); axis image; colormap gray; colorbar; drawnow;
caxis([-0.06 0.06]);
   
tmp                 = zeros([round(2*rotAxis) numProj Np(3)]);
tmp(1:Np(1),:,:)    = sinoThis(1:Np(1),:,:);
sinoThis            = tmp; clear tmp;
Np                  = size(sinoThis);
% -------------------------------------------------------------------------
% FBP
output_size         = round(2*rotAxis)%N(1);
% frequency_scaling   = 0.7;

imThis      = iradon(squeeze(sinoThis(:,:,1)),theta,'linear','none',output_size);
figure; imagesc(imThis); colormap gray; colorbar; 

if 0
    figure;
    for ip = 1:Np(3)
        sinoSlice   = squeeze(sinoThis(:,:,ip));
        imThis      = iradon(squeeze(sinoSlice),theta,'linear','none',output_size);
        imAll(:,:,ip) = imThis;
        imagesc(imThis); colormap gray; colorbar; colormap gray; colorbar; axis image; title(num2str(ip)); drawnow;
    end
    
    imAllMax = max(imAll,[],3);
    figure; imagesc(imAllMax); colormap gray; colorbar;
end

% -------------------------------------------------------------------------
% Undersampling
R           = zeros(Np(1:2));
indSamp     = 1:50:Np(2); length(indSamp)
R(:,indSamp) = 1;

sinoUnd     = sinoThis(:,:,1).*R;
figure; imagesc(sinoUnd); colormap gray; colorbar; 
caxis([-0.06 0.06]);
figure, spy(R)

imUnd      = iradon(squeeze(sinoUnd(:,:,1)),theta,'linear','none',output_size);
figure; imagesc(imUnd); colormap gray; colorbar; 

indx    = 1000:1800;
indy    = 800:1600;
figure; 
subplot(2,1,1); imagesc(imThis(indx,indy)); colorbar; axis image;
subplot(2,1,2);imagesc(imUnd(indx,indy)); colorbar; axis image;
colormap gray;

% TV impainting
mu      = 2;
lambda  = 1;
nin     = 1;
nout    = 20;
sinoImp = TV_SB_impainting_2D(sinoUnd,mu,lambda,nin,nout,R);
imImp   = iradon(squeeze(sinoImp(:,:,1)),theta,'linear','none',output_size);
% figure; imagesc(imImp); colormap gray; colorbar; 
figure; 
subplot(2,2,1); imagesc(imThis(indx,indy)); colorbar; axis image; 
subplot(2,2,2);imagesc(imUnd(indx,indy)); colorbar; axis image;
subplot(2,2,3);imagesc(imImp(indx,indy)); colorbar; axis image;
colormap gray;


% -------------------------------------------------------------------------
% TV denoising
imNoisy     = imThis+abs(min(imThis(:)));
figure; imagesc(imNoisy); colormap gray; colorbar; 
mu      = 0.3;
lambda  = 2*mu;
nin     = 5;
nout    = 2;
tic,imDen   = TV_SB_denoising_2D(imNoisy,mu,lambda,nin,nout);toc

indx    = 1000:1800;
indy    = 800:1600;
figure; 
subplot(2,1,1); imagesc(imNoisy(indx,indy)); colorbar; axis image;
subplot(2,1,2);imagesc(imDen(indx,indy)); colorbar; axis image;
colormap gray;


% -------------------------------------------------------------------------
% SB REC
N           = [Np(1) Np(1)];
A           = @(x)(radon(x,theta));  
AT          = @(x)(iradon(x,theta,'linear','None',1.0,N(1)));
        
mu          = 1; 
lambda      = 1;
gamma       = 1; % gammaAll(im);
alpha       = 1;
nInner      = 1;
nBreg       = 100;
[u,err] = TV_SB_2D_Gen(A,AT,sinoThis(:,:,1),N,mu, lambda, gamma, alpha, nInner, nBreg);

% -------------------------------------------------------------------------
% ART
d           = sinoThis(:,:,1);
thetaThis   = 0;
A           = @(x)(radon(x,thetaThis));  
u           = zeros(Np(1),Np(1));
relaxParam  = 0.1;
numIterART  = 10;
u           = ARTReconstructionRd(A,d,relaxParam,numIterART,u(:),theta);



% -------------------------------------------------------------------------

%