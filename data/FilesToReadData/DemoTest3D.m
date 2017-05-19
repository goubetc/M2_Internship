cd('~/Documents/M2_Internship/data/FilesToReadData');
addpath('~/Documents/M2_Internship/matLabImplem/ART/Efficient-ART-Split-Bregman-Reconstruction');
addpath('~/Documents/M2_Internship/data/FilesToReadData/demoSBtronc');
addpath('~/Documents/M2_Internship/matLabImplem/sift-0.9.19/sift');

nameMain    = '../multips_166_2012_femR_1L_120nm_tomo1_/multips_166_2012_femR_1L_120nm_tomo1_';
pathSave    = '~/Documents/M2_Internship/data//multips_166_2012_femR_1L_120nm_tomo1_';
nameMainSave    = 'multips_166_2012_femR_1L_120nm_stack';

pathResSave= '~/Documents/M2_Internship/data/res/SB_Reconstruction/';
nameResSave= 'Tronc_Images';


I = phantom('Modified Shepp-Logan',64);




I3D = zeros(5,64,64);

for i=1:5
    I3D(i,:,:) = I(:,:);
end

numProj     = ceil(64/3);
thetas      = (1:numProj)*180/numProj;
d           = radon(I,thetas);



data = zeros([5, size(d)]);
for i=1:5
%     epsilon = 5; % corrects for negative vlaues in the pre logged sinogram data that would otherwise results in infinite projection values
%     Noise = noise*exp(-d);
%     N_noise = poissrnd(Noise);
%     d2 = -log(N_noise/noise);
%     idx = isinf(d2);
%     d2(idx) = -log(epsilon/noise);    
%     data(i,:,:) = d2(:,:);
    data(i,:,:) = d+(0.05*max(d(:)))*randn(size(d));
end

nBreg       = 1500;


mu          = 1;
lambda      = 1;
gamma       = 1e-4;
alpha       = 1;
        
[recImg2D, errStruct2D, ctrst2D, exTime2D] = SB_3D_2D_Call_Low_Dose(I3D, numProj, nBreg, 0, gamma, alpha, lambda, mu, data);
[recImg3D, errStruct3D, ctrst3D, exTime3D] = SB_3D_Call_Low_Dose(I3D, numProj, nBreg, 0, gamma, alpha, lambda, mu, data);


%plot errors
figure;
plot(errStruct2D.errAll, 'b');
hold on;
plot(errStruct3D.errAll, 'm');
legend('2D','3D')
hold off;

%plot images
%%% 2D
figure('Name', '2D');
for i= 1:size(recImg2D,1)
    subplot(ceil(size(recImg2D,1)/2),ceil(size(recImg2D,1)/ceil(size(recImg2D,1)/2)),i);
    imagesc(squeeze(recImg2D(i,3,:,:))); axis image; colormap gray; colorbar; drawnow;
    %caxis([-25 5]); drawnow;
end

%%% 3D
figure('Name', '3D');
for i= 1:size(recImg3D,1)
    subplot(ceil(size(recImg3D,1)/2),ceil(size(recImg3D,1)/ceil(size(recImg3D,1)/2)),i);
    imagesc(squeeze(recImg3D(i,3,:,:))); axis image; colormap gray; colorbar; drawnow;
    %caxis([-25 5]); drawnow;
end
