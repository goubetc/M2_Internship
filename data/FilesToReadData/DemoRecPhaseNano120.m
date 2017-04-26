

cd('C:\CreatisMain\DATA\PhaseNano120');

nameMain    = 'multips_166_2012_femR_1L_120nm_tomo1_';
pathSave    = 'F:\PhaseNano120';

numProj     = 2000;

if 0
    ind         = 0;
    projAll     = zeros(128,1024,200);
    pixelX      = 2048/2-64+1:2048/2+64;
    pixelY      = 2048/2-512+1:2048/2+512;
    figure;
    for ip = 1:numProj
        ind       = ind + 1;
        filename  = [nameMain num2str(ip,'%04.0f')];
        %   im        = edfread([filename '.edf']);
        load(fullfile(pathSave,filename),'im');
        projAll(:,:,ind)  = im(pixelX,pixelY);
        imagesc(projAll(:,:,ind)); colormap gray; axis image; colorbar; title(num2str(ip)); drawnow;
        clear im;
    end    
    
    figure; imagesc(projAll(:,:,1)); colormap gray; colorbar; drawnow;
    
    save('multips_166_2012_femR_1L_120nm_tomo_1_2000_128x1024','projAll');
else
    load('multips_166_2012_femR_1L_120nm_tomo_1_2000_128x1024','projAll');
end

N           = size(projAll);
theta       = (1:numProj)*360/numProj

% FBP
output_size         = N(1);
frequency_scaling   = 0.7;

% Loop across slices in z
figure;
for iz = 1:N(1)
    sinoThis    = squeeze(projAll(:,iz,:));
    subplot(2,1,1); 
    imagesc(sinoThis); colormap gray; colorbar; axis image; drawnow; 
    imThis      = iradon(sinoThis,theta,'linear',frequency_scaling,output_size);
    title(['z=' num2str(iz)]);
    subplot(2,1,2); 
    imagesc(imThis); colormap gray; colorbar; axis image; drawnow; 
    imAll(:,:,iz) = imThis;
end

imSum   = sum(imAll,3);
figure; imagesc(imThis); colormap gray; colorbar; drawnow; axis image

%