

cd('C:\CreatisMain\DATA\PhaseNano120');

nameMain    = 'multips_166_2012_femR_1L_120nm_tomo1_';
pathSave    = 'F:\PhaseNano120';

if 0
    numProj     = 2000;
    ind         = 0;
    projAll     = zeros(2048,2048,200);

    figure;
    for ip = 1:10:numProj
        ind       = ind + 1;
        filename  = [nameMain num2str(ip,'%04.0f')];
        %   im        = edfread([filename '.edf']);
        load(fullfile(pathSave,filename),'im');
        projAll(:,:,ind)  = im;
        imagesc(im); colormap gray; colorbar; title(num2str(ip)); drawnow;
        clear im;
    end    
    
    figure; imagesc(projAll(:,:,1)); colormap gray; colorbar; drawnow;
    
else
    load('multips_166_2012_femR_1L_120nm_tomo_1_10_2000');
end

N           = size(projAll);
numProj     = 2000;
theta       = (1:10:numProj)*360/numProj;

% FBP
output_size         = N(1);
frequency_scaling   = 0.7;

% Loop across slices in z
figure;
for iz = 1:N(1)
    sinoThis    = squeeze(projAll(iz,:,:));
    imagesc(sinoThis); colormap gray; colorbar; drawnow;
    imThis      = iradon(sinoThis,theta,'linear',frequency_scaling,output_size);
    imagesc(imThis); colormap gray; colorbar; drawnow;
    imAll(:,:,iz) = imThis;
end

