

cd('C:\CreatisMain\DATA\PhaseNano120');

nameMain    = 'multips_166_2012_femR_1L_120nm_tomo1_';
pathSave    = 'F:\PhaseNano120';

if 1
    numProj     = 2000;
%     ind         = 0;
    indSelect   = 1:124;
    projAll     = zeros(2048,length(indSelect),2000);
    figure;
    parfor ip = 1:numProj
%         ind       = ind + 1;
        filename  = [nameMain num2str(ip,'%04.0f')];
        %   im        = edfread([filename '.edf']);
        data=load(fullfile(pathSave,filename),'im');
        projAll(:,:,ip)  = squeeze(data.im(:,indSelect)); %data.im(:,indSelect);
%         imagesc(projAll(:,:,ip)); colormap gray; colorbar; axis image; title(num2str(ip)); drawnow;
    end    
    
else
    load('multips_166_2012_femR_1L_120nm_tomo_1_10_2000');
end


figure; imagesc(squeeze(projAll(:,1,:))); axis image; colormap gray; colorbar; drawnow;
caxis([-0.06 0.06]);
    
% ANGLE_BETWEEN_PROJECTIONS = 0.090000
N           = size(projAll);
numProj     = 2000;
theta       = (1:numProj)*180/numProj;
rotAxis     = 1043.5;

% FBP
output_size         = round(2*rotAxis)%N(1);
frequency_scaling   = 0.7;

projThis            = zeros([round(2*rotAxis) numProj]);
projThis(1:N(1),:)= projAll(1:N(1),1,:);
imThis      = iradon(squeeze(projThis),theta,'linear','none',output_size);
figure; imagesc(imThis); colormap gray; colorbar; 
caxis([-0.001 0.0017]); drawnow;

projSum     = squeeze(sum(projAll,2));
imThis      = iradon(projSum,theta,'linear','none',output_size);
figure; imagesc(imThis); colormap gray; colorbar; 

figure;
for ip = 1:size(projAll,2)
    projThis    = squeeze(projAll(:,ip,:));
    imThis      = iradon(squeeze(projAll(:,ip,:)),theta,'linear','none',output_size);
    imAll(:,:,ip) = imThis;
    imagesc(imThis); colormap gray; colorbar; colormap gray; colorbar; axis image; title(num2str(ip)); drawnow;
end

imAllMax = max(imAll,[],3);
figure; imagesc(imAllMax); colormap gray; colorbar; 

% % Loop across slices in z
% figure;
% for iz = 1:N(1)
%     sinoThis    = squeeze(projAll(iz,:,:));
%     imagesc(sinoThis); colormap gray; colorbar; drawnow;
%     imThis      = iradon(sinoThis,theta,'linear',frequency_scaling,output_size);
%     imagesc(imThis); colormap gray; colorbar; drawnow;
%     imAll(:,:,iz) = imThis;
% end

