
cd('/home/goubet/Documents/M2_Internship/data/FilesToReadData');

nameMain    = '../multips_166_2012_femR_1L_120nm_tomo1_/multips_166_2012_femR_1L_120nm_tomo1_';
pathSave    = '/home/goubet/Documents/M2_Internship/data//multips_166_2012_femR_1L_120nm_tomo1_';
nameMainSave    = 'multips_166_2012_femR_1L_120nm_stack';



ip = 800;
filename  = [nameMain num2str(ip,'%04.0f')];
%   im        = edfread([filename '.edf']);
data=load(fullfile(pathSave,filename),'im');
proj1 = data.im(:,:);
imagesc(proj1(:,:)); colormap gray; colorbar; axis image; title(num2str(ip)); drawnow;
%caxis([-0.04 0.04]); drawnow;

if 1
    numProj     = 2000;
%     ind         = 0;
    slice = 1000;
    indSelect   = 1:100;
    projAll     = zeros(2048,2000,100);
    %figure;
    parfor ip = 1:numProj
%       ind       = ind + 1;
        filename  = [nameMain num2str(ip,'%04.0f')];
        %   im        = edfread([filename '.edf']);
        data=load(fullfile(pathSave,filename),'im');
        projAll(:,ip,indSelect)  = data.im(:,indSelect); %im(:,indSelect); %  why squeeze?
        %imagesc(projAll(:,:)); colormap gray; colorbar; axis image; title(num2str(ip)); drawnow;
    end    
    imagesc(projAll(:,:,100)); colormap gray; colorbar; axis image; title(num2str(ip)); drawnow;
    caxis([-37 13]); drawnow;
    save(fullfile(pathSave,[nameMainSave '_sinoOneSlice1000_100slices']),'projAll','-v7.3');
else
    %load slice 1000
    load(fullfile(pathSave,[nameMainSave '_sinoOneSlice1000_100slices']));
end


%%%% plot of sinogram 1  25  50  100
figure('Name', 'Sinograms');
subplot(2,2,1);
imagesc(squeeze(projAll(:,:,selected4(1)))); axis image; colormap gray; colorbar; drawnow;
caxis([-25 5]); drawnow;
title('Slice 1000');

subplot(2,2,2);
imagesc(squeeze(projAll(:,:,selected4(2)))); axis image; colormap gray; colorbar; drawnow;
caxis([-25 5]); drawnow;
title('Slice 1025');

subplot(2,2,3);
imagesc(squeeze(projAll(:,:,selected4(3)))); axis image; colormap gray; colorbar; drawnow;
caxis([-25 5]); drawnow;
title('Slice 1050');

subplot(2,2,4);
imagesc(squeeze(projAll(:,:,selected4(4)))); axis image; colormap gray; colorbar; drawnow;
caxis([-25 5]); drawnow;
title('Slice 1100');





%%%%% without zero padding
selected5 = [1:30];%;25;50;100];

N           = size(projAll);
N(2)/2;
frequency_scaling   = 0.7;
output_size         = N(1);


%2000 projection
numProj     = 2000;
theta       = (1:numProj)*180/numProj;
projThisAll            = zeros([length(selected5) output_size numProj]);
for i= 1:length(selected5)
    %projThis2000            = zeros([output_size numProj]);
    projThisAll(i,1:N(2),:)= projAll(1:N(2),:, selected5(i));
    %before1=clock;
    imThis2000      = iradon(squeeze(projThisAll(i,:,:)),theta);%,'linear','none',output_size);
    %timeRad2000= etime(clock,before1);
    figure; imagesc(imThis2000); colormap gray; colorbar; 
    caxis([-0.010 0.001]); drawnow;
    title(['Slice ' num2str(1000 + selected5(i),'%04.0f')]); %selected4(i)]);
end


% with zero padding

selected5 = [1:30];%;25;50;100];

N           = size(projAll);
rotAxis     = 1043.5;
N(2)/2;

% FBP
output_size         = round(2*rotAxis)%N(1);
frequency_scaling   = 0.7;
numProj     = 2000;
theta       = (1:numProj)*180/numProj;
projThisAll            = zeros([length(selected5) output_size numProj]);
imThis2000 = zeros(length(selected5), 2048, 2048);
for i= 1:length(selected5)
    projThisAll(i,1:N(2),:)= projAll(1:N(2),:, selected5(i));
    %before1=clock;
    imThis2000(i,:,:)     = iradon(squeeze(projThisAll(i,:,:)),theta,2048);
    figure; imagesc(squeeze(imThis2000(i,:,:))); colormap gray; colorbar; 
    caxis([-0.010 0.001]); title(['Slice ' num2str(1000 + selected5(i),'%04.0f')]); %selected4(i)]);
    drawnow;
end

imThis25 = imThis2000(1:25,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%
%troncating images
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 3D %%%%%%
targets1 = zeros(length(selected5), 156, 156);
for i= 1:length(selected5)
    targets1(i,:,:) = imThis2000(i,1410:(1410+155), 1500:1655);
    figure; imagesc(squeeze(targets1(i,:,:))); colormap gray; colorbar; 
        caxis([-0.010 0.001]); 
        title(['Slice ' num2str(1000 + selected5(i),'%04.0f')]); 
        drawnow;
end

save(fullfile(pathSave,[nameMainSave '_target_3D_img']),'targets1','-v7.3');


targets1 = zeros(25, 156, 156);
for i= 1:length(selected5-5)
    targets1(i,:,:) = imThis25(i,1410:(1410+155), 1500:1655);
    figure; imagesc(squeeze(targets1(i,:,:))); colormap gray; colorbar; 
        caxis([-0.010 0.001]); 
        title(['Slice ' num2str(1000 + selected5(i),'%04.0f')]); 
        drawnow;
end
        
save(fullfile(pathSave,[nameMainSave '_target_3D_img_25']),'targets1','-v7.3');


targets2 = zeros(length(selected5), 156, 156);
for i= 1:length(selected5)
    targets2(i,:,:) = imThis2000(i,1200:(1200+155), 900:(900+155));
    figure; imagesc(squeeze(targets2(i,:,:))); colormap gray; colorbar; 
        caxis([-0.010 0.001]); 
        title(['Slice ' num2str(1000 + selected5(i),'%04.0f')]); 
        drawnow;
end
  
save(fullfile(pathSave,[nameMainSave '_target2_3D_img']),'targets2','-v7.3');


targets3 = zeros(length(selected5), 156, 156);
for i= 1:length(selected5)
    targets3(i,:,:) = imThis2000(i,1300:(1300+155), 325:(325+155));
    figure; imagesc(squeeze(targets3(i,:,:))); colormap gray; colorbar; 
        caxis('auto'); 
        title(['Slice ' num2str(1000 + selected5(i),'%04.0f')]); 
        drawnow;
end

save(fullfile(pathSave,[nameMainSave '_target3_3D_img']),'targets3','-v7.3');


%%% 2D %%%


target = squeeze(imThis2000(1,1410:(1410+155), 1500:(1500+155)));
figure; imagesc(target); colormap gray; colorbar; 
    caxis([-0.010 0.001]); drawnow;
    title(['Slice ' num2str(1000 + selected5(i),'%04.0f')]); 
        
save(fullfile(pathSave,[nameMainSave '_target_img']),'target','-v7.3');

target2 = squeeze(imThis2000(1,1200:(1200+155), 900:(900+155)));
figure; imagesc(target2); colormap gray; colorbar; 
    caxis([-0.010 0.001]); drawnow;
    title(['Slice ' num2str(1000 + selected5(i),'%04.0f')]); 
        
save(fullfile(pathSave,[nameMainSave '_target2_img']),'target2','-v7.3');

target3 = squeeze(imThis2000(1,1300:(1300+155), 325:(325+155)));
figure; imagesc(target3); colormap gray; colorbar; 
    caxis([-0.010 0.001]); drawnow;
    title(['Slice ' num2str(1000 + selected5(i),'%04.0f')]); 
        
save(fullfile(pathSave,[nameMainSave '_target3_img']),'target3','-v7.3');


%%%%%focus on one slice and multiple freq scal

frequency_scaling = [1;.7;.5;.2;.1];

% with zero padding
N           = size(projAll);
rotAxis     = 1043.5;
N(2)/2;

% FBP
output_size         = round(2*rotAxis)%N(1);
numProj     = 2000;
theta       = (1:numProj)*180/numProj;
for i= 4:5
    projThis2000            = zeros([output_size numProj]);
    projThis2000(1:N(2),:)= projAll(1:N(2),:, selected4(1));
    %before1=clock;
    %add freq in the iradon
    imThis2000      = iradon(squeeze(projThis2000),theta, 'linear',frequency_scaling(i));
    %timeRad2000= etime(clock,before1);
    figure; imagesc(imThis2000); colormap gray; colorbar; 
    caxis([-0.0150 .001]); drawnow;
    title(['Frequency scale ' sprintf('%03.0g', frequency_scaling(i))]); %selected4(i)]);

end


%%%%%focus on one slice and reduce change interpolation

selected_interpolation = {'linear'; 'nearest'; 'spline'; 'pchip'; 'v5cubic'};

N           = size(projAll);
rotAxis     = 1043.5;
N(2)/2;

% FBP
output_size         = round(2*rotAxis)%N(1);
frequency_scaling   = 0.7;
numProj     = 2000;
theta       = (1:numProj)*180/numProj;
time = zeros(5);
for i= 1:5
    projThis2000            = zeros([output_size numProj]);
    projThis2000(1:N(2),:)= projAll(1:N(2),:, selected4(1));
    before1=clock;
    imThis2000      = iradon(squeeze(projThis2000),theta, selected_interpolation{i});
    time(i)= etime(clock,before1);
    figure; imagesc(imThis2000); colormap gray; colorbar; 
    caxis([-0.0150 .001]); drawnow; 
    title(selected_interpolation{i});
end





%%%%%focus on one slice and reduce change filter

selected_filter = {'None'; 'Hann'; 'Hamming'; 'Cosine'; 'Shepp-Logan'; 'Ram-Lak'};

N           = size(projAll);
rotAxis     = 1043.5;
N(2)/2;

% FBP
output_size         = round(2*rotAxis)%N(1);
frequency_scaling   = 0.7;
numProj     = 2000;
theta       = (1:numProj)*180/numProj;
time = zeros(5);
for i= 1:6
    projThis2000            = zeros([output_size numProj]);
    projThis2000(1:N(2),:)= projAll(1:N(2),:, selected4(1));
    before1=clock;
    imThis2000      = iradon(squeeze(projThis2000),theta, selected_filter{i});
    time(i)= etime(clock,before1);
    figure; imagesc(imThis2000); colormap gray; colorbar; 
    caxis([-0.0150 .001]); drawnow; 
    title(selected_filter{i});
end




%%%%%focus on 1 slice and reduce number of projections

%2000

numProj     = 2000;
theta       = (1:numProj)*180/numProj;
projThis2000            = zeros([round(2*rotAxis) numProj]);
projThis2000(1:N(2),:)= projAll(1:N(2),:, selected4(1));
before1=clock;
imThis2000      = iradon(squeeze(projThis2000),theta)%,'linear','none',output_size);
timeRad2000= etime(clock,before1);
figure; imagesc(imThis2000); colormap gray; colorbar; 
title('2000 projections');
 caxis([-0.0150 .001]); drawnow; 

save(fullfile(pathSave,[nameMainSave '_img_2000v_Slice1000']),'imThis2000','-mat');
    
%1000

numProj     = 1000;
theta       = (1:numProj)*180/numProj;
projThis1000 = zeros([round(2*rotAxis) numProj]);
projThis1000(1:N(2),:) = projAll(1:N(2),1:2:2000, selected4(1));
before1=clock;
imThis1000      = iradon(squeeze(projThis1000),theta);%,'linear','none',output_size);
timeRad1000= etime(clock,before1);
figure; imagesc(imThis1000); colormap gray; colorbar; 
caxis([-0.0150 0.001]); drawnow;
title('1000 projections');
peaksnr1000 = psnr(imThis2000,imThis1000) 

save(fullfile(pathSave,[nameMainSave '_img_1000v_Slice1000']),'imThis1000','-mat');


%400
numProj     = 400;
theta       = (1:numProj)*180/numProj;
projThis400 = zeros([round(2*rotAxis) numProj]);
projThis400(1:N(2),:) = projAll(1:N(2),1:5:2000, selected4(1));
before1=clock;
imThis400      = iradon(squeeze(projThis400),theta);%,'linear','none',output_size);
timeRad400= etime(clock,before1);
figure; imagesc(imThis400); colormap gray; colorbar; 
caxis([-0.0150 0.001]); drawnow;
title('400 projections');
peaksnr400 = psnr(imThis2000,imThis400) 

save(fullfile(pathSave,[nameMainSave '_img_400v_Slice1000']),'imThis400','-mat');


%200
numProj     = 200;
theta       = (1:numProj)*180/numProj;
projThis200 = zeros([round(2*rotAxis) numProj]);
projThis200(1:N(2),:) = projAll(1:N(2),1:10:2000, selected4(1));
before1=clock;
imThis200      = iradon(squeeze(projThis200),theta);%,'linear','none',output_size);
timeRad200= etime(clock,before1);
figure; imagesc(imThis200); colormap gray; colorbar; 
caxis([-0.0150 0.001]); drawnow;
title('200 projections');
peaksnr200 = psnr(imThis2000,imThis200) 

save(fullfile(pathSave,[nameMainSave '_img_200v_Slice1000']),'imThis200','-mat');

%100
numProj     = 100;
theta       = (1:numProj)*180/numProj;
projThis100 = zeros([round(2*rotAxis) numProj]);
projThis100(1:N(2),:) = projAll(1:N(2),1:20:2000, selected4(1));
before1=clock;
imThis100      = iradon(squeeze(projThis100),theta);%,'linear','none',output_size);
timeRad100= etime(clock,before1);
figure; imagesc(imThis100); colormap gray; colorbar; 
caxis([-0.0150 0.001]); drawnow;
title('100 projections');
peaksnr100 = psnr(imThis2000,imThis100) 

save(fullfile(pathSave,[nameMainSave '_img_100v_Slice1000']),'imThis100','-mat');




projThis100full = zeros([round(2*rotAxis) 2000]);

for i=1:20:2000
    projThis100full(1:N(2),i) = projAll(1:N(2),i,1);
end

before1=clock;
imThis100full      = iradon(squeeze(projThis100),theta);%,'linear','none',output_size);
timeRad100full= etime(clock,before1);
figure; imagesc(imThis100full); colormap gray; colorbar; 
caxis([-0.0150 0.001]); drawnow;
title('100 projections full');

diffTime = timeRad100full - timeRad100;

errorEval100 = norm(imThis2000 - imThis100)/norm(imThis2000);
errorEval100full = norm(imThis2000 - imThis100full)/norm(imThis2000);

diffErrors = errorEval100full - errorEval100;


%%%% different filters 100 proj

selected_filter = {'None'; 'Hann'; 'Hamming'; 'Cosine'; 'Shepp-Logan'; 'Ram-Lak'};

N           = size(projAll);
rotAxis     = 1043.5;
N(2)/2;

% FBP

numProj     = 100;
theta       = (1:numProj)*180/numProj;
projThis100 = zeros([round(2*rotAxis) numProj]);
projThis100(1:N(2),:) = projAll(1:N(2),1:20:2000, selected4(1));

output_size         = round(2*rotAxis);%N(1);
frequency_scaling   = 0.7;
time = zeros(5);
for i= 1:6
    before1=clock;
    imThis2000      = iradon(squeeze(projThis100),theta, selected_filter{i});
    time(i)= etime(clock,before1);
    figure; imagesc(imThis2000); colormap gray; colorbar; 
    caxis([-0.0150 .001]); drawnow; 
    title(selected_filter{i});
end


%%%%% different rec on 100 proj

selected_interpolation = {'linear'; 'nearest'; 'spline'; 'pchip'; 'v5cubic'};

numProj     = 100;
theta       = (1:numProj)*180/numProj;
projThis100 = zeros([round(2*rotAxis) numProj]);
projThis100(1:N(2),:) = projAll(1:N(2),1:20:2000, selected4(1));

output_size         = round(2*rotAxis);%N(1);
frequency_scaling   = 0.7;
time = zeros(5);
for i= 1:5
    before1=clock;
    imThis2000      = iradon(squeeze(projThis100),theta, selected_interpolation{i});
    time(i)= etime(clock,before1);
    figure; imagesc(imThis2000); colormap gray; colorbar; 
    caxis([-0.0150 .001]); drawnow; 
    title(selected_interpolation{i});
end



% load projections 
% load(fullfile(pathSave,[nameMainSave '_img_100v_Slice1000']), '-mat');
% load(fullfile(pathSave,[nameMainSave '_img_200v_Slice1000']), '-mat');
% load(fullfile(pathSave,[nameMainSave '_img_400v_Slice1000']), '-mat');
% load(fullfile(pathSave,[nameMainSave '_img_1000v_Slice1000']), '-mat');
% load(fullfile(pathSave,[nameMainSave '_img_2000v_Slice1000']), '-mat');

peaksnr100 = psnr(imThis100,imThis2000);
peaksnr200 = psnr(imThis200,imThis2000);
peaksnr400 = psnr(imThis400,imThis2000);
peaksnr1000 = psnr(imThis1000,imThis2000);
peaksnr2000 = psnr(imThis2000,imThis2000);

ssim100 = ssim(imThis100,imThis2000);
ssim200 = ssim(imThis200,imThis2000);
ssim400 = ssim(imThis400,imThis2000);
ssim1000 = ssim(imThis1000,imThis2000);
ssim2000 = ssim(imThis2000,imThis2000);

ssim100 = immse(imThis100,imThis2000);
ssim200 = immse(imThis200,imThis2000);
ssim400 = immse(imThis400,imThis2000);
ssim1000 = immse(imThis1000,imThis2000);
ssim2000 = immse(imThis2000,imThis2000);


%save(fullfile(pathSave,[nameMainSave '_img_100v_Slice1000']),'imThis100','-mat');
%iterative
%numProj     = 2000;
%theta       = (1:numProj)*180/numProj;
%projThis = zeros([round(2*rotAxis) numProj]);
%projThis(1:N(2),:) = projAll(1:N(2),:);
%for is = 1:N(2)
%  if (mod(is,5) == 0) 
%    projThis(1:N(2),is) = projAll(1:N(2),is);
%  end
%end
%imThis      = iradon(squeeze(projThis),theta);%,'linear','none',output_size);
%figure; imagesc(imThis); colormap gray; colorbar; 
%caxis([-0.010 0.001]); drawnow;


projSum     = squeeze(sum(projAll,2));
imThis      = iradon(projSum,theta,'linear',frequency_scaling,output_size);
figure; imagesc(imThis); colormap gray; colorbar; 

figure;
for ip = 1:30%size(projAll,2)
    projThis    = squeeze(projAll(:,ip,:));
    imThis      = iradon(squeeze(projAll(:,ip,:)),theta,'linear',frequency_scaling,output_size);
    imAll(:,:,ip) = imThis;
    imagesc(imThis); colormap gray; colorbar; colormap gray; colorbar; axis image; title(num2str(ip)); drawnow;
end
% Loop across slices in z
figure;
for iz = 1:N(1)
    sinoThis    = squeeze(projAll(iz,:,:));
    imagesc(sinoThis); colormap gray; colorbar; drawnow;
    imThis      = iradon(sinoThis,theta,'linear',frequency_scaling,output_size);
    imagesc(imThis); colormap gray; colorbar; drawnow;
    imAll(:,:,iz) = imThis;
end

