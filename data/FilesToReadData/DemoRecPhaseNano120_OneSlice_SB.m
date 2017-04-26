cd('/home/goubet/Documents/data/FilesToReadData');
addpath('/home/goubet/Documents/matLabImplem/ART/Efficient-ART-Split-Bregman-Reconstruction');
nameMain    = '../multips_166_2012_femR_1L_120nm_tomo1_/multips_166_2012_femR_1L_120nm_tomo1_';
pathSave    = '/home/goubet/Documents/data//multips_166_2012_femR_1L_120nm_tomo1_';
nameMainSave    = 'multips_166_2012_femR_1L_120nm_stack';



ip = 800;
filename  = [nameMain num2str(ip,'%04.0f')];
%   im        = edfread([filename '.edf']);
data=load(fullfile(pathSave,filename),'im');
proj1 = data.im(:,:);
imagesc(proj1(:,:)); colormap gray; colorbar; axis image; title(num2str(ip)); drawnow;
%caxis([-0.04 0.04]); drawnow;

if 0
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
    load(fullfile(pathSave,[nameMainSave '_small_p']));
    
end



%selected = [1;25;50;100];
selected = [1:5];
ip = 800;

%%%% plot of sinogram 1  25  50  100

figure('Name', 'Sinograms');
for i= 1:length(selected)
    subplot(ceil(length(selected)/2),ceil(length(selected)/ceil(length(selected)/2)),i);
    imagesc(squeeze(projAll(:,:,selected(1)))); axis image; colormap gray; colorbar; drawnow;
    caxis([-25 5]); drawnow;
    title(['Slice ' num2str(1000 + selected(i),'%04.0f')]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FBP with zero padding
N           = size(projAll);
rotAxis     = 1043.5;
N(2)/2;

% FBP
output_size         = round(2*rotAxis)%N(1);
frequency_scaling   = 0.7;
numProj     = 2000;
thetas       = (1:numProj)*180/numProj;
%imThis2000 = zeros(length(selected), output_size,output_size);
for i= 1:length(selected)
    projThis2000            = zeros([output_size numProj]);
    projThis2000(1:N(2),:)= projAll(1:N(2),:, selected(i));
    %before1=clock;
    imThis2000     = iradon(squeeze(projThis2000),thetas)%,'linear','none',output_size);
    %timeRad2000= etime(clock,before1);
    figure; imagesc(imThis2000); colormap gray; colorbar; 
    caxis([-0.010 0.001]); 
    title(['Slice ' num2str(1000 + selected(i),'%04.0f')]); drawnow; %selected(i)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split Brebman

N1 = size(imThis2000);
sizeProj = size(projThis2000);
% Define Fwd and adjoint operators
A           = @(x)(radon(x,thetas, sizeProj(1)));
AT          = @(x)(iradon(x,thetas,'linear','None',1.0, N1(1)));%%%%,output_size CHANGED

flagNorm    = 0;


% N = [output_size output_size];
% if flagNorm == 1
%     tmp1        = ones(N);
%     Atmp1       = A(tmp1);
%     tmp2        = reshape(1:prod(N),N);
%     Atmp2       = A(tmp2);
%     AAtmp2      = AT(Atmp2);
%     a           = abs(tmp1(:)'*AAtmp2(:)/(Atmp1(:)'*Atmp2(:)));
%     b           = norm(tmp2(:))/norm(AAtmp2(:));
%     scaleA      = sqrt(a*b);
%     scaleAT     = sqrt(b/a);
%     
%     A           = @(x)(A(x)*scaleA);
%     AT          = @(x)(AT(x)*scaleAT);
%     
%     % % Test adjoint definition
%     % tmp1        = ones(N);
%     % Atmp1       = J(tmp1);
%     % tmp2        = reshape(1:prod(N),N);
%     % Atmp2       = J(tmp2);
%     % AAtmp2      = JT(Atmp2);
%     % abs(tmp1(:)'*AAtmp2(:))/(Atmp1(:)'*Atmp2(:))
%     % norm(tmp2(:))/norm(AAtmp2(:))
%     
%     uTarget         = uTarget/scaleA;
% end

% Recover the image
N           = size(projAll);
mu          = 1;
lambda      = 1;
gamma       = 1;
alpha       = 0.1;
nInner      = 1;
nBreg       = 200;
arr_idx     = ones(N);
rotAxis     = 1043.5;
output_size         = N1(1);

projThis2000(1:N(2),:)= projAll(1:N(2),:, 1);

[u,errAll, ctrst]  = TV_SB_2D_Gen(A,AT,projThis2000,[output_size output_size],mu, lambda, gamma, alpha, nInner, nBreg,imThis2000);




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % debug add
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% addpath('/home/goubet/Documents/data/FilesToReadData/fundebug');
% 
% %function [u,err] = TV_SB_2D_Gen(A,AT,f, N,mu, lambda, gamma, alpha, nInner, nBreg,varargin)
% f = p400;
% N = [output_size output_size];
% varargin = ref400;
% 
% tolKrylov   = 1e-2;
% 
% flagNorm    = 1;
% flagDisplay = 1;
% 
% % Normalize data
% %normFactor  = getNormalizationFactor(f,N(1));
% normFactor = 1/norm(f(:)/N(1));
% %  normFactor  = 1
% f           = normFactor*f;
% 
% % Stopiing criteria
% errMin      = 5e-4; 
% errPrevMin  = 1e-4;
% 
% % switch nargin
% %     case 11
%          uTarget     = varargin;
%          uTarget     = normFactor*uTarget;
% %     case 12
% %         uTarget     = varargin{1};
% %         uTarget     = normFactor*uTarget;
% %         optSolver   = varargin{1};
% %         tolKrylov   = optSolver.tolKrylov;
% %     case 13
% %         uTarget     = varargin{1};
% %         uTarget     = normFactor*uTarget;
% %         optSolver   = varargin{2};
% %         if isempty(optSolver) == 0
% %             tolKrylov   = optSolver.tolKrylov;
% %         end
% %         flagNorm    = varargin{3};
% %     case 14
% %         uTarget     = varargin{1};
% %         uTarget     = normFactor*uTarget;
% %         optSolver   = varargin{2};
% %         if isempty(optSolver) == 0
% %             tolKrylov   = optSolver.tolKrylov;
% %         end
%        % flagNorm    = varargin{3};
%        % stopCrit    = varargin{4};
%        % errMin      = stopCrit.errMin; 
%        % errPrevMin  = stopCrit.errPrevMin;
% %end % nargin
% 
% errAll      = zeros(nBreg,1);
% obj0        = zeros(nBreg,1);
% rate        = zeros(nBreg,1);
% errPrev     = zeros(nBreg,1);
% nnzx        = zeros(nBreg,1);
% 
% % Reserve memory for the auxillary variables
% rows        = N(1);
% cols        = N(2);
% f0          = f;
% u           = zeros(rows,cols);
% uPrev       = u;
% x           = zeros(rows,cols);
% y           = zeros(rows,cols);
% 
% bx          = zeros(rows,cols);
% by          = zeros(rows,cols);
% 
% % Normalize the forward and adjoint operators to comply with adjoint
% % operator definition and |J'Jx|~|x|
% if flagNorm == 1
%     tmp1        = ones(N);
%     Atmp1       = A(tmp1);
%     tmp2        = reshape(1:prod(N),N);
%     Atmp2       = A(tmp2);
%     AAtmp2      = AT(Atmp2);
%     a           = abs(tmp1(:)'*AAtmp2(:)/(Atmp1(:)'*Atmp2(:)));
%     b           = norm(tmp2(:))/norm(AAtmp2(:));
%     scaleA      = sqrt(a*b);
%     scaleAT     = sqrt(b/a);
%     
%     A           = @(x)(A(x)*scaleA);
%     AT          = @(x)(AT(x)*scaleAT);
%     
%     % % Test adjoint definition
%     % tmp1        = ones(N);
%     % Atmp1       = J(tmp1);
%     % tmp2        = reshape(1:prod(N),N);
%     % Atmp2       = J(tmp2);
%     % AAtmp2      = JT(Atmp2);
%     % abs(tmp1(:)'*AAtmp2(:))/(Atmp1(:)'*Atmp2(:))
%     % norm(tmp2(:))/norm(AAtmp2(:))
%     
%     uTarget         = uTarget/scaleA;
% end
% 
% murf        = mu*AT(f);
% outer       = 0;
% flagIt      = 1;
% uBest       = u;
% errBest     = inf;
% uDP         = u;
% errDPBest   = inf;
% iterDP      = 1;
% it          = 1;
% %  Do the reconstruction
% for outer = 1:nBreg
%     for inner = 1:nInner;
%         % update u
%         rhs         = murf+lambda*Dxt(x-bx)+lambda*Dyt(y-by)+gamma*u;
%         
%         u           = reshape(krylov(rhs(:)),N);
%         
%         dx          = Dx(u);
%         dy          = Dy(u);
%         
%         % update x and y
%         [x,y]       = shrink2(dx+bx,dy+by,alpha/lambda);
%         
%         % update bregman parameters
%         bx          = bx+dx-x;
%         by          = by+dy-y;                
%     end   % inner loop
%     
%     fForw           = A(u);
%     f               = f + f0-fForw;
%     murf            = mu*AT(f);
%     
%     nnzx(it)        = nnz(x(:));
%     if nargin >= 11
% %         if and(outer == 1, inner ==1)
% %             it = it + 1;
% %         end
%         it          = outer;
%         
%         % Solution error norm, data fidelity and rate of convergence
%         errAll(it)    = norm(uTarget(:)-u(:))/norm(uTarget(:));
%         %obj0(it)      = 0.5*norm(reshape(A(u)-f0,[],1));
%         obj0(it)      = norm(reshape((A(u)-f0)/normFactor,[],1));        
%         rate(it)      = norm(u(:)-uTarget(:))/norm(uTarget(:)-uPrev(:));
%         errPrev(it)   = norm(u(:)-uPrev(:))/norm(uPrev(:)+1e-2);
%         
%         if flagDisplay == 1 && any([it ==1, it == 10, rem(it, 25) == 0])
%             close;
%             h=figure;
%             subplot(3,2,1);
%             imagesc(abs(u)); title(['u, iter. ' num2str(it)]); colorbar;
%             subplot(3,2,2);
%             imagesc(abs(x));  colorbar; title('x');
%             subplot(3,2,3);
%             plot(errAll(1:it)); axis tight; title(['Sol. error' ]);
%             colormap gray;
%             subplot(3,2,4);
%             plot(obj0(1:it)); axis tight; title(['obj0' ]);
%             subplot(3,2,5);
%             plot(rate(1:it)); axis tight; title(['Rate' ]);
%             subplot(3,2,6);
%             plot(errPrev(1:it)); axis tight; title(['Err prev' ]);
%             axis([0 it 0 1]);
%             drawnow;
%         end % rem
%     end % nargin
%     
%     uPrev       = u;
% end % outer
% 
% if flagDisplay == 1 && nargin >= 11
%     close;
%     h=figure;
%     subplot(3,2,1);
%     imagesc(abs(u)); title(['u, iter. ' num2str(it)]); colorbar;
%     subplot(3,2,2);
%     imagesc(abs(x));  colorbar; title('x');
%     subplot(3,2,3);
%     plot(errAll(1:it)); axis tight; title(['Sol. error' ]);
%     colormap gray;
%     subplot(3,2,4);
%     plot(obj0(1:it)); axis tight; title(['obj0' ]);
%     subplot(3,2,5);
%     plot(rate(1:it)); axis tight; title(['Rate' ]);
%     subplot(3,2,6);
%     plot(errPrev(1:it)); axis tight; title(['Err prev' ]);
%     axis([0 it 0 1]);
%     drawnow;    
% end
% 
% if  nargin >= 11
%     err.errAll  = errAll;
%     err.errPrev = errPrev;
%     err.rate    = rate;
%     err.obj0    = obj0;
%     err.nnzx    = nnzx;
% end
% err.x       = x;
% err.y       = y;
% 
% % undo the normalization so that results are scaled properly
% % if nargin >= 11
% %     u = uBest/normFactor;
% % else
% %     u = u/normFactor;
% % end
% u = u*(scaleA*scaleAT)/(normFactor); % Check if correct!!!


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % end debug add
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 







figure; 
imagesc(u); colorbar; axis image; colormap gray; 

% %%%%% without zero padding
% 
% N           = size(projAll);
% N(2)/2;
% frequency_scaling   = 0.7;
% output_size         = N(1);
% 
% 
% %2000 projection
% numProj     = 2000;
% theta       = (1:numProj)*180/numProj;
% 
% for i= 1:4
%     projThis2000            = zeros([output_size numProj]);
%     projThis2000(1:N(2),:)= projAll(1:N(2),:, selected(i));
%     %before1=clock;
%     imThis2000      = iradon(squeeze(projThis2000),theta)%,'linear','none',output_size);
%     %timeRad2000= etime(clock,before1);
%     figure; imagesc(imThis2000); colormap gray; colorbar; 
%     caxis([-0.010 0.001]); drawnow;
%     title(['Slice ' num2str(1000 + selected(i),'%04.0f')]); %selected(i)]);
% end
% 
% 
% % with zero padding
% N           = size(projAll);
% rotAxis     = 1043.5;
% N(2)/2;
% 
% % FBP
% output_size         = round(2*rotAxis)%N(1);
% frequency_scaling   = 0.7;
% numProj     = 2000;
% theta       = (1:numProj)*180/numProj;
% for i= 1:4
%     projThis2000            = zeros([output_size numProj]);
%     projThis2000(1:N(2),:)= projAll(1:N(2),:, selected(i));
%     %before1=clock;
%     imThis2000      = iradon(squeeze(projThis2000),theta)%,'linear','none',output_size);
%     %timeRad2000= etime(clock,before1);
%     figure; imagesc(imThis2000); colormap gray; colorbar; 
%     caxis([-0.010 0.001]); drawnow;
%     title(['Slice ' num2str(1000 + selected(i),'%04.0f')]); %selected(i)]);
% end
% 
% 
% %%%%%focus on one slice and multiple freq scal
% 
% frequency_scaling = [1;.7;.5;.2;.1];
% 
% % with zero padding
% N           = size(projAll);
% rotAxis     = 1043.5;
% N(2)/2;
% 
% % FBP
% output_size         = round(2*rotAxis)%N(1);
% numProj     = 2000;
% theta       = (1:numProj)*180/numProj;
% for i= 4:5
%     projThis2000            = zeros([output_size numProj]);
%     projThis2000(1:N(2),:)= projAll(1:N(2),:, selected(1));
%     %before1=clock;
%     %add freq in the iradon
%     imThis2000      = iradon(squeeze(projThis2000),theta, 'linear',frequency_scaling(i));
%     %timeRad2000= etime(clock,before1);
%     figure; imagesc(imThis2000); colormap gray; colorbar; 
%     caxis([-0.0150 .001]); drawnow;
%     title(['Frequency scale ' sprintf('%03.0g', frequency_scaling(i))]); %selected(i)]);
% 
% end
% 
% 
% %%%%%focus on one slice and reduce change interpolation
% 
% selected_interpolation = {'linear'; 'nearest'; 'spline'; 'pchip'; 'v5cubic'};
% 
% N           = size(projAll);
% rotAxis     = 1043.5;
% N(2)/2;
% 
% % FBP
% output_size         = round(2*rotAxis)%N(1);
% frequency_scaling   = 0.7;
% numProj     = 2000;
% theta       = (1:numProj)*180/numProj;
% time = zeros(5);
% for i= 1:5
%     projThis2000            = zeros([output_size numProj]);
%     projThis2000(1:N(2),:)= projAll(1:N(2),:, selected(1));
%     before1=clock;
%     imThis2000      = iradon(squeeze(projThis2000),theta, selected_interpolation{i});
%     time(i)= etime(clock,before1);
%     figure; imagesc(imThis2000); colormap gray; colorbar; 
%     caxis([-0.0150 .001]); drawnow; 
%     title(selected_interpolation{i});
% end
% 
% 
% 
% 
% 
% %%%%%focus on one slice and reduce change filter
% 
% selected_filter = {'None'; 'Hann'; 'Hamming'; 'Cosine'; 'Shepp-Logan'; 'Ram-Lak'};
% 
% N           = size(projAll);
% rotAxis     = 1043.5;
% N(2)/2;
% 
% % FBP
% output_size         = round(2*rotAxis)%N(1);
% frequency_scaling   = 0.7;
% numProj     = 2000;
% theta       = (1:numProj)*180/numProj;
% time = zeros(5);
% for i= 1:6
%     projThis2000            = zeros([output_size numProj]);
%     projThis2000(1:N(2),:)= projAll(1:N(2),:, selected(1));
%     before1=clock;
%     imThis2000      = iradon(squeeze(projThis2000),theta, selected_filter{i});
%     time(i)= etime(clock,before1);
%     figure; imagesc(imThis2000); colormap gray; colorbar; 
%     caxis([-0.0150 .001]); drawnow; 
%     title(selected_filter{i});
% end
% 
% 
% 
% 
% %%%%%focus on 1 slice and reduce number of projections
% 
% %2000
% 
% numProj     = 2000;
% theta       = (1:numProj)*180/numProj;
% projThis2000            = zeros([round(2*rotAxis) numProj]);
% projThis2000(1:N(2),:)= projAll(1:N(2),:, selected(1));
% before1=clock;
% imThis2000      = iradon(squeeze(projThis2000),theta)%,'linear','none',output_size);
% timeRad2000= etime(clock,before1);
% figure; imagesc(imThis2000); colormap gray; colorbar; 
% title('2000 projections');
%  caxis([-0.0150 .001]); drawnow; 
% 
% save(fullfile(pathSave,[nameMainSave '_img_2000v_Slice1000']),'imThis2000','-mat');
%     
% %1000
% 
% numProj     = 1000;
% theta       = (1:numProj)*180/numProj;
% projThis1000 = zeros([round(2*rotAxis) numProj]);
% projThis1000(1:N(2),:) = projAll(1:N(2),1:2:2000, selected(1));
% before1=clock;
% imThis1000      = iradon(squeeze(projThis1000),theta);%,'linear','none',output_size);
% timeRad1000= etime(clock,before1);
% figure; imagesc(imThis1000); colormap gray; colorbar; 
% caxis([-0.0150 0.001]); drawnow;
% title('1000 projections');
% peaksnr1000 = psnr(imThis2000,imThis1000) 
% 
% save(fullfile(pathSave,[nameMainSave '_img_1000v_Slice1000']),'imThis1000','-mat');
% 
% 
% %400
% numProj     = 400;
% theta       = (1:numProj)*180/numProj;
% projThis400 = zeros([round(2*rotAxis) numProj]);
% projThis400(1:N(2),:) = projAll(1:N(2),1:5:2000, selected(1));
% before1=clock;
% imThis400      = iradon(squeeze(projThis400),theta);%,'linear','none',output_size);
% timeRad400= etime(clock,before1);
% figure; imagesc(imThis400); colormap gray; colorbar; 
% caxis([-0.0150 0.001]); drawnow;
% title('400 projections');
% peaksnr400 = psnr(imThis2000,imThis400) 
% 
% save(fullfile(pathSave,[nameMainSave '_img_400v_Slice1000']),'imThis400','-mat');
% 
% 
% %200
% numProj     = 200;
% theta       = (1:numProj)*180/numProj;
% projThis200 = zeros([round(2*rotAxis) numProj]);
% projThis200(1:N(2),:) = projAll(1:N(2),1:10:2000, selected(1));
% before1=clock;
% imThis200      = iradon(squeeze(projThis200),theta);%,'linear','none',output_size);
% timeRad200= etime(clock,before1);
% figure; imagesc(imThis200); colormap gray; colorbar; 
% caxis([-0.0150 0.001]); drawnow;
% title('200 projections');
% peaksnr200 = psnr(imThis2000,imThis200) 
% 
% save(fullfile(pathSave,[nameMainSave '_img_200v_Slice1000']),'imThis200','-mat');
% 
% %100
% numProj     = 100;
% theta       = (1:numProj)*180/numProj;
% projThis100 = zeros([round(2*rotAxis) numProj]);
% projThis100(1:N(2),:) = projAll(1:N(2),1:20:2000, selected(1));
% before1=clock;
% imThis100      = iradon(squeeze(projThis100),theta);%,'linear','none',output_size);
% timeRad100= etime(clock,before1);
% figure; imagesc(imThis100); colormap gray; colorbar; 
% caxis([-0.0150 0.001]); drawnow;
% title('100 projections');
% peaksnr100 = psnr(imThis2000,imThis100) 
% 
% save(fullfile(pathSave,[nameMainSave '_img_100v_Slice1000']),'imThis100','-mat');
% 
% 
% 
% 
% 
% 
% 
% %%%% different filters 100 proj
% 
% selected_filter = {'None'; 'Hann'; 'Hamming'; 'Cosine'; 'Shepp-Logan'; 'Ram-Lak'};
% 
% N           = size(projAll);
% rotAxis     = 1043.5;
% N(2)/2;
% 
% % FBP
% 
% numProj     = 100;
% theta       = (1:numProj)*180/numProj;
% projThis100 = zeros([round(2*rotAxis) numProj]);
% projThis100(1:N(2),:) = projAll(1:N(2),1:20:2000, selected(1));
% 
% output_size         = round(2*rotAxis);%N(1);
% frequency_scaling   = 0.7;
% time = zeros(5);
% for i= 1:6
%     before1=clock;
%     imThis2000      = iradon(squeeze(projThis100),theta, selected_filter{i});
%     time(i)= etime(clock,before1);
%     figure; imagesc(imThis2000); colormap gray; colorbar; 
%     caxis([-0.0150 .001]); drawnow; 
%     title(selected_filter{i});
% end
% 
% 
% %%%%% different rec on 100 proj
% 
% selected_interpolation = {'linear'; 'nearest'; 'spline'; 'pchip'; 'v5cubic'};
% 
% numProj     = 100;
% theta       = (1:numProj)*180/numProj;
% projThis100 = zeros([round(2*rotAxis) numProj]);
% projThis100(1:N(2),:) = projAll(1:N(2),1:20:2000, selected(1));
% 
% output_size         = round(2*rotAxis);%N(1);
% frequency_scaling   = 0.7;
% time = zeros(5);
% for i= 1:5
%     before1=clock;
%     imThis2000      = iradon(squeeze(projThis100),theta, selected_interpolation{i});
%     time(i)= etime(clock,before1);
%     figure; imagesc(imThis2000); colormap gray; colorbar; 
%     caxis([-0.0150 .001]); drawnow; 
%     title(selected_interpolation{i});
% end
% 
% 
% 
% % load projections 
% % load(fullfile(pathSave,[nameMainSave '_img_100v_Slice1000']), '-mat');
% % load(fullfile(pathSave,[nameMainSave '_img_200v_Slice1000']), '-mat');
% % load(fullfile(pathSave,[nameMainSave '_img_400v_Slice1000']), '-mat');
% % load(fullfile(pathSave,[nameMainSave '_img_1000v_Slice1000']), '-mat');
% % load(fullfile(pathSave,[nameMainSave '_img_2000v_Slice1000']), '-mat');
% 
% peaksnr100 = psnr(imThis100,imThis2000);
% peaksnr200 = psnr(imThis200,imThis2000);
% peaksnr400 = psnr(imThis400,imThis2000);
% peaksnr1000 = psnr(imThis1000,imThis2000);
% peaksnr2000 = psnr(imThis2000,imThis2000);
% 
% ssim100 = ssim(imThis100,imThis2000);
% ssim200 = ssim(imThis200,imThis2000);
% ssim400 = ssim(imThis400,imThis2000);
% ssim1000 = ssim(imThis1000,imThis2000);
% ssim2000 = ssim(imThis2000,imThis2000);
% 
% ssim100 = immse(imThis100,imThis2000);
% ssim200 = immse(imThis200,imThis2000);
% ssim400 = immse(imThis400,imThis2000);
% ssim1000 = immse(imThis1000,imThis2000);
% ssim2000 = immse(imThis2000,imThis2000);
% 
% 
% %save(fullfile(pathSave,[nameMainSave '_img_100v_Slice1000']),'imThis100','-mat');
% %iterative
% %numProj     = 2000;
% %theta       = (1:numProj)*180/numProj;
% %projThis = zeros([round(2*rotAxis) numProj]);
% %projThis(1:N(2),:) = projAll(1:N(2),:);
% %for is = 1:N(2)
% %  if (mod(is,5) == 0) 
% %    projThis(1:N(2),is) = projAll(1:N(2),is);
% %  end
% %end
% %imThis      = iradon(squeeze(projThis),theta);%,'linear','none',output_size);
% %figure; imagesc(imThis); colormap gray; colorbar; 
% %caxis([-0.010 0.001]); drawnow;
% 
% 
% projSum     = squeeze(sum(projAll,2));
% imThis      = iradon(projSum,theta,'linear',frequency_scaling,output_size);
% figure; imagesc(imThis); colormap gray; colorbar; 
% 
% figure;
% for ip = 1:30%size(projAll,2)
%     projThis    = squeeze(projAll(:,ip,:));
%     imThis      = iradon(squeeze(projAll(:,ip,:)),theta,'linear',frequency_scaling,output_size);
%     imAll(:,:,ip) = imThis;
%     imagesc(imThis); colormap gray; colorbar; colormap gray; colorbar; axis image; title(num2str(ip)); drawnow;
% end
% % Loop across slices in z
% figure;
% for iz = 1:N(1)
%     sinoThis    = squeeze(projAll(iz,:,:));
%     imagesc(sinoThis); colormap gray; colorbar; drawnow;
%     imThis      = iradon(sinoThis,theta,'linear',frequency_scaling,output_size);
%     imagesc(imThis); colormap gray; colorbar; drawnow;
%     imAll(:,:,iz) = imThis;
% end

