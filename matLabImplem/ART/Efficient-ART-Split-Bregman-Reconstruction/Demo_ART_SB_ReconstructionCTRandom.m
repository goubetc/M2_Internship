%Shepp-Logan phantom
N   = 64;
uin = phantom('Modified Shepp-Logan', N);
figure, imagesc(uin), colorbar, colormap gray

Nthetas=100;
thetas = 1:floor(360/Nthetas):360;

numIterART = 200;
relaxParam = 0.5;

%experimental projections
d           = radon(uin, thetas);
% indZeros    = randperm(Nthetas);
% d(:,indZeros(1:200)) = 0;
figure, imagesc(d), colorbar, colormap gray

NData   = prod(size(d));

uRetro      = iradon(d, thetas,'linear','none',N);
figure, imagesc(uRetro), colorbar, colormap gray

uFBP        = iradon(d, thetas,'spline','Ram-Lak',1,N);
figure, imagesc(uFBP), colorbar, colormap gray

% u   = zeros(N);
u   = uFBP;
% Ru  = radon(u, thetas);
% figure, imagesc(abs(d-Ru)), colorbar, colormap gray

close all;
% ---
% img = phantom(64);
% bw = double(im2bw(img));
h = figure; 
% Loop across iterations
for it = 1:numIterART
    indData     = 0;
    indDataRd   = randperm(NData);
    % Loop across angles
    for id = 1:NData
%      for ir = 1:Nthetas
         if rem(id,10000) == 0
             figure(h); imagesc(u); colormap gray; colorbar; axis image; 
             title(['Iter ' num2str(it) ', \theta=' num2str(id)]);
             drawnow;
         end
        % Loop across u
%          for i = 1:size(d,1)
            %for j = 1:64
            %bw = zeros(64,64);
            %bw(i,:)=1;
            bw = zeros(size(d));
            bw(indDataRd(id))=1;
            [temp, thetaThis] = find(bw);
            Ai      = iradon(bw, thetaThis,'linear','none',N);  
            if norm(Ai(:)) ~= 0
                ratio   = (d(indDataRd(id))-Ai(:)'*u(:))/(Ai(:)'*Ai(:));
                uUpdate = relaxParam*ratio*Ai;
                u       = u + uUpdate;
            end            
    %             Ai = radon(bw, thetas(ir));
%                 Ai = radon(bw, thetas(ir));
%                 d2 = radon(u(:,:), thetas(ir));
                %Air = radon(bw(:,:), thetas(ir));
%                 num         = d(ir)-d2;
                %u(i,j)           = u(i,j) + mtimes((relaxParam*num/((Ai') * Ai)),Ai);
%                 u           = u + mtimes((.9*num/((bw(i,:)') * bw(i,:))),bw(i,:));
    end %i
        
%     end % ir
end % it
% ---


urec = zeros(64,64);
urec = ARTReconstructionCT(thetas, d, .9, 14, urec);
plot(urec);
delta1 = zeros(64,64);
delta1(1) = 1;


A1 = radon(delta1, 0);

test = A1 * u(1,1);

