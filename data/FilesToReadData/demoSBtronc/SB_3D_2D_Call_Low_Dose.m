function [u,errAll, ctrst, timeExec, uBest] = SB_3D_2D_Call_Low_Dose(target, numProj, nBreg, noise, gamma, alpha, lambda, mu, f)

    sizeImg= size(target);
    output_size = sizeImg(2);
    height = size(target, 1);

    thetas      = (1:numProj)*180/numProj;
    tmpImg(:,:) = target(1,:,:);
    p           = radon(tmpImg, thetas);
    N           = size(p);
    if nargin < 9
        f       = radon3D(target);
    end
    if nargin >= 4 && noise ~= 0
        % add gaussian noise
         fnoisy = zeros(size(f));
        for im=1:size(f,1)
            d = squeeze(f(im,:,:));
%             d2 = d + max(d(:));
%             d2 = d2/min(nonzeros(d2(:)));
%             v = noise*var(d2(:));
%             dnoise = (imnoise(d2, 'gaussian', 0, v));
%             mins = min(nonzeros(d(:)));
%             fnoisy(im,:,:) = mins.* dnoise;
            fnoisy(im,:,:) = d+(noise*min(nonzeros(d(:)))*randn(size(d)));
        end
    else 
        fnoisy = f;
    end

    
    % Define Fwd and adjoint operators
    A           = @(x)(radon3D(x));
    AT          = @(x)(iradon3D(x));%%%%,output_size CHANGED
    % 2D operators for normalization
    ANorm       = @(x)(radon(x,thetas, N(1)));
    ATNorm      = @(x)(iradon(x,thetas,'linear','None',1.0, output_size));
    
    uretro      = iradon(squeeze(fnoisy(1,:,:)),thetas);
    figure; imagesc(uretro); colormap gray;colormap(flipud(colormap)); colorbar; 
    
    flagNorm    = 0;
    if nargin <= 4
        mu          = 5;
        lambda      = 5;
        gamma       = 0;
        alpha       = 0.05;
    end
    nInner      = 1;
    %nBreg       = 1000;
    %arr_idx     = ones(N);

    utest = AT(A(target));
    before1=clock;
    [u,errAll, ctrst, uBest]  = TV_SB_3D_Gen_2D(A,AT,ANorm,ATNorm,fnoisy,[size(target,1) output_size output_size],mu, lambda, gamma, alpha, nInner, nBreg,target);
    timeExec= etime(clock,before1);

    function d3D = radon3D(u)
        d3D = zeros([size(u,1) N]);
        %%for octave
        %numProj         = 156/10;
        %height          = 5;
        %thetas          = (1:numProj)*180/numProj;
        %tmp             = radon(u(1), thetas);
        %d3D             = zeros([height size(tmp)]);
        for i = 1:height
            imtmp(:,:) = u(i,:,:);
            d3D(i,:,:) = radon(imtmp,thetas);
        end
    end

    function u3D = iradon3D(p3D)
    %thetas      = (1:numProj)*180/numProj;
        u3D = zeros(size(target,1),size(target,2), size(target,3));
        % changed for octave
        %imtmp(:,:) = p3D(1,:,:);
        %tmp = iradon(imtmp,thetas,'linear','None');
        %u3D = zeros([height size(tmp)]);
        for i = 1:height
            imtmp(:,:) = p3D(i,:,:);
            u3D(i,:,:) = iradon(imtmp,thetas,'linear','None',1.0, output_size);
        end
    end
end
