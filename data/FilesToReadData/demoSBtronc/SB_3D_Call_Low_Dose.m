function [u,errAll, ctrst, timeExec] = SB_3D_Call_Low_Dose(target, numProj, nBreg, noise)

    sizeImg= size(target);
    output_size = sizeImg(2);
    height = size(target, 1);

    thetas      = (1:numProj)*180/numProj;
    tmpImg(:,:) = target(1,:,:);
    p           = radon(tmpImg, thetas);
    N           = size(p);
    f           = radon3D(target);
    if nargin == 4
        %Add poisson noise
        epsilon = 5; % corrects for negative vlaues in the pre logged sinogram data that would otherwise results in infinite projection values
        Noise = noise*exp(-f);
        N_noise = poissrnd(Noise);
        fnoisy = -log(N_noise/noise);
        idx = isinf(fnoisy);
        fnoisy(idx) = -log(epsilon/noise);
    %     figure; imagesc(squeeze(f(1,:,:))); colormap gray;colormap(flipud(colormap)); colorbar; 
    %     caxis('auto'); 
    %     title('Projections with poisson noise '); drawnow;
    else 
        fnoisy = f;
    end

    % Define Fwd and adjoint operators
    A           = @(x)(radon3D(x));
    AT          = @(x)(iradon3D(x));%%%%,output_size CHANGED
    % 2D operators for normalization
    ANorm       = @(x)(radon(x,thetas, N(1)));
    ATNorm      = @(x)(iradon(x,thetas,'linear','None',1.0, output_size));
    
    uretro      = iradon(squeeze(f(1,:,:)),thetas);
    figure; imagesc(uretro); colormap gray;colormap(flipud(colormap)); colorbar; 

    uretro      = iradon(squeeze(fnoisy(1,:,:)),thetas);
    figure; imagesc(uretro); colormap gray;colormap(flipud(colormap)); colorbar; 

    flagNorm    = 0;

    mu          = 1;
    lambda      = 1;
    gamma       = 1e-2;
    alpha       = 1;
    nInner      = 1;
    %nBreg       = 1000;
    %arr_idx     = ones(N);

    utest = AT(A(target));
    before1=clock;
    [u,errAll, ctrst]  = TV_SB_3D_Gen(A,AT,ANorm,ATNorm,fnoisy,[size(target,1) output_size output_size],mu, lambda, gamma, alpha, nInner, nBreg,target);
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
