function [errors, ePsnr, tv] = plotErrorProj(target, iter, err, numProj, nBreg, recImages, targets, numImg)
    lengthErros = length(numProj);
    errors = zeros(lengthErros,1);
    ePsnr = zeros(lengthErros,1);
    tv = zeros(lengthErros,1);
    %sumDiff = zeros(lengthErros,1);
    for i = 1:lengthErros
        errImg = err{target, i, iter}{1, 1};
        errors(i) = errImg(length(errImg));
        %img = zeros(size(targets{1,target}));
        img(:,:) = recImages(target,i,iter,numImg,:,:);
        imgTarg(:,:) = targets{1,target}(numImg,:,:);
        ePsnr(i) = psnr(imgTarg, img);
        % tv
        t(:,:) = targets{1,target}(numImg,:,:);
        sumDiff = t(:,:) - img(:,:);
        [gDiff_x, gDiff_y] = gradient(sumDiff);
        tv(i) = sum(sum(sqrt(gDiff_x.^2 + gDiff_y.^2)));
%         t = targets{1,target};
%         sumDiff(i) = t(:,:) - img(:,:);
    end
%     h = figure;
%     plot(numProj(:),errors(:),numProj(:), ePsnr(:)); axis tight; title(['Error image ' num2str(target) ' nbIter ' num2str(nBreg(iter)) ]);
%     legend(['Error rate', 'PSNR']);
%     drawnow;
   
   figure
    ax1 = gca;
    hold on
    plot(numProj(:),errors(:))  
    legend('Error');
    ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
    linkaxes([ax1 ax2],'x');
    hold on
    plot(numProj(:),ePsnr(:),'Parent',ax2); 
    legend('PSNR');
    
end