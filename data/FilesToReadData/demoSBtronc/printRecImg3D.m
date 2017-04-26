function img = printRecImg3D(targetId, numP, nbIt, recImages, ImgSize, nbProj, nbNBreg, imgNum)
    img         =zeros(ImgSize, ImgSize);
    img(:,:)    = recImages(targetId,numP,nbIt,imgNum,:,:);
    figure; imagesc(img); colormap gray;colormap(flipud(colormap)); colorbar; 
        caxis('auto'); 
        title(['Reconstructed image from target ' num2str(targetId) ', nbProj: ' num2str(nbProj(numP)) ', nbIt: ' num2str(nbNBreg(nbIt))] ); drawnow; 
end