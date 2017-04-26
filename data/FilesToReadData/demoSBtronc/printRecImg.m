function img = printRecImg(targetId, numP, nbIt, recImages, ImgSize, nbProj, nbNBreg)
    img         =zeros(ImgSize);
    img(:,:)    = recImages(targetId,numP,nbIt,:,:);
    figure; imagesc(img); colormap gray;colormap(flipud(colormap)); colorbar; 
        caxis('auto'); 
        title(['Reconstructed image from target ' num2str(targetId) ', nbProj: ' num2str(nbProj(numP)) ', nbIt: ' num2str(nbNBreg(nbIt))] ); drawnow; 
end