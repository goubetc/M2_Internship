function img = printRecImg3D(targetId, numP, recImages, ImgSize, nbProj, nbNBreg, imgNum)
    img         = squeeze(recImages(targetId,numP,imgNum,:,:));
    figure; imagesc(img); colormap gray;colormap(flipud(colormap)); colorbar; 
        caxis('auto'); 
        title(['Reconstructed image from target ' num2str(targetId) ', nbProj: ' num2str(nbProj(numP))] ); drawnow; 
end