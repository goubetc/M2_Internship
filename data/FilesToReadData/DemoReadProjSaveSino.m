

cd('/home/goubet/Documents/data/FilesToReadData');

nameMain    = '../multips_166_2012_femR_1L_120nm_tomo1_';
pathSave    = '/home/goubet/Documents/data/FilesToReadData';
nameMainSave    = 'multips_166_2012_femR_1L_120nm_stack';

if 1
    numProj     = 2000;
    Np          = [2048 2048 numProj];
    numPStack   = 100;
    indLoop     = 1:numPStack:Np(1);
    numStack    = length(indLoop);
    for is = 1:numStack-1
        indSelect   = indLoop(is):indLoop(is+1)-1;
        sinoThis    = zeros(Np(1),numPStack,Np(3));
        for ip = 1:Np(3)
            %         ind       = ind + 1;
            filename  = [nameMain num2str(ip,'%04.0f')];
            %   im        = edfread([filename '.edf']);
            data=load(fullfile(pathSave,filename),'im');
            %sinoThis(:,:,ip)  = squeeze(data.im(:,indSelect)); %data.im(:,indSelect);
            %         imagesc(projAll(:,:,ip)); colormap gray; colorbar; axis image; title(num2str(ip)); drawnow;
        end
        save(fullfile(pathSave,[nameMainSave '_' num2str(is)]),'sinoThis');
    end    
end

figure; imagesc(squeeze(projAll(:,1,:))); axis image; colormap gray; colorbar; drawnow;
caxis([-0.06 0.06]);
    
