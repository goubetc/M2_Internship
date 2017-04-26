
cd('/home/goubet/Documents/data/FilesToReadData');

nameMain    = '../multips_166_2012_femR_1L_120nm_tomo1_/multips_166_2012_femR_1L_120nm_tomo1_';
pathSave    = '/home/goubet/Documents/data/FilesToReadData';

cd(pathSave);
numProj     = 2000;
ind         = 0;
%projAll     = zeros(2048,2048,200);
for ip = 1:numProj
  ind       = ind + 1;
  filename  = [nameMain num2str(ip,'%04.0f')];
  im        = edfread([filename '.edf']);
  %projAll(:,:,ind)  = im;
  filename  = [filename '.mat']
  %save('-mat7-binary',filename,'im')
  save(filename,'im')
end