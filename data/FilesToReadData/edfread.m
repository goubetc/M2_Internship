% edfmread
% function im=edfmread(filename,varargin)
%      reads an image in edf format into the matrix im
%      uses pmedf_read by Petr Mikulik if only filename is given
%      if both header and image are needed, use [hd,im]=pmedf_read(filename)
%      does not transpose the image, 1st dimension remains 1st dimension of original file
%      
%      arguments:
%      argument 1: filename
%      argument 2 (optional): vector of adjacent lines to be read
%      argument 3 (optional): vector of adjacent colums to be read
%      argument 4 (optional): vector of layers to be read
%      
%      examples:
%      im = edfread('dark.edf');
%          reads complete image
%      im = edfread('dark.edf',10:200,30:50);
%          reads subimage corresponding to row (Dim_1) 10 until 200 and column (Dim_2) 30 until 50
%      im = edfread('dark.edf',10:200,30:50,[1 4 5]);
%          reads layers 1, 4 and 5 of a 3-dimensional image (see Dim_3 and W. Ludwig)
%
%      See also: edfwrite, pmedf_read, pmedf_write

% Author: P. Cloetens <cloetens@esrf.fr>
%
% 2007-03-03 P. Cloetens <cloetens@esrf.fr>
% * Initial revision called edfread
% 2015-06-16 W. Augustin <awerner@esrf.fr>
% * implementation of multi edf files
% * using actual header size value instead of fixed value 1024
% * current implementation still wrongly assumes same headersize for all images!
% 2016-03-07 C. Nemoz <nemoz@esrf.fr>
% * rename into edfread
% * add SignedInteger to deal with Lima accumulation mode

function im=edfread(filename,varargin)

usepm = 0;

switch nargin
    case 1
        usepm = isoctave&exist('pmedf_read');
    case 3
        lines = varargin{1};
        columns = varargin{2};
    case 4
        lines = varargin{1};
        columns = varargin{2};
        layers =varargin{3};  
end


if usepm
    [hd,im] = pmedf_read(filename);
else
    fid=fopen(filename,'r');

    if fid==-1
	    im=fid;
	    disp(sprintf('!!! Error opening file %s !!!',filename))
    else
	    hd=readedfheader(fid);

        byteorder=findheader(hd,'ByteOrder','string');
        headersize=ftell(fid); % store present position in file
        fclose(fid); 
        if strcmp(byteorder,'HighByteFirst')
            byteorder='b';
        else
            byteorder='l';
        end
        fid=fopen(filename,'r',byteorder); % re-open with good format
        fseek(fid,headersize,0); % re-position at end of header

        xsize=findheader(hd,'Dim_1','integer');
	    ysize=findheader(hd,'Dim_2','integer');
        %zsize=findheader(hd,'Dim_3','integer');
        %if isempty(zsize)
            %zsize=1;
        %end
				if exist('layers','var')
					zsize=length(layers);
				end

        datatype=findheader(hd,'DataType','string');
	    switch datatype
		    case 'UnsignedByte'
			    datatype='uint8';
			    nbytes=1;
		    case {'UnsignedShort', 'UnsignedShort;', 'SignedShort'}
			    datatype='uint16';
			    nbytes=2;
		    case {'UnsignedInteger','SignedInteger','UnsignedLong'}
			    datatype='uint32';
			    nbytes=4;
		    case {'Float','FloatValue','FLOATVALUE','Real'}
			    datatype='float32';
			    nbytes=4;
		    case 'DoubleValue'
			    datatype='float64';
			    nbytes=8;
    %etcetera
	    end

	    if ~exist('lines','var') || isempty(who('lines'))
		    lines=1:xsize;
	    end
	    if ~exist('columns','var') || isempty(who('columns'))
		    columns=1:ysize;
			end
        %if isempty(who('layers'))
					%layers=1:zsize;
        %end

        if 0 % zsize==1
     	    fseek(fid,nbytes*(lines(1)-1+(columns(1)-1)*xsize),0);
	        im=fread(fid,[length(lines),length(columns)],sprintf('%d*%s',length(lines),datatype),nbytes*(xsize-length(lines)));
        else
					j=1;
					for i=layers
					%disp('blah')
							%orig:fseek(fid,headersize+nbytes*(lines(1)-1+(columns(1)-1)*xsize+xsize*ysize*(i-1)),-1);
							fseek(fid,headersize*i+nbytes*(xsize*ysize*(i-1)),'SEEK_SET'); % SEEK_SET not ok!
							%fseek(fid,headersize*i+nbytes*(xsize*ysize*(i-1)),-1);
							%fseek(fid,headersize*i+nbytes*(lines(1)-1+(columns(1)-1)*xsize+xsize*ysize*(i-1)),'SEEK_SET');
							% absolute (ok up to 2GB file size):
							%fseek(fid,headersize*i+nbytes*(lines(1)-1+(columns(1)-1)*xsize+xsize*ysize*(i-1)),-1);
							% relative (OK):
							if j==1
								fseek(fid,headersize*i+nbytes*xsize*ysize*(i-1),-1);
							else	
								fseek(fid,headersize*(i-layers(j-1))+nbytes*xsize*ysize*(i-layers(j-1)-1),0); %'SEEK_CUR' doesn't work either!!!
							end
							%
							%ftell(fid)  % for debugging
							%im(:,:,j)=fread(fid,[length(lines),length(columns)],sprintf('%d*%s',length(lines),datatype),nbytes*(xsize-length(lines)));        
							im(:,:,j)=fread(fid,[xsize,ysize],sprintf('%d*%s',length(lines),datatype),nbytes*(xsize-length(lines)));        
							j=j+1;
					end
			end    
			st=fclose(fid);

	end
end % if usepm
