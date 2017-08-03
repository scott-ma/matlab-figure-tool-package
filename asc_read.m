function dat = asc_read(varargin)
%read ascii file with data head
%  asc_read('filename',start_line)
% =============================================
if (nargin < 1)
    error('Too few input arguments');
elseif (nargin == 1)
    start = 1;
else
    start = varargin{2};
end
fname = varargin{1};

% open file
frd   = fopen(fname,'rt');
if frd==-1
    disp('File does not exist'); 
end

% read file
dat = [];
lth = 0;
count = 0;
nline = 1;
while ~feof(frd)        % end of file?
    count = count+1;
    tline = fgetl(frd);   % read a line
    dat0  = str2num(tline);
    xx = length(dat0);
    if (count>start && length(dat0)==lth)
        dat(nline,1:lth) = dat0;
        nline = nline+1;
    elseif (count==start)
        lth  = length(dat0);
        dat(nline,1:lth) = dat0;
        nline = nline+1;
    end
end

fclose(frd);
