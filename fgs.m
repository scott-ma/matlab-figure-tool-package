function r_out = fgs(varargin)
%mainly for setting fig parameters  
% Copyright 2009 lshuily@gmail.com
%  fgs('lxy',[xl xb yl yb])
%  fgs('pfg','fnm','ftp','fco')
% =============================================
if (nargin < 1)
  error('Too few input arguments');
end
swh = varargin{1};    % case

switch swh
% ---------------------------------------------
case 'lxy'   % fgs('lxy',[xl xb yl yb])
    if (nargin < 2)
        error('Too few input arguments');
    else
        opt = varargin{2};
        lxy(opt); % change xy label position
    end 
case 'pfg'   % pfg('fnm','ftp','fco')
    if nargin < 2
        error('Too few input arguments');
    else
        par = varargin(2:end);
        opt{1} = varargin{2};
        opt{2} = 'eps';
        opt{3} = 'c';
        if length(par)==2
             if strcmp(par{2},'tiff') % tiff
                 opt{2} = 'tiff';
             else
                 opt{3} = par{2};
             end
        elseif length(par)==3
             opt{2} = par{2};
             opt{3} = par{3};
        end
        pfg(opt);     % print fig 
        if strcmp(opt{2},'eps')
             ltp(opt{1});  % change line type
        end
    end   
end % switch


%
%  Local Functions
% ----------------
function lxy(xyp) 
% change xy label position
% set
a = gca;
hxlbl = get(a,'xlabel');
hylbl = get(a,'ylabel');
xylunits = get(hxlbl,'units');
% xlbl
xlblv = get(hxlbl,'string');
if ~isempty(xlblv)
    set(hxlbl,'units','points');
    posxl = get(hxlbl,'position');
    posxl = posxl(1:2)+xyp(1:2);
    set(hxlbl,'position',posxl);
end
% ylbl
ylblv = get(hylbl,'string');
if ~isempty(ylblv)
    set(hylbl,'units','points');
    posyl = get(hylbl,'position');
    posyl = posyl(1:2)+xyp(3:4);
    set(hylbl,'position',posyl);
end
set([hxlbl hylbl],'units',xylunits);

% ----------------
function pfg(opt) 
% print fig 
fnm = opt{1};
ftp = opt{2};
fco = opt{3};
if strcmp(ftp,'eps')      % eps
    if strcmp(fco,'k')
        printargs = {'-deps2', '-loose', '-r300'};
    elseif strcmp(fco,'c')
        printargs = {'-depsc2', '-loose', '-r300'};
    elseif strcmp(fco,'ki')
        printargs = {'-deps', '-loose', '-r300'};
    elseif strcmp(fco,'ci')
        printargs = {'-depsc', '-loose', '-r300'};
    end
elseif strcmp(ftp,'tiff') % tiff
    if strcmp(fco,'k')
        set(findobj(gca,'Type','line'),'Color','k');
    end
    printargs = {'-dtiff', '-loose', '-r300'};
end
print(gcf,fnm,printargs{:});

% ----------------
function ltp(filename) 
% change line type of fig
stra = '/SO { [] 0 setdash } bdef';
str1 = '/DO { [ 1.5 dpi2point mul 2 dpi2point mul] 0 setdash} bdef';
str2 = '/DA { [ 5 dpi2point mul 2 dpi2point mul] 0 setdash} bdef';
str3 = '/DD { [ 1 dpi2point mul 2 dpi2point mul 4 dpi2point mul 2';
figname = strcat(filename,'.eps');
tmpname = strcat('t',figname);
copyfile(figname,tmpname);
fwd = fopen(figname,'wt');
frd = fopen(tmpname,'rt');
while ~feof(frd)        % end of file?
    tline = fgetl(frd); % read a line
    x = strmatch(stra,tline,'exact'); % match
    if ~isempty(x)      % the line of stra
        break; 
    else
        fprintf(fwd,'%s',tline); 
        fprintf(fwd,'\n'); 
    end
end
fprintf(fwd,'%s',tline); % stra
fprintf(fwd,'\n'); 
fprintf(fwd,'%s',str1);  % str1
fprintf(fwd,'\n'); 
fprintf(fwd,'%s',str2);  % str2
fprintf(fwd,'\n'); 
fprintf(fwd,'%s',str3);  % str3
fprintf(fwd,'\n'); 
tline = fgetl(frd); % read a line str1
tline = fgetl(frd); % read a line str2
tline = fgetl(frd); % read a line str3
while ~feof(frd)        % end of file?
    tline = fgetl(frd); % read a line
    fprintf(fwd,'%s',tline); 
    fprintf(fwd,'\n'); 
end
fclose(frd);
fclose(fwd);
delete(tmpname);

