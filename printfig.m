function varargout = printfig(varargin)
% set figure parameters and export as eps or tiff
%
%   EXPORTFIG(H, FILENAME, OPT) writes the figure H to FILENAME.  H is
%   a figure handle and FILENAME is a string that specifies the
%   name of the output file.  OPT is setting parameters.
%
%   'lbru' minor adjust
%  opts = struct('lbrt',[3 1 0 0],'figsize',[width hight delete],'tickxylblfs',[ticklblfs xylblfs],
%                'ticksize',[xfactor yfactor minortickornot],'xlbl',[x1 xstep xn],
%                'ylbl',[y1 ystep yn],'axis',[x1 xn y1 yn],'figcolor',[c OR b],'rhtlbl',2,
%                'legend',[Hl pos(1) pos(2) fs],'relegend',[row column ds],'bkxy','x OR y');
%   minortickornot: 0 no minortick; 1 one minor tick, 
%   'relegend',[-line -column dist]: sign(-line), line first; sign(-column), legend boxon
%
% EXAMPLES
%     opts = struct('figsize',[10 7.5],'tickxylblfs',[8 9],'ticksize',[0.01 0.025 1],...
%                   'xlbl',[0 0.1 1],'ylbl',[0 1 10],'axis',[0 1 0 10],'figcolor','c',...
%                   'legend',[hl 0.3 0.4],'relegend',[3 2],'lbx','1 2 3','lby','2 2 4');
%     printfig('test',opts);
%
%     printfig('e:\test','figsize',[10 7.5]);
%
%  Copyright 2007 lshuily@gmail.com
%  blank calculation codes are obtained from "expotfig.m" 
%  functions 'wysiwyg' and 'MinorTick' are used,
%  function 'changelegend' of ziliu is modified as a subfunction.
% I am gratefully acknowldge these authors for their codes.

if (nargin < 1)
  error('Too few input arguments');
end

% printfig(filename, [options,] ...)
H  = get(0,'CurrentFigure');
if isempty(H)
  error('There is no CurrentFigure');
end
Ha = get(H,'CurrentAxes'); % figure axes
% if ~LocalIsHG(H,'figure') %% H is the handle of a figure or not
%   error('First argument must be a handle to a figure.');
% end
filename = varargin{1};   % file name
if ~ischar(filename)      %% file name must be characters
  error('First argument must be a string.');
end
paramPairs = {varargin{2:end}}; %% opt = paramPairs
if nargin > 1
  if isstruct(paramPairs{1})    %% is a struct
    pcell = LocalToCell(paramPairs{1}); %% LocalToCell convert a struct to {field1,val1,field2,val2,...}
    paramPairs = {pcell{:}, paramPairs{2:end}};
  end
end
% verstr = version;
% majorver = str2num(verstr(1));
% defaults = [];
% if majorver > 5 %???
%   if ispref('exportfig','defaults') %% use default opt
%     defaults = getpref('exportfig','defaults');
%   end
% elseif exist('getappdata')
%   defaults = getappdata(0,'exportfigdefaults');
% end
% if ~isempty(defaults)
%   dcell = LocalToCell(defaults);
%   paramPairs = {dcell{:}, paramPairs{:}};
% end

% Do some validity checking on param-value pairs
if (rem(length(paramPairs),2) ~= 0)
  error(['Invalid input syntax. Optional parameters and values' ...
	 ' must be in pairs.']);
end

auto.lbrt        = [0 0 0 0];       % 0 0 0 0
auto.aps         = [0 0 0 0];
auto.figsize     = [];              % 8 6
auto.tickxylblfs = [8 9];           % 8 9
auto.ticksize    = [];              % 0.025 0.025 0
auto.xlbl        = [];
auto.ylbl        = [];
auto.axis        = [];
auto.legend      = [];
auto.relegend    = [];
auto.bkxy        = '';
auto.lbx         = '';
auto.lby         = '';
auto.rhtlbl      = 0;
auto.figtp       = 'eps';
auto.figcolor    = 'c';   % print color xxxx
opts             = auto;
% xxx here
axesLwdth = 0.5;  % set linewidth of axes

% Process param-value pairs
args = {};
for k = 1:2:length(paramPairs)
  param = lower(paramPairs{k});
  if ~ischar(param)        %% parameter name is string
    error('Optional parameter names must be strings');
  end
  value = paramPairs{k+1}; %% value of the parameter
  
  switch (param)
  case 'lbrt'      %% lbrt must be a numeric vector
    opts.lbrt = LocalToNum(value, auto.lbrt);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsScalars(opts.lbrt,4)
	error('lbrt must be a 4 numeric vector');
      end
    end
  case 'aps'      %% aps must be a numeric matrix
    opts.aps = LocalToNum(value, auto.aps);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsScalars(opts.aps(1,:),4)
	error('lbrt must be a 4 numeric vector');
      end
    end  
  case 'figsize'      %% figsize must be a numeric vector
    opts.figsize = LocalToNum(value, auto.figsize);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsPositiveScalars(opts.figsize,2) & ~LocalIsPositiveScalars(opts.figsize,3)
	error('figsize must be a 2 or 3 numeric vector > 0');
      end
    end
  case 'tickxylblfs'  %% tickxylblfs must be a numeric vector
    opts.tickxylblfs = LocalToNum(value, auto.tickxylblfs);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsPositiveScalars(opts.tickxylblfs,2)
	error('tickxylblfs must be a 2 numeric vector > 0');
      end
    end
  case 'ticksize'    %% ticksize must be a numeric vector
    opts.ticksize = LocalToNum(value,auto.ticksize);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsScalars(opts.ticksize,3)
	error('ticksize must be a 3 numeric vector');
      end
    end
  case 'xlbl'       %% xlbl must be a numeric vector
    opts.xlbl = LocalToNum(value,auto.xlbl);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsScalars(opts.xlbl,3)
	error('xlbl must be a 3 numeric vector');
      end
    end
  case 'ylbl'       %% ylbl must be a numeric vector
    opts.ylbl = LocalToNum(value,auto.ylbl);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsScalars(opts.ylbl,3)
	error('ylbl must be a 3 numeric vector');
      end
    end
  case 'rhtlbl'      %% rhtlbl must be a numeric value
    opts.rhtlbl = LocalToNum(value, auto.rhtlbl);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsScalars(opts.rhtlbl,1)
	error('lbrt must be a numeric value');
      end
    end    
  case 'axis'       %% axis must be a numeric vector
    opts.axis = LocalToNum(value,auto.axis);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsScalars(opts.axis,4)
	error('axis must be a 4 numeric vector');
      end
    end
  case 'figcolor'   %% figcolor selection
    opts.figcolor = LocalCheckAuto(lower(value),auto.figcolor);
    if ~strcmp(opts.figcolor,{'c','k','ci','ki'})
      error('figcolor must be ''c'' or ''k''.');
    end 
  case 'legend'    %% legend must be a numeric vector
    opts.legend = LocalToNum(value,auto.legend);
    Hl = opts.legend(1);
    if ~strcmpi(get(Hl,'Tag'),'legend') %% Hl is the handle of a legend or not
  error('First legend argument must be a handle to a legend.');
    end
    opts.legend = opts.legend(2:end);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsPositiveScalars(opts.relegend,2) & ~LocalIsPositiveScalars(opts.legend,3)
	error('Legend position and fontsize must be a 2 or 3 numeric scalar > 0');
      end
    end
  case 'relegend'    %% relegend must be a numeric vector
    opts.relegend = LocalToNum(value,auto.relegend);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsScalars(opts.relegend,2) & ~LocalIsScalars(opts.relegend,3)
	error('relegend must be a 2 or 3 numeric scalar');
      end
    end
  case 'figtp'    %% figtp must be characters
    opts.figtp = LocalCheckAuto(lower(value),auto.figtp);
    if ~strcmp(opts.figtp,{'eps','tiff'})
      error('figtp must be ''eps'' or ''tiff''.');
    end 
  case 'lbx'       %% lbx
    if ischar(value) & ~strcmp(value,'')   
      opts.lbx = value;
    else 
	  error('lbx must be a string');
    end
  case 'lby'       %% lby
    if ischar(value) & ~strcmp(value,'')      
      opts.lby = value;
    else 
	  error('lby must be a string');
    end 
  case 'bkxy'    %% bkxy must be characters
    opts.bkxy = LocalCheckAuto(lower(value),auto.bkxy);
    if ~strcmp(opts.bkxy,{'x','y'})
      error('bkxy must be ''x'' or ''y''');
    end 
  otherwise
    error(['Unrecognized option ' param '.']);
  end
end


%  Set figure parameters
 
 % figsize
 if ~isempty(opts.figsize)
   width = opts.figsize(1);                       
   hight = opts.figsize(2); 
   funit = get(H,'Units');
   set(H,'PaperUnits','centimeters');          % set paper units
   set(H,'Units','centimeters');               % set figure units
   ppos = [0.2*width,0.2*hight,width,hight];
   fpos = get(H,'Position');
   fpos(3:4) = [width,hight];
   set(H,'PaperPosition',ppos);                % set paper position
   set(H,'Position',fpos);                     % set figure position
   set(H,'Units',funit);
 end
 
 if ~isempty(Ha)
   set(Ha,'linewidth',axesLwdth);                          % set linewidth of axes
 % tickxylblfs   
   set(Ha,'fontsize',opts.tickxylblfs(1));                 % tick labels 8 point
   set(get(Ha,'XLabel'),'FontSize',opts.tickxylblfs(2));   % xaxis label 9 point
   set(get(Ha,'YLabel'),'FontSize',opts.tickxylblfs(2));   % yaxis label 9 point 
end

 % axis xlbl ylbl
 if ~isempty(opts.axis)
   axis(opts.axis);
 end
 if ~isempty(opts.xlbl)
   xlbls = opts.xlbl(1):opts.xlbl(2):opts.xlbl(3);  
   set(Ha,'XTick',xlbls);
   if ~strcmp(opts.lbx,'')
       xlbls = strrep(opts.lbx,' ','|');
       set(Ha,'xticklabel',xlbls); % set xtick
   else
       set(Ha,'xticklabel',sprintf(LocalTransformChar(opts.xlbl),xlbls)); 
   end
 end
 if ~isempty(opts.ylbl)
   ylbls = opts.ylbl(1):opts.ylbl(2):opts.ylbl(3);     
   set(Ha,'YTick',ylbls);
   if ~strcmp(opts.lby,'')
       ylbls = strrep(opts.lby,' ','|');
       set(Ha,'yticklabel',ylbls); % set ytick
   else
       set(Ha,'yticklabel',sprintf(LocalTransformChar(opts.ylbl),ylbls)); 
   end 
 end     
 
 % ticksize
 if ~isempty(opts.ticksize)
     set(gca,'TickLength',opts.ticksize(1:2));    % set TickLength  [0.01 0.025] 
     nTick = abs(fix(opts.ticksize(3)));
     if nTick~=0                                  % minor tick
         MinorTick(nTick) %,opts.ticksize(2));
     end
 end
 
 % legend
 if ~isempty(opts.legend)
     legend(Hl,'boxoff');
     posl = get(Hl,'position');
     posl(1) = opts.legend(1);                         
     posl(2) = opts.legend(2);                         
     set(Hl,'Position',posl);    % set legend Position
     if length(opts.legend)==3
         set(Hl,'FontSize',opts.legend(3)); 
     else
         set(Hl,'FontSize',8); 
     end
 end
 % relegend %% must at the end of other setting, as 'axis' will be affected
 if ~isempty(opts.legend) & ~isempty(opts.relegend) 
     Hftp = get(0,'CurrentFigure');
     Hatp = get(Hftp,'CurrentAxes');
     relegend(Hl,opts.relegend);
     set(0,'CurrentFigure',Hftp);
     set(Hftp,'CurrentAxes',Hatp);
 end
 
 if strcmp(opts.figtp,'eps')
      set(gca,'Color','none');      % change axis background to transparent
 end

 if ~strcmp(filename,'set')   % not set parameters
  if strcmp(opts.figtp,'eps')
     set(H,'Color','none');      % xxxchange background to transparent, not white [1,1,1]
     set(H,'InvertHardcopy','off'); % needed!
 else
     set(H,'Color',[1,1,1]);      % xxxchange background to white [1,1,1]
 end
  % correct the positions of x,y labels
  allAxes  = findall(H,'type','axes');
  for k = 1:length(allAxes)
      if strcmp(opts.figtp,'eps')
          set(allAxes(k),'Color','none');      % change axis background to transparent
      else
          set(H,'Color',[1,1,1]); 
      end
      if ~strcmpi(get(allAxes(k),'Tag'),'legend')
          ax = allAxes(k);
          set(H,'CurrentAxes',ax);    
          xylblPos(ax); % move pos of xylbl     
      end
  end 
  
  wysiwyg;                                        % what you see is what you get
  resizefig(H,opts.lbrt,opts.aps,opts.rhtlbl);    % determine figure size and crop blank area
  wysiwyg; 
  
  % break x,y-axis
  if ~isempty(opts.bkxy)
     breakxy(opts.bkxy);
  end

% print color xxxx
if strcmp(opts.figtp,'eps')      % eps
    if strcmp(opts.figcolor,'k')
        printargs = {'-deps2', '-loose', '-r300'};
    elseif strcmp(opts.figcolor,'c')
        printargs = {'-depsc2', '-loose', '-r300'};
    elseif strcmp(opts.figcolor,'ki')
        printargs = {'-deps', '-loose', '-r300'};
    elseif strcmp(opts.figcolor,'ci')
        printargs = {'-depsc', '-loose', '-r300'};
    end
elseif strcmp(opts.figtp,'tiff') % tiff
    if strcmp(opts.figcolor,'k')
        set(findobj(gca,'Type','line'),'Color','k');
    end
    printargs = {'-dtiff', '-loose', '-r300'};
end
legendTag(H,0);
print(H,filename,printargs{:});
legendTag(H,1);

% correct line type
if strcmp(opts.figtp,'eps')  % only for eps fig
    relntp(filename);
end
end
  
%
%  Local Functions
% ----------------
function cellArray = LocalGetAsCell(fig,prop,allowemptycell);
cellArray = get(fig,prop);
if nargin < 3
  allowemptycell = 0;
end
if ~iscell(cellArray) & (allowemptycell | ~isempty(cellArray))
  cellArray = {cellArray};
end
% ----------------
function bool = LocalIsPositiveScalars(value,n)
bool = isnumeric(value) & ...
       prod(size(value)) == n & value > 0;
% ----------------
function bool = LocalIsScalars(value,n)
bool = isnumeric(value) & ...
       prod(size(value)) == n;
% ----------------   
function value = LocalToNum(value,auto)
if ischar(value)
  if strcmp(value,'auto')
    value = auto;
  else
    value = str2num(value);
  end
end

% ----------------
%convert a struct to {field1,val1,field2,val2,...}
function c = LocalToCell(s)
f = fieldnames(s);
v = struct2cell(s);
opts = cell(2,length(f));
opts(1,:) = f;
opts(2,:) = v;
c = {opts{:}};

% ----------------
function c = LocalIsHG(obj,hgtype)
c = 0;
if (length(obj) == 1) & ishandle(obj) 
  c = strcmp(get(obj,'type'),hgtype);
end

% ----------------
function r = LocalUnionRect(r1,r2)
if isempty(r1)
  r = r2;
elseif isempty(r2)
  r = r1;
elseif max(r2(3:4)) > 0
  left = min(r1(1),r2(1));
  bot = min(r1(2),r2(2));
  right = max(r1(1)+r1(3),r2(1)+r2(3));
  top = max(r1(2)+r1(4),r2(2)+r2(4));
  r = [left bot right-left top-bot];
else
  r = r1;
end

% ----------------
function c = LocalLabelsMatchTicks(labs,ticks)
c = 0;
try
%   t1 = num2str(ticks(1));
%   n = length(ticks);
%   tend = num2str(ticks(n));
%   c = strncmp(labs(1),t1,length(labs(1))) & ...
%       strncmp(labs(n),tend,length(labs(n)));
  n = length(ticks);            % I am here
  tend = num2str(ticks(n));     % only compare the last one 
  c = strncmp(labs(n),tend,length(labs(n)));
end

% ----------------
function r = LocalAxesTightBoundingBox(axesR,a)
r = [];
atext = findall(a,'type','text','visible','on');
if ~isempty(atext) %% get min cover area of all text in axes(k)
  txtunits = get(atext,'units');
  set(atext,'units','points');
  res=LocalGetAsCell(atext,'extent');
  for n=1:length(atext)
    r = LocalUnionRect(r,res{n});
  end
  set(atext,'units',txtunits{1});
end
if strcmp(get(a,'visible'),'on')
  r = LocalUnionRect(r,[0 0 axesR(3:4)]);
  oldunits = get(a,'fontunits');
  set(a,'fontunits','points');
  label = text(0,0,'','parent',a,...
	       'units','points',...
	       'fontsize',get(a,'fontsize'),...
	       'fontname',get(a,'fontname'),...
	       'fontweight',get(a,'fontweight'),...
	       'fontangle',get(a,'fontangle'),...
	       'visible','off');
  fs = get(a,'fontsize');

  % handle y axis tick labels
  ry = [0 -fs/2 0 axesR(4)+fs]; %% bottom 1/2*fs up 1/2*fs
  ylabs = get(a,'yticklabel');
  yticks = get(a,'ytick');
  maxw = 0;
  if ~isempty(ylabs)
    for n=1:size(ylabs,1)
      set(label,'string',ylabs(n,:));
      ext = get(label,'extent');
      maxw = max(maxw,ext(3)); %% get the max ylabs width
    end
%     if ~LocalLabelsMatchTicks(ylabs,yticks) & ... %% if labs & ticks not matches
% 	  strcmp(get(a,'xaxislocation'),'bottom')     %% xaxis locats at bottom
%       ry(4) = ry(4) + 1.5*ext(4);                 %% there is exponent at up, add 1.5*ylabel(height)
%     end
    if strcmp(get(a,'yaxislocation'),'left')
      ry(1) = -(maxw+5);                          %% left, maxwidth+5
    else
      ry(1) = axesR(3);                           %% right, no space at right
    end
    ry(3) = maxw+5;
    r = LocalUnionRect(r,ry);
  end

  % handle x axis tick labels
  rx = [0 0 0 fs+5];
  xlabs = get(a,'xticklabel');
  xticks = get(a,'xtick');
  if ~isempty(xlabs)
    if strcmp(get(a,'xaxislocation'),'bottom')
      rx(2) = -(fs+5);          %% bottom, fs+5
      if ~LocalLabelsMatchTicks(xlabs,xticks);
	rx(4) = rx(4) + 2*fs;
	rx(2) = rx(2) - 2*fs;
      end
    else
      rx(2) = axesR(4);
      % exponent is still below axes
%       if ~LocalLabelsMatchTicks(xlabs,xticks);
% 	rx(4) = rx(4) + axesR(4) + 2*fs;
% 	rx(2) = -2*fs;
%       end
    end
    set(label,'string',xlabs(1,:));
    ext1 = get(label,'extent');
    rx(1) = -ext1(3)/2;
    set(label,'string',xlabs(size(xlabs,1),:));
    ext2 = get(label,'extent');
    rx(3) = axesR(3) + (ext2(3) + ext1(3))/2;
    r = LocalUnionRect(r,rx);
  end
  set(a,'fontunits',oldunits);
  delete(label);
end

% ----------------
function val = LocalCheckAuto(val, auto)
if ischar(val) & strcmp(val,'auto')
  val = auto;
end

% ----------------
function strings = LocalTransformChar(vector)
ninteger = 0;
npoint   = 0;
for i=vector(1):vector(2):vector(3)
    ninteger = max(ninteger,length(num2str(round(abs(i)))));
    ntemp = 0;
    if i<0.1&i>-0.1&i~=0
        ntemp = ntemp+1;
        i = 10*i;
    end
    while round(10*i)/10~=round(i)
        ntemp = ntemp+1;
        i = 10*i;
    end
    npoint = max(npoint,ntemp);
end
strings = strcat('%',num2str(ninteger),'.',num2str(npoint),'f|');

% ----------------
function  resizefig(H,lbrt,aps,rhtlbl) % change figure size and blank area
%% determine figure bounds to surround axes
allAxes  = findall(H,'type','axes');
set(allAxes,'units','points');
apos     = LocalGetAsCell(allAxes,'Position');
oldunits = get(H,'Units');
set(H,'units','points');
fpos     = get(H,'Position');
% calculate (1) lmin,bmin; (2) lblank,bblank,rblank,tblank; (3) RatoX,RatoY
lmin = fpos(3);
bmin = fpos(4);
rmax = 0;
tmax = 0;
fr = [];
for k=1:length(allAxes)
    if ~strcmpi(get(allAxes(k),'Tag'),'legend')
        axesR = apos{k};
        lmin = min(lmin,axesR(1));
        bmin = min(bmin,axesR(2));
        rmax = max(rmax,axesR(1)+axesR(3));
        tmax = max(tmax,axesR(2)+axesR(4));
        r = LocalAxesTightBoundingBox(axesR, allAxes(k));
        r(1:2) = r(1:2) + axesR(1:2);
        fr = LocalUnionRect(fr,r); 
    end
end
if isempty(fr)
    lmin = 0;
    bmin = 0;
    rmax = fpos(3);
    tmax = fpos(4);
    fr = [0 0 fpos(3:4)];
end
fr(1) = fr(1)+lbrt(1)+4;     % I am here xxxx 4 2 3 0; 3 -1.5 -1 2
fr(2) = fr(2)+lbrt(2)-0;
fr(3) = fr(3)-lbrt(3)-lbrt(1)-2.4;     
fr(4) = fr(4)-lbrt(4)-lbrt(2)-0.5;
lbk   = lmin-fr(1);        % left blank
bbk   = bmin-fr(2);
rbk   = fr(1)+fr(3)-rmax; 
tbk   = fr(2)+fr(4)-tmax;
ratioX = (fpos(3)-lbk-rbk)/(fr(3)-lbk-rbk);
ratioY = (fpos(4)-bbk-tbk)/(fr(4)-bbk-tbk);
% change axes position & size
j = 0;
for k = 1:length(allAxes)
    j = j+1;
    ax = allAxes(k);
    if aps(1,4)>0 % set axes position
        if ~strcmpi(get(allAxes(k),'Tag'),'legend')
            xyp = xylblGetPos(ax);
            oldunits = get(ax,'Units');
            set(ax,'units','normalized');
            if length(allAxes) < size(aps,1)
                set(ax,'Position',aps(j+1,:));  % xxxx
            else
                set(ax,'Position',aps(j,:));
            end
            set(ax,'units',oldunits);
            xylblRePos(ax,xyp,rhtlbl);
        else
            j = j-1;      
            r = apos{k};
            r(1) = (r(1)-lmin)*ratioX+lbk;
            r(2) = (r(2)-bmin)*ratioY+bbk;
            set(ax,'Position',r);
        end
    else      % calculate axes position
        r = apos{k};
        r(1) = (r(1)-lmin)*ratioX+lbk;
        r(2) = (r(2)-bmin)*ratioY+bbk;
        if ~strcmpi(get(allAxes(k),'Tag'),'legend')
            r(3) = r(3)*ratioX;
            r(4) = r(4)*ratioY;
            xyp = xylblGetPos(ax);
            set(ax,'Position',r);
            xylblRePos(ax,xyp,rhtlbl);
        else
            set(ax,'Position',r);
        end
    end
end
set(H, 'Units',oldunits);
set(allAxes,'units','normalized');

% ----------------
function wysiwyg
%       Dan(K) Braithwaite, Dept. of Hydrology U.of.A  11/93
unis = get(gcf,'units');
ppos = get(gcf,'paperposition');
set(gcf,'units',get(gcf,'paperunits'));
pos = get(gcf,'position');
pos(3:4) = ppos(3:4);
set(gcf,'position',pos);
set(gcf,'units',unis);

% ----------------
function MinorTick(N) %,exy)
% N is the number of MinorTicks
if nargin < 1
    return;
end
N = N+1;
xx = max(xlim);
xn = min(xlim);
yx = max(ylim);
yn = min(ylim);
xt = get(gca,'xtick');
yt = get(gca,'ytick');
if length(xt) < 3 | length(yt) < 3
    return;
elseif strcmpi(get(gca,'Tag'),'legend')
    return;
else
    if abs(xt(1)+xt(3)-2*xt(2)) > 1.0E-16  % xt(3)/xt(2)==10 & xt(2)/xt(1)==10
        sgnx = 1;
    else
        sgnx = 0;
    end
    if abs(yt(1)+yt(3)-2*yt(2)) > 1.0E-16 % yt(3)/yt(2)==10 & yt(2)/yt(1)==10
        sgny = 1;
    else
        sgny = 0;
    end
end
L  = strcmp(get(gca,'box'),'on');
dx = (xt(2)-xt(1))/N;
dy = (yt(2)-yt(1))/N;
% if exy > 0
%     dx = exy;
% elseif exy < 0
%     dy = abs(exy);
% end
xl = fliplr(xt(1)-dx:-dx:xn);
xh = xt(end)+dx:dx:xx;
xs = xt(1):dx:xt(end);
xs = setdiff([xl xs xh],xt);
yl = fliplr(yt(1)-dy:-dy:yn);
yh = yt(end)+dy:dy:yx;
ys = yt(1):dy:yt(end);
ys = setdiff([yl ys yh],yt);

tlen = get(gca,'Ticklength');
aps  = get(gca,'Position');
fps  = get(gcf,'Position');
wth  = aps(3)*fps(3);
hgt  = aps(4)*fps(4);
if wth > hgt
    fx = wth/hgt;
    fy = 1;
else
    fx = 1;
    fy = hgt/wth;
end
if sgny==1 %%| exy < 0  % x tick
    Lx = tlen(1)*log10(yx/yn)*fx*0.5;
    Lxl= yn*10^Lx;
    Lxh= yx/10^Lx;
else
    Lx = tlen(1)*(yx-yn)*fx*0.5;
    Lxl= yn+Lx;
    Lxh= yx-Lx;
end
if sgnx==1 %%| exy > 0 % y tick
    Ly = tlen(1)*log10(xx/xn)*fy*0.5;
    Lyl= xn*10^Ly;
    Lyh= xx/10^Ly;
else
    Ly = tlen(1)*(xx-xn)*fy*0.5;
    Lyl= xn+Ly;
    Lyh= xx-Ly;
end
axesLwdth = get(gca,'linewidth');  
hold on;
if ~strcmpi(get(gca,'Tag'),'y') %% only plot x
    if sgnx==0                    %% not exp
        if strcmp(get(gca,'xaxislocation'),'top')
            Lxll = [yx,Lxh];
            Lxhh = [yn,Lxl];
        else
            Lxll = [yn,Lxl];
            Lxhh = [yx,Lxh];
        end
        for f = xs;
            plot([f,f],Lxll,'k','linewidth',axesLwdth);
            if L==1;plot([f,f],Lxhh,'k','linewidth',axesLwdth);end; 
        end
    end
end
if ~strcmpi(get(gca,'Tag'),'x') %% only plot y
    if sgny==0                    %% not exp
        if strcmp(get(gca,'yaxislocation'),'right')
            Lyll = [xx,Lyh];
            Lyhh = [xn,Lyl];
        else
            Lyll = [xn,Lyl];
            Lyhh = [xx,Lyh];
        end
        for f = ys;
            plot(Lyll,[f,f],'k','linewidth',axesLwdth);
            if L==1;plot(Lyhh,[f,f],'k','linewidth',axesLwdth);end;
        end
    end
end
hold off;

% ----------------
function breakxy(sxy)
r  = 2;
d0 = 0.006;
% get axes & pos
allAxes = findall(gcf,'type','axes');
ax = allAxes;
j  = 0;
for k = 1:length(allAxes)
    if strcmpi(get(allAxes(k),'Tag'),sxy)
        j = j+1;
        ax(j) = allAxes(k);
    elseif strcmpi(get(allAxes(k),'Tag'),setdiff({'x','y'},sxy))
        tHand = allAxes(k);
    end
end
ax = ax(1:j);
if j<2 | j>4
    error('Axes number error');
else
    ps = zeros(j,4);
    for k = 1:j
        ps(k,:) = get(ax(k),'position');
    end
    if strcmpi(sxy,'x')
        tp = ps(:,1);
    else
        tp = ps(:,2);
    end
    allAxes = ax;   % sort
    mxv = max(tp);
    axps = ps;
    for k = 1:j
        [mval mdex] = min(tp);
        allAxes(k) = ax(mdex);
        axps(k,:) = ps(mdex,:);
        tp(mdex) = 2*mxv;
    end
end
% width & hight
fps = get(gcf,'position');
if strcmpi(sxy,'x')
    width = sum(axps(:,3));
    hight = axps(1,4);
else
    width = axps(1,3);
    hight = sum(axps(:,4));
end
if width*fps(3) > hight*fps(4)
    dx = d0;
    dy = d0*width*fps(3)/(hight*fps(4));
else
    dx = d0*hight*fps(4)/(width*fps(3));
    dy = d0;
end
width = width+axps(1,1);
hight = hight+axps(1,2);
if strcmp(sxy,'x')
    axh = axes('position',[0 0 width 1]);
    set(axh,'color','none','xcolor','w','ycolor','w','box','off');
    axis([0 1 0 1])
    % x break
    y1 = axps(1,2);
    y2 = axps(1,2)+axps(1,4);
    for i = 1:length(allAxes)-1
        LWth = get(allAxes(i),'linewidth');
        Lbox = strcmp(get(allAxes(i),'box'),'on');
        xi = (sum(axps(1:i,3))+axps(1,1))/width;
        xl = [xi-dx-r*dx xi-dx+r*dx];
        xh = [xi+dx-r*dx xi+dx+r*dx];
        yl = [y1-r*dy y1+r*dy];
        yh = [y2-r*dy y2+r*dy];
        hold on
        plot([xi-dx xi+dx],[y1 y1],'-w','linewidth',4*LWth);
        plot(xl,yl,'-k','linewidth',LWth);
        plot(xh,yl,'-k','linewidth',LWth);
        if Lbox
            plot([xi-dx xi+dx],[y2 y2],'-w','linewidth',4*LWth);
            plot(xl,yh,'-k','linewidth',LWth);
            plot(xh,yh,'-k','linewidth',LWth);
        end
        hold off
    end
    set(axh,'xtick',[]);
    set(axh,'ytick',[]);
    set(tHand,'xtick',[]);
else
    axh = axes('position',[0 0 1 hight]);
    set(axh,'color','none','xcolor','w','ycolor','w','box','off');
    axis([0 1 0 1])
    % y break
    x1 = axps(1,1);
    x2 = axps(1,1)+axps(1,3);
    for i = 1:length(allAxes)-1
        LWth = get(allAxes(i),'linewidth');
        Lbox = strcmp(get(allAxes(i),'box'),'on');
        yi = (sum(ps(1:i,4))+axps(1,2))/hight;
        xl = [x1-r*dx x1+r*dx];
        xh = [x2-r*dx x2+r*dx];
        yl = [yi-dy-r*dy yi-dy+r*dy];
        yh = [yi+dy-r*dy yi+dy+r*dy];
        hold on
        plot([x1 x1],[yi-dy yi+dy],'-w','linewidth',4*LWth);
        plot(xl,yl,'-k','linewidth',LWth);
        plot(xl,yh,'-k','linewidth',LWth);
        if Lbox
            plot([x2 x2],[yi-dy yi+dy],'-w','linewidth',4*LWth);
            plot(xh,yl,'-k','linewidth',LWth);
            plot(xh,yh,'-k','linewidth',LWth);
        end
        hold off
    end
    set(axh,'xtick',[]);
    set(axh,'ytick',[]);
    set(tHand,'ytick',[]);
end

% ----------------
function relegend(Hl,matx)
% change the form of the legend
ml = abs(matx(1));
nl = abs(matx(2));
if length(matx) < 3
    ds = 0.2;
else
    ds = matx(3);
end
[legh,objh,outh,outm] = legend(Hl);
suml  = (length(objh)-1)/2;
if suml > 1
    Hline = findall(Hl,'type','line');
    xdat1 = get(Hline(2),'xdata');
    ydat1 = get(Hline(1),'ydata');
    ydat2 = get(Hline(3),'ydata');
    deltm = ydat2-ydat1;
    spacm = 1-deltm*length(Hline)/2;   % border text distance
    newdeltm = deltm/(spacm+ml*deltm);
    newspacm = spacm/(spacm+ml*deltm);
    boxsz = get(objh(1),'extent');     % text position
                  
    oldlpos = get(Hl,'Position');
    sH      = get(Hl,'Children');
    newlpos = [oldlpos(1:2),1.0*nl*oldlpos(3),0.9*(ml*deltm+spacm)*oldlpos(4)]; % lgdis 1.0(hspace) 0.85(vspace)
    axh = axes('Position',newlpos);
    hold on
    axis([0,1,0,1]);
    set(gca,'xtick',[],'ytick',[]);
    box on;
for m = 1:ml
    for n = 1:nl
          if sign(matx(1))==-1 % line first
              k = (m-1)*nl+n;
          else
              k = (n-1)*ml+m;
          end
        if k<=suml
         vtp = ds*(n-1);      % distance between columns 
         xd1 = xdat1-xdat1(1);
         llg = 1.26;          % length of legend line for 8font
         gap = 0.70;          % ?
         tlx = llg*(xd1/nl+(n-1)/nl-vtp)+gap*xdat1(1);
         llx = llg*([mean(xd1)/nl mean(xd1)/nl]+(n-1)/nl-vtp)+gap*xdat1(1);
         ltx = llg*((boxsz(1)-xdat1(1))/nl+(n-1)/nl-vtp)+gap*xdat1(1)-0.02/nl; % distance of symbol & text
         posy = 1-(newspacm+newdeltm)/2-(m-1)*newdeltm;
         tline(m,n) = plot(tlx,[1,1]*posy);
         lline(m,n) = plot(llx,[1,1]*posy,'color',[1,1,1]);
         ltext(m,n) = text(ltx,posy,char(outm(k)));
         set(lline(m,n),'marker',get(sH(2*suml+1-2*k),'marker'),'markersize',get(sH(2*suml+1-2*k),'markersize'),...
             'MarkerFaceColor',get(sH(2*suml+1-2*k),'MarkerFaceColor')); % marker face
         set(lline(m,n),'color', get(sH(2*suml+2-2*k),'color'), 'LineStyle', 'none');
         set(tline(m,n),'color', get(sH(2*suml+2-2*k),'color'),...
             'LineStyle', get(sH(2*suml+2-2*k),'Linestyle'),...
             'LineWidth',get(sH(2*suml+2-2*k),'LineWidth'));
         set(ltext(m,n),'fontsize',get(sH(2*suml+1),'fontsize'));
        end
    end
end
   hold off
   delete(Hl);
   set(axh,'Tag','legend');
   if sign(matx(2))==1 % box off
       set(axh,'visible','off');
   end
end

% ----------------
function relntp(filename)
% change line type
% stra = '/SO { [] 0 setdash } bdef';
% str1 = '/DO { [ 1 dpi2point mul 1.5 dpi2point mul] 0 setdash} bdef';
% str2 = '/DA { [ 4 dpi2point mul 1.5 dpi2point mul] 0 setdash} bdef';
% str3 = '/DD { [ 1 dpi2point mul 2 dpi2point mul 4 dpi2point mul 2';
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

% ----------------
function xylblPos(a) % correct the positions of x,y labels
dxy = [2 4.5]; % 2 5 distance of xylbl  zzzz I am here
xyp = [0 0];
% x label
if strcmp(get(a,'xaxislocation'),'bottom')
    xyp(1) =  dxy(1);     %% bottom, ht+5+5
else
    xyp(1) = -dxy(1);    
end
% y label
if strcmp(get(a,'yaxislocation'),'left')
    xyp(2) =  dxy(2);   %% left, wd+5+5
else
    xyp(2) = -dxy(2);    
end
% set
hxlbl = get(a,'xlabel');
hylbl = get(a,'ylabel');
xylunits = get(hxlbl,'units');
% xlbl
xlblf = get(hxlbl,'fontsize');
xlblv = get(hxlbl,'string');
if ~isempty(xlblv)
    delete(hxlbl);
    hxlbl = XLabel(xlblv);
    set(hxlbl,'fontsize',xlblf);
    set(hxlbl,'units','points');
    posxl = get(hxlbl,'position');
    posxl(2) = posxl(2)+xyp(1);
    set(hxlbl,'position',posxl);
end
% ylbl
ylblf = get(hylbl,'fontsize');
ylblv = get(hylbl,'string');
if ~isempty(ylblv)
    delete(hylbl);
    hylbl = YLabel(ylblv);
    set(hylbl,'fontsize',ylblf);
    set(hylbl,'units','points');
    posyl = get(hylbl,'position');
    posyl(1) = posyl(1)+xyp(2);
    set(hylbl,'position',posyl);
end
set([hxlbl hylbl],'units',xylunits);

% ----------------
function xyp = xylblGetPos(a) % get the positions of x,y labels
xyp = [0 0];
% set
hxlbl = get(a,'xlabel');
hylbl = get(a,'ylabel');
xylunits = get(hxlbl,'units');
set([hxlbl hylbl],'units','points');
% xlbl
xlblv = get(hxlbl,'string');
if ~isempty(xlblv)
    posxl = get(hxlbl,'position');
    xyp(1) = posxl(2);
end
% ylbl
ylblv = get(hylbl,'string');
if ~isempty(ylblv)
    posyl = get(hylbl,'position');
    xyp(2) = posyl(1);
end
set([hxlbl hylbl],'units',xylunits);

% ----------------
function xylblRePos(a,xyp,rhtlbl) % correct the positions of x,y labels again
set(gcf,'CurrentAxes',a);
% set
hxlbl = get(a,'xlabel');
hylbl = get(a,'ylabel');
xylunits = get(hxlbl,'units');
% xlbl
xlblf = get(hxlbl,'fontsize');
xlblv = get(hxlbl,'string');
if ~isempty(xlblv)
    delete(hxlbl);
    hxlbl = XLabel(xlblv);
    set(hxlbl,'fontsize',xlblf);
    set(hxlbl,'units','points');
    posxl = get(hxlbl,'position');
    posxl(2) = xyp(1);
    set(hxlbl,'position',posxl);
end
% ylbl
ylblf = get(hylbl,'fontsize');
ylblv = get(hylbl,'string');
if ~isempty(ylblv)
    delete(hylbl);
    hylbl = YLabel(ylblv);
    set(hylbl,'fontsize',ylblf);
    set(hylbl,'units','points');
    posyl = get(hylbl,'position');
    if strcmp(get(a,'yaxislocation'),'left')
        posyl(1) = xyp(2);
    else % tmp correct!
        oldu = get(gcf,'units');
        set(gcf,'units','points');
        fps = get(gcf,'position');
        posyl(1) = fps(3)-42.3+rhtlbl;    % xxxx rtdis 39.5
        set(gcf,'units',oldu);
    end
    set(hylbl,'position',posyl);
end
set([hxlbl hylbl],'units',xylunits);

% ----------------
function legendTag(H,tags) % set tag of legend
if tags==0
    allAxes  = findall(H,'type','axes');
    for k=1:length(allAxes)
        if strcmpi(get(allAxes(k),'Tag'),'legend')
            set(allAxes(k),'Tag','notleg');
            break;
        end
    end
else
    allAxes  = findall(H,'type','axes');
    for k=1:length(allAxes)
        if strcmpi(get(allAxes(k),'Tag'),'notleg')
            set(allAxes(k),'Tag','legend');
            break;
        end
    end
end