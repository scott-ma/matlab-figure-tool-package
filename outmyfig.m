function varargout = printfig(varargin)
% output figure only with the gca reserved
%
%   EXPORTFIG(H, FILENAME, OPT) writes the figure H to FILENAME.  H is
%   a figure handle and FILENAME is a string that specifies the
%   name of the output file.  OPT is setting parameters.
%
%  opts = struct('lbrt',[3 1 0 0],'figsize',[width height],'tickxylblfs',[ticklblfs xlblfs ylblfs],
%                'ticksize',[xfactor yfactor minortickornot],'xlbl',[x1 xstep xn],
%                'ylbl',[y1 ystep yn],'axis',[x1 xn y1 yn],figcolor',[c OR b],
%                'legend',[Hl pos(1) pos(2) fs],'relegend',[row column]);

%if (nargin < 2)
if (nargin < 1)
  error('Too few input arguments');
end

%H  = varargin{1};          % figure hand
H  = get(0,'CurrentFigure');
Ha = get(H,'CurrentAxes'); % figure axes
if ~LocalIsHG(H,'figure') %% H is the handle of a figure or not
  error('First argument must be a handle to a figure.');
end
filename = varargin{1};    % file name
if ~ischar(filename)      %% file name must be characters
  error('Second argument must be a string.');
end
paramPairs = {varargin{2:end}}; %% opt = paramPairs
if nargin > 1
  if isstruct(paramPairs{1})    %% is a struct
    pcell = LocalToCell(paramPairs{1}); %% LocalToCell convert a struct to {field1,val1,field2,val2,...}
    paramPairs = {pcell{:}, paramPairs{2:end}};
  end
end
verstr = version;
majorver = str2num(verstr(1));
defaults = [];
if majorver > 5 % version of Matlab
  if ispref('exportfig','defaults') %% use default opt
    defaults = getpref('exportfig','defaults');
  end
elseif exist('getappdata')
  defaults = getappdata(0,'exportfigdefaults');
end
if ~isempty(defaults)
  dcell = LocalToCell(defaults);
  paramPairs = {dcell{:}, paramPairs{:}};
end

% Do some validity checking on param-value pairs
if (rem(length(paramPairs),2) ~= 0)
  error(['Invalid input syntax. Optional parameters and values' ...
	 ' must be in pairs.']);
end

auto.lbrt        = [0 0 0 0];
auto.figsize     = [10 7.5];
auto.tickxylblfs = [7 8 8];
auto.ticksize    = [0.01 0.025 0];
auto.xlbl        = [];
auto.ylbl        = [];
auto.axis        = [];
auto.figcolor    = 'c';
auto.legend      = [];
auto.relegend    = [];
auto.figtp       = 'eps';
opts             = auto;

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
      if ~LocalIs4Scalars(opts.lbrt)
	error('Width must be a 4 integer vector');
      end
    end   
   case 'figsize'      %% figsize must be a numeric vector
    opts.figsize = LocalToNum(value, auto.figsize);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsPositiveScalar2(opts.figsize)
	error('Width must be a 2 numeric vector> 0');
      end
    end
   case 'tickxylblfs'  %% tickxylblfs must be a numeric vector
    opts.tickxylblfs = LocalToNum(value, auto.tickxylblfs);
    if ~ischar(value) | ~strcmp(value,'auto')
      if(~LocalIsPositiveScalars(opts.tickxylblfs))
	error('Tickxylblfs must be a 3 numeric vector > 0');
      end
    end
   case 'ticksize'    %% ticksize must be a numeric vector
    opts.ticksize = LocalToNum(value,auto.ticksize);
    opts.ticksize(3) = opts.ticksize(3)+1;
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsPositiveScalars(opts.ticksize)
	error('Ticksize must be a 3 numeric scalar >= 0');
      end
    end
    opts.ticksize(3) = opts.ticksize(3)-1;
    case 'xlbl'       %% xlbl must be a numeric vector
    opts.xlbl = LocalToNum(value,auto.xlbl);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIs3Scalars(opts.xlbl)
	error('Xlbl must be a 3 numeric vector');
      end
    end
    case 'ylbl'       %% ylbl must be a numeric vector
    opts.ylbl = LocalToNum(value,auto.ylbl);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIs3Scalars(opts.ylbl)
	error('Xlbl must be a 3 numeric vector');
      end
    end
    case 'axis'       %% axis must be a numeric vector
    opts.axis = LocalToNum(value,auto.axis);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIs4Scalars(opts.axis)
	error('Xlbl must be a 4 numeric vector');
      end
    end
    case 'figcolor'   %% figcolor selection
    opts.figcolor = LocalCheckAuto(lower(value),auto.figcolor);
    if ~strcmp(opts.figcolor,{'c','k'})
      error('Color must be ''c'' or ''k''.');
    end 
    case 'legend'    %% legend must be a numeric vector
    opts.legend = LocalToNum(value,auto.legend);
    Hl = opts.legend(1);
    if ~strcmpi(get(Hl,'Tag'),'legend') %% Hl is the handle of a legend or not
  error('First legend argument must be a handle to a legend.');
    end
    opts.legend = opts.legend(2:end);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsPositiveScalars(opts.legend)
	error('Legend position and fontsize must be a 3 numeric scalar > 0');
      end
    end
    case 'relegend'    %% relegend must be a numeric vector
    opts.relegend = LocalToNum(value,auto.relegend);
    if ~ischar(value) | ~strcmp(value,'auto')
      if ~LocalIsPositiveScalar2(opts.relegend)
	error('Relegend must be a 2 numeric scalar > 0');
      end
    end
    case 'figtp'    %% figtp must be characters
    opts.figtp = LocalCheckAuto(lower(value),auto.figtp);
    if ~strcmp(opts.figtp,{'eps','tiff'})
      error('Figure type must be ''eps'' or ''tiff''.');
    end 
   otherwise
    error(['Unrecognized option ' param '.']);
  end
end


%  Set figure parameters
  
 % figsize
 width = opts.figsize(1);                       
 hight = opts.figsize(2);                   
 
 set(H,'Color',[1,1,1])                      % change background to white
 set(Ha,'Units','centimeters'); 
 set(Ha,'Position',[0,0,width,hight]); 
 set(Ha,'Units','point');
 apos = get(Ha,'Position');
 pl = apos(1)+opts.lbrt(1);
 pw = apos(2)+opts.lbrt(2);
 pd = apos(3)+opts.lbrt(3);
 ph = apos(4)+opts.lbrt(4);
 set(Ha,'Position',[pl,pw,pd,ph]);  % ==h
 set(H,'PaperUnits','centimeters'); 
 set(H, 'PaperPosition', [0,0,width,hight]);

 % tickxylblfs   
   set(Ha,'fontsize',opts.tickxylblfs(1));                 % tick labels 8 point
   set(get(Ha,'XLabel'),'FontSize',opts.tickxylblfs(2));   % xaxis label 9 point
   set(get(Ha,'YLabel'),'FontSize',opts.tickxylblfs(3));   % yaxis label 9 point 

 % axis xlbl ylbl
 if ~isempty(opts.axis)
   axis(opts.axis);
 end
 if ~isempty(opts.xlbl)
   xlbls = opts.xlbl(1):opts.xlbl(2):opts.xlbl(3);  
   set(Ha,'XTick',xlbls);
   set(Ha,'xticklabel',sprintf(LocalTransformChar(opts.xlbl),xlbls));  
 end
 if ~isempty(opts.ylbl)
   ylbls = opts.ylbl(1):opts.ylbl(2):opts.ylbl(3);     
   set(Ha,'YTick',ylbls);
   set(Ha,'yticklabel',sprintf(LocalTransformChar(opts.ylbl),ylbls));  
 end     
 
 % ticksize
 set(gca,'TickLength',opts.ticksize(1:2));    % set TickLength  [0.01 0.025] 
 if opts.ticksize(3)~=0                       % minor tick
     MinorTick(abs(fix(opts.ticksize(3))));
 end
 
 % legend
 if ~isempty(opts.legend)
     legend(Hl,'boxoff');
     posl = get(Hl,'position');
     posl(1) = opts.legend(1);                         
     posl(2) = opts.legend(2);                         
     set(Hl,'Position',posl);    % set legend Position
     set(Hl,'FontSize',opts.legend(3)); 
 end
 % relegend %% must at the end of other setting, as 'axis' will be affected
 if ~isempty(opts.legend) & ~isempty(opts.relegend) 
     Hftp = get(0,'CurrentFigure');
     Hatp = get(Hftp,'CurrentAxes');
     relegend(Hl,opts.relegend);
     set(0,'CurrentFigure',Hftp);
     set(Hftp,'CurrentAxes',Hatp);
 end

wysiwyg;                      % what you see is what you get

% correct the positions of x,y labels
xlblv = get(get(Ha,'xlabel'),'string');
ylblv = get(get(Ha,'ylabel'),'string');
delete(get(Ha,'xlabel'));
delete(get(Ha,'ylabel'));
XLabel(xlblv);
YLabel(ylblv);
if ~isempty(opts.tickxylblfs)
  set(get(Ha,'XLabel'),'FontSize',opts.tickxylblfs(2));   % xaxis label 9 point
  set(get(Ha,'YLabel'),'FontSize',opts.tickxylblfs(3));   % yaxis label 9 point 
end

% print color xxxx
opts.figtp = 'eps';  % tmp use
opts.figcolor = 'k'; % tmp use
if strcmp(opts.figtp,'eps') % eps
    if strcmp(opts.figcolor,'c')
        printargs = {'-depsc2', '-loose', '-r300'};
    else
        printargs = {'-deps', '-loose', '-r300'};
    end
elseif strcmp(opts.figtp,'tiff') % tiff
    printargs = {'-dtiff', '-loose', '-r300'};
end
print(H,filename,printargs{:});


%
%  Local Functions
%
function cellArray = LocalGetAsCell(fig,prop,allowemptycell);
cellArray = get(fig,prop);
if nargin < 3
  allowemptycell = 0;
end
if ~iscell(cellArray) & (allowemptycell | ~isempty(cellArray))
  cellArray = {cellArray};
end

function bool = LocalIsPositiveScalar2(value)
bool = isnumeric(value) & ...
       prod(size(value)) == 2 & ...
       value > 0;
   
function bool = LocalIsPositiveScalars(value)
bool = isnumeric(value) & ...
       prod(size(value)) == 3 & ...
       value > 0;

function bool = LocalIs3Scalars(value)
bool = isnumeric(value) & ...
       prod(size(value)) == 3;

function bool = LocalIs4Scalars(value)
bool = isnumeric(value) & ...
       prod(size(value)) == 4;
   
function value = LocalToNum(value,auto)
if ischar(value)
  if strcmp(value,'auto')
    value = auto;
  else
    value = str2num(value);
  end
end

%convert a struct to {field1,val1,field2,val2,...}
function c = LocalToCell(s)
f = fieldnames(s);
v = struct2cell(s);
opts = cell(2,length(f));
opts(1,:) = f;
opts(2,:) = v;
c = {opts{:}};

function c = LocalIsHG(obj,hgtype)
c = 0;
if (length(obj) == 1) & ishandle(obj) 
  c = strcmp(get(obj,'type'),hgtype);
end

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

function r = LocalAxesTightBoundingBox(axesR, a)
r = [];
atext = findall(a,'type','text','visible','on');
if ~isempty(atext) %% get min cover area of all text in axes(k)
  set(atext,'units','points');
  res=LocalGetAsCell(atext,'extent');
  for n=1:length(atext)
    r = LocalUnionRect(r,res{n});
  end
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
  ry = [0 -fs/2 0 axesR(4)+fs]; %% bottom 1/2*fs up fs
  ylabs = get(a,'yticklabel');
  yticks = get(a,'ytick');
  maxw = 0;
  if ~isempty(ylabs)
    for n=1:size(ylabs,1)
      set(label,'string',ylabs(n,:));
      ext = get(label,'extent');
      maxw = max(maxw,ext(3)); %% get the max ylabs width
    end
    if ~LocalLabelsMatchTicks(ylabs,yticks) & ... %% if labs & ticks not matches
	  strcmp(get(a,'xaxislocation'),'bottom')     %% xaxis locats at bottom
      ry(4) = ry(4) + 1.5*ext(4);                 %% there is exponent at up, add 1.5*ylabel(height)
    end
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
      if ~LocalLabelsMatchTicks(xlabs,xticks);
	rx(4) = rx(4) + axesR(4) + 2*fs;
	rx(2) = -2*fs;
      end
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

function val = LocalCheckAuto(val, auto)
if ischar(val) & strcmp(val,'auto')
  val = auto;
end

function strings = LocalTransformChar(vector)
ninteger = 0;
npoint   = 0;
for i=vector(1):vector(2):vector(3)
    ninteger = max(ninteger,length(num2str(fix(abs(i)))));
    ntemp = 0;
    while fix(10*i)/10~=fix(i)
        ntemp = ntemp+1;
        i = 10*i;
    end
    npoint = max(npoint,ntemp);
end
strings = strcat('%',num2str(ninteger),'.',num2str(npoint),'f|');

function  resizefig(H,Ha,figsize) % change figure size and blank area
width  = figsize(1);
height = figsize(2);
  % adjust figure bounds to surround axes
      oldapos   = get(Ha,'Position');
      allAxes   = findall(H, 'type', 'axes');
      set(allAxes,'units','points');
      apos      = LocalGetAsCell(allAxes,'Position');
      oldunits  = get(H,'Units');
      set(H,'units','points');
      oldfpos   = get(H,'Position');
      fr = [];
      for k=1:length(allAxes)
	if ~strcmpi(get(allAxes(k),'Tag'),'legend')
	  axesR = apos{k};
	  r = LocalAxesTightBoundingBox(axesR, allAxes(k));
	  r(1:2) = r(1:2) + axesR(1:2);
	  fr = LocalUnionRect(fr,r); 
	end
      end
      if isempty(fr)
	fr = [0 0 oldfpos(3:4)];
      end
  fr(1)   = fr(1)+1; % minor change            == I am here ==
  fr(2)   = fr(2)+2;
%   fr(3)   = fr(3)-2;
      for k=1:length(allAxes)
	ax = allAxes(k);
	r = apos{k};
	r(1:2) = r(1:2) - fr(1:2);
	set(ax,'Position',r);
      end
      r = [oldfpos(1) oldfpos(2)+oldfpos(4)-fr(4) fr(3:4)];
      set(H,'Position',r);
      set(H,'units','centimeters')
      tmpfpos = get(H,'Position');
      
%       fr(1)   = fr(1)+3; % minor change            == I am here ==
      aleft   = oldapos(1)-fr(1)/fr(3)*tmpfpos(3); % left blank
      abottom = oldapos(2)-fr(2)/fr(3)*tmpfpos(3); % bottom blank
      aright  = tmpfpos(3)-oldapos(3)-aleft;       % right blank
      aup     = tmpfpos(4)-oldapos(4)-abottom;     % up blank
      awidth  = width-aleft-aright;       % axes width
      aheight = height-abottom-aup;       % axes height
      newapos = [aleft abottom awidth aheight];

      newfpos = [tmpfpos(1:2) 1.25*oldapos(3) 1.25*oldapos(4)];
      set(H, 'PaperUnits', 'centimeters');  
      set(H,'Paperposition',newfpos);
      set(Ha,'Units','centimeters');
      set(Ha,'Position',newapos);

function wysiwyg
%       Dan(K) Braithwaite, Dept. of Hydrology U.of.A  11/93
unis = get(gcf,'units');
ppos = get(gcf,'paperposition');
set(gcf,'units',get(gcf,'paperunits'));
pos = get(gcf,'position');
pos(3:4) = ppos(3:4);
set(gcf,'position',pos);
set(gcf,'units',unis);

function MinorTick(N);
% N is the number of parts of division;
% 'box off' should be used before you use this function;
if nargin==0;
    N = 1;
end
N = N+1;
set(gca,'TickLength',[0.015 0.05]);
xx=max(xlim);xn=min(xlim);
yx=max(ylim);yn=min(ylim);
xt=get(gca,'xtick');
yt=get(gca,'ytick');
L=strcmp(get(gca,'box'),'on');
xs=xn:(xt(2)-xt(1))/N:xx;
ys=yn:(yt(2)-yt(1))/N:yx;
Lx=[yx-yn]*0.009;
Ly=[xx-xn]*0.008;
hold on;
for f=xs;
    plot([f,f],[yn,yn+Lx],'k');
    if L==1;plot([f,f],[yx,yx-Lx],'k');end;
end
for f=ys;
    plot([xn,xn+Ly],[f,f],'k');
    if L==1;plot([xx,xx-Ly],[f,f],'k');end; 
end
hold off;

function relegend(Hl,matx)
% change the form of the legend
ml = matx(1);
nl = matx(2);
[legh,objh,outh,outm] = legend(Hl);
suml = (length(objh)-1)/2;
Hline  = findall(Hl,'type','line');
xdat1 = get(Hline(2),'xdata');
ydat1 = get(Hline(1),'ydata');
ydat2 = get(Hline(3),'ydata');
deltm = ydat2-ydat1;
spacm = 1-deltm*length(Hline)/2;
newdeltm = deltm/(spacm+ml*deltm);
newspacm = spacm/(spacm+ml*deltm);
boxsz = get(objh(1),'extent');

oldlpos = get(Hl,'Position');
sH      = get(Hl,'Children');
% newlpos = [oldlpos(1:2),nl*oldlpos(3),(ml*deltm+spacm)*oldlpos(4)]; % position 
newlpos = [oldlpos(1:2),2*nl*oldlpos(3),1.5*(ml*deltm+spacm)*oldlpos(4)]; % 2* I am here
set(Hl,'Position',newlpos);
axh = axes('Position',newlpos);
hold on
axis([0,1,0,1]);
set(gca,'xtick',[],'ytick',[]);
box on;
for m=1:ml
    for n=1:nl
        k = (m-1)*nl+n;
        if k<=suml
         % row
%          tline(m,n) = plot(xdat1/nl+(m-1)/nl,[1,1]*(1-(newspacm+newdeltm)/2-(n-1)*newdeltm));
%          lline(m,n) = plot([mean(xdat1)/nl mean(xdat1)/nl]+(m-1)/nl,[1,1]*(1-(newspacm+newdeltm)/2-(n-1)*newdeltm));
%          ltext(m,n) = text(boxsz(1)/nl+(m-1)/nl,1-(newspacm+newdeltm)/2-(n-1)*newdeltm,char(outm(k)));
%          set(lline(m,n),'marker',get(sH(2*suml+1-2*k),'marker'),'markersize',get(sH(2*suml+1-2*k),'markersize')...
%                );
%             ,'MarkerFaceColor',[1 1 1]); % marker face invisible
%          set(lline(m,n),'color', get(sH(2*suml+2-2*k),'color'), 'LineStyle', 'none');
%          set(tline(m,n),'color', get(sH(2*suml+2-2*k),'color'), 'LineStyle', get(sH(2*suml+2-2*k),'Linestyle'));
%          set(ltext(m,n),'fontsize',get(sH(2*suml+1),'fontsize'));
         % column
         tline(m,n) = plot(xdat1/nl+(n-1)/nl,[1,1]*(1-(newspacm+newdeltm)/2-(m-1)*newdeltm));
         lline(m,n) = plot([mean(xdat1)/nl mean(xdat1)/nl]+(n-1)/nl,[1,1]*(1-(newspacm+newdeltm)/2-(m-1)*newdeltm),'color',[1,1,1]);
         ltext(m,n) = text(boxsz(1)/nl+(n-1)/nl,1-(newspacm+newdeltm)/2-(m-1)*newdeltm,char(outm(k)));
         set(lline(m,n),'marker',get(sH(2*suml+1-2*k),'marker'),'markersize',get(sH(2*suml+1-2*k),'markersize')...
         );
%          ,'MarkerFaceColor',[1 1 1]); % marker face invisible
         set(lline(m,n),'color', get(sH(2*suml+2-2*k),'color'), 'LineStyle', 'none');
         set(tline(m,n),'color', get(sH(2*suml+2-2*k),'color'), 'LineStyle', get(sH(2*suml+2-2*k),'Linestyle'));
         set(ltext(m,n),'fontsize',get(sH(2*suml+1),'fontsize'));
        end
    end
end
hold off
delete(Hl);
set(axh,'visible','off');