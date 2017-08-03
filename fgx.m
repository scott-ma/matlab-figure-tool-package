function r_out = fgx(varargin)
%mainly for plotting different types of figures 
% Copyright 2009 lshuily@gmail.com
%  fgx(n)
%  fgx('apos')
%  fgx('fgx1',val,xyl,xyt,tks,clr)
%  fgx('fgy2',opt)
%  fgx('fgbx',opt)
%  fgx('fgby',opt)
% =============================================
% r_out = [];
if (nargin < 1)
  error('Too few input arguments');
end
swh = varargin{1};    % case
if isnumeric(swh)
    fgn = swh;
    swh = 'fgn';
end

switch swh
% ---------------------------------------------
case 'fgn'    % fgx(1)
    figure(fgn);    % create a figure with default linewidth
    set(gcf,'DefaultLineLineWidth',0.8);   % 0.8
case 'apos'   % fgx('apos')
    if nargout > 0
        r_out = apos;          % get positions of all axeses
    end
case 'fgx1'   % fgx('fgx1',val,xyl,xyt,tks,clr)
    if (nargin < 2)
        error('Too few input arguments');
    else
        opt = varargin(2:end);
        if length(opt)~=7
             error('Too few input arguments');
        else
             fgx1(opt);
        end
    end
case 'fgy2'    % fgx('fgy2','seq',val,xyl,xyt,tks)
    if (nargin < 6)
        error('Too few input arguments');
    else
        fgy2(varargin(2:end));
    end
case 'fgbx'    % fgx('fgbx','seq',val,xyl,xyt,tks)
    if (nargin < 2)
        error('Too few input arguments');
    else
        fgbx(varargin(2:end));
    end    
case 'fgby'    % fgx('fgby','seq',val,xyl,xyt,tks)
    if (nargin < 2)
        error('Too few input arguments');
    else
        fgby(varargin(2:end));
    end     
end % switch


%
%  Local Functions
% ----------------
function r_t = apos()
r_t = [];
lps = [];
H = get(0,'CurrentFigure');
if ~isempty(H)
allAxes = findall(H,'type','axes');
if ~isempty(allAxes)
    lgthAxe = length(allAxes);
    aps = zeros(lgthAxe,4);
    set(allAxes,'units','normalized');
    i = 1;
    for k = 1:lgthAxe
        ax = allAxes(k);
        if ~strcmpi(get(ax,'Tag'),'legend')
            aps(i,:) = get(ax,'position');
            i = i+1;
        else
            lps = get(ax,'position');
        end
    end
    if i>1
        aps = aps(1:i-1,:);
    else
        aps = [];
    end
    if isempty(lps)
        r_t = aps;
    else
        r_t = [aps; lps];
    end
end
end

% ----------------
function fgx1(opt)
seq = opt{1}; % select
fsz = opt{2}; % fig size
xyl = opt{3}; % xy ticks
xyt = opt{4}; % xy label
tks = opt{5}; % ticksize
dis = opt{6}; % distance
fn  = opt{7};

% axes position
mn  = [2 1 0 0];                  % subfigs; m n lingap colgap 
aps = zeros(1,4);
aps(3) = (0.8+mn(4))/mn(2)-mn(4); % aw;
aps(4) = (0.8+mn(3))/mn(1)-mn(3); % ah;
aps(1) = 0.1;                     % al
aps(2) = 0.1+(fn-1)*(aps(4)+dis); % ab

%% different operations
switch seq
case 'a'     % for plot a
    printfig('set','figsize',fsz);
    axes('position',aps); % pos of subfig
case 'ab'    % set a & for b
    fn = fn-1;
    % set xylbl
    if fn==1 % set xlbl
        XLabel(xyl{1});
    end
    YLabel(xyl(fn+1));
    % set other parameters
    opts = struct('ticksize',tks,...
        'xlbl',xyt(1,1:3),'ylbl',xyt(fn+1,1:3),...
        'axis',[xyt(1,4:5) xyt(fn+1,4:5)]);
    printfig('set',opts);
    if fn~=1 % delete xticklbl
        set(gca,'XTickLabel',[]);
    end
    axes('position',aps); % pos of subfig
case 'b'     % set b
    % set xylbl
    if fn==1 % set xlbl
        XLabel(xyl{1});
    end
    YLabel(xyl(fn+1));
    % set other parameters
    opts = struct('ticksize',tks,...
        'xlbl',xyt(1,1:3),'ylbl',xyt(fn+1,1:3),...
        'axis',[xyt(1,4:5) xyt(fn+1,4:5)]);
    printfig('set',opts);
    if fn~=1 % delete xticklbl
        set(gca,'XTickLabel',[]);
    end
end  % switch      

% ----------------
function fgy2(opt)
    seq = opt{1};
    fsz = opt{2}; % fig size
    xyl = opt{3}; % xy ticks
    xyt = opt{4}; % xy label
    tks = opt{5}; % ticksize
    if length(opt) > 5
        sgnbox = 0;
    else
        sgnbox = 1;
    end
%% different operations
  switch seq
case 'a'    % for plot a
    printfig('set','figsize',fsz);
    newplot;
%     axes('position',[0.13 0.12 0.9 0.81]);
    set(gcf,'nextplot','add')
case 'ab'    % set a & for b
    set(gca,'box','off','Tag','y')
    % set xylbl
    XLabel(xyl{1});
    YLabel(xyl{2});
    % set other parameters
    opts = struct('ticksize',tks,...
        'xlbl',xyt(1,1:3),'ylbl',xyt(2,1:3),...
        'axis',[xyt(1,4:5) xyt(2,4:5)]);
    printfig('set',opts);
    % axes 2
    axes('position',get(gca,'position'));
case 'b'    % set b
    set(gca,'YAxisLocation','right','color','none', ...
          'xgrid','off','ygrid','off','box','off','Tag','y');  
    % set xylbl
    YLabel(xyl{3});
    % set other parameters
    opts = struct('ticksize',tks,...
        'xlbl',xyt(1,1:3),'ylbl',xyt(3,1:3),...
        'axis',[xyt(1,4:5) xyt(3,4:5)]);
    printfig('set',opts);
    % temp axes for box on
    if sgnbox==1
        axes('position',get(gca,'position'));
        set(gca,'YAxisLocation','right','color','none', ...
            'xgrid','off','ygrid','off','box','on','Tag','x');   
        opts = struct('ticksize',tks,...
            'xlbl',xyt(1,1:3),'ylbl',xyt(3,1:3),...
            'axis',[xyt(1,4:5) xyt(3,4:5)]);
        printfig('set',opts);
        set(gca,'ytick',[]);
    end
end  % switch

% ----------------
function fxy = tlbx(lst) % get ratio for tick length set: fgbx
fxy = 1;
allAxes = findall(gcf,'type','axes');
for k = 1:length(allAxes)
    if strcmpi(get(allAxes(k),'Tag'),'t')
        tpos = get(allAxes(k),'position');
        apos = get(gca,'position');
        fpos = get(gcf,'position');
        tw = fpos(3)*tpos(3);
        th = fpos(4)*tpos(4);
        aw = fpos(3)*apos(3);
        ah = fpos(4)*apos(4);
        fxy = max(tw,th)/max(aw,ah);
        if nargin > 0
            delete(allAxes(k));
        end
        break
    end
end

% ----------------
function fgbx(opt)
seq = opt{1};
%% different operations
  switch seq
case 'a'    % for plot a
    if (length(opt) < 3)
        error('Too few input arguments');
    else
        fsz = opt{2}; % fig size
        rxy = opt{3}; % x-axis ratio
        printfig('set','figsize',fsz);
        newplot;
        set(gcf,'nextplot','add')
        apos = get(gca,'position');
        set(gca,'Tag','t');
        apos(3) = rxy(1);
        axes('position',apos);
    end
case 'ab'    % set a & for b
    if (length(opt) < 5)
        error('Too few input arguments');
    else
        xyt = opt{2}; % xy ticks
        tks = opt{3}; % ticksize
        rxy = opt{4}; % x-axis ratio
        axn = opt{5}; % ratio number
        if length(opt) > 5
            sgnbox = 0;
        else
            sgnbox = 1;
        end
        % set y-axis color & tag
        set(gca,'YColor','w','Tag','x');
        % tick length
        fxy = tlbx;
        % set other parameters
        opts = struct('ticksize',tks*fxy,...
            'xlbl',xyt(axn,1:3),'ylbl',xyt(end,1:3),...
            'axis',[xyt(axn,4:5) xyt(end,4:5)]);
        printfig('set',opts);
        set(gca,'ytick',[]);
        % axes b
        apos = get(gca,'position');
        apos(1) = apos(1)+apos(3);
        apos(3) = rxy(axn+1);
        axes('position',apos);
    end
case 'b'    % set b
    if (length(opt) < 4)
        error('Too few input arguments');
    else
        xyt = opt{2}; % xy ticks
        tks = opt{3}; % ticksize
        xyl = opt{4}; % xy label
        if length(opt) > 4
            sgnbox = 0;
        else
            sgnbox = 1;
        end
        if sgnbox==1
            set(gca,'box','on')
        else
            set(gca,'box','off')
        end
        % set y-axis color & tag
        set(gca,'YColor','w','Tag','x');
        % tick length
        fxy = tlbx(1);
        % set other parameters
        opts = struct('ticksize',tks*fxy,...
            'xlbl',xyt(end-1,1:3),'ylbl',xyt(end,1:3),...
            'axis',[xyt(end-1,4:5) xyt(end,4:5)]);
        printfig('set',opts);
        set(gca,'ytick',[]);
        % temp axes for box on
        ax  = findall(gcf,'type','axes');
        apos = get(ax(1),'position');
        apos(3:4) = apos(3:4)+apos(1:2);
        for k=2:length(ax)
            pos = get(ax(k),'position');
            pos(3:4) = pos(1:2)+pos(3:4);
            apos(1) = min(apos(1),pos(1));
            apos(2) = min(apos(2),pos(2));
            apos(3) = max(apos(3),pos(3));
            apos(4) = max(apos(4),pos(4));
        end
        apos(3:4) = apos(3:4)-apos(1:2);
        axes('position',apos);
        set(gca,'color','none','xgrid','off','ygrid','off','Tag','y');  
        xlabel(xyl{1});
        ylabel(xyl{2});
        if sgnbox==1
            set(gca,'box','on')
        else
            set(gca,'box','off')
        end
        % set
        opts = struct('ticksize',tks,...
            'xlbl',[0 0.2 1],'ylbl',xyt(end,1:3),...
            'axis',[0 1 xyt(end,4:5)]);
        printfig('set',opts);
    end
end  % switch

% ----------------
function fgby(opt)
seq = opt{1};
%% different operations
  switch seq
case 'a'    % for plot a
    if (length(opt) < 3)
        error('Too few input arguments');
    else
        fsz = opt{2}; % fig size
        rxy = opt{3}; % y-axis ratio
        printfig('set','figsize',fsz);
        newplot;
        set(gcf,'nextplot','add')
%         ax = gca;
        apos = get(gca,'position');
        set(gca,'Tag','t');
%         delete(ax);
        apos(4) = rxy(1);
        axes('position',apos);
    end
case 'ab'    % set a & for b
    if (length(opt) < 5)
        error('Too few input arguments');
    else
        xyt = opt{2}; % xy ticks
        tks = opt{3}; % ticksize
        rxy = opt{4}; % x-axis ratio
        axn = opt{5}; % ratio number
        if length(opt) > 5
            sgnbox = 0;
        else
            sgnbox = 1;
        end
        % set x-axis color & tag
        set(gca,'XColor','w','Tag','y');
        % tick length
        fxy = tlbx;
        % set other parameters
        opts = struct('ticksize',tks*fxy,...
            'xlbl',xyt(1,1:3),'ylbl',xyt(axn+1,1:3),...
            'axis',[xyt(1,4:5) xyt(axn+1,4:5)]);
        printfig('set',opts);
        set(gca,'xtick',[]);
        % axes b
        apos = get(gca,'position');
        apos(2) = apos(2)+apos(4);
        apos(4) = rxy(axn+1);
        axes('position',apos);
    end
case 'b'    % set b
    if (length(opt) < 4)
        error('Too few input arguments');
    else
        xyt = opt{2}; % xy ticks
        tks = opt{3}; % ticksize
        xyl = opt{4}; % xy label
        if length(opt) > 4
            sgnbox = 0;
        else
            sgnbox = 1;
        end
        if sgnbox==1
            set(gca,'box','on')
        else
            set(gca,'box','off')
        end
        % set x-axis color & tag
        set(gca,'XColor','w','Tag','y');
        % tick length
        fxy = tlbx(1);
        % set other parameters
        opts = struct('ticksize',tks*fxy,...
            'xlbl',xyt(1,1:3),'ylbl',xyt(end,1:3),...
            'axis',[xyt(1,4:5) xyt(end,4:5)]);
        printfig('set',opts);
        set(gca,'xtick',[]);
        % temp axes for box on
        ax  = findall(gcf,'type','axes');
        apos = get(ax(1),'position');
        apos(3:4) = apos(3:4)+apos(1:2);
        for k=2:length(ax)
            pos = get(ax(k),'position');
            pos(3:4) = pos(1:2)+pos(3:4);
            apos(1) = min(apos(1),pos(1));
            apos(2) = min(apos(2),pos(2));
            apos(3) = max(apos(3),pos(3));
            apos(4) = max(apos(4),pos(4));
        end
        apos(3:4) = apos(3:4)-apos(1:2);
        handt = axes('position',apos);
        set(gca,'color','none','xgrid','off','ygrid','off','Tag','x');  
        xlabel(xyl{1});
        ylabel(xyl{2});
        if sgnbox==1
            set(gca,'box','on')
        else
            set(gca,'box','off')
        end
        % get axes & x ticks
        allAxes = findall(gcf,'type','axes');
        ax = allAxes;
        j  = 0;
        for k = 1:length(allAxes)
            if strcmpi(get(allAxes(k),'Tag'),'y')
                j = j+1;
                ax(j) = allAxes(k);
            end
        end
        ax = ax(1:j);
        if j<2 | j>4
            error('Axes number error');
        else
            posmin = 10000;
            k = 0;
            for i = 1:j
                set(gcf,'CurrentAxes',ax(i));
                hylbl = ylabel(xyl{2});
                posyl = get(hylbl,'position');
                delete(hylbl);
                if posmin > posyl(1)
                    posmin = posyl(1);
                    k = i;
                end
            end
            allytick = get(ax(k),'ytick');
            k = 0;
            for i = 1:j
                if allytick(1) < xyt(i+1,1) | allytick(1) == xyt(i+1,1)
                    k = i;
                    break;
                end
            end
        end
        % set
        set(gcf,'CurrentAxes',handt);
        opts = struct('ticksize',tks,...
            'xlbl',xyt(1,1:3),'ylbl',xyt(k+1,1:3),...
            'axis',[xyt(1,4:5) xyt(k+1,4:5)]);
        printfig('set',opts);
    end
end  % switch

