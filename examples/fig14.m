% fig14

clear all
% load syme_data
file = 'E:\experProc\pulsac\081103\'; 
syme_file = strcat(file,'intss15c5l');  
syme_temp = load(syme_file);                     
syme_cell = struct2cell(syme_temp);
[syme_data] = deal(syme_cell{:});
clear syme_file syme_temp syme_cell
syme_data = double(syme_data);
lm = size(syme_data,2);
for i = 1:lm
    dt = squeeze(syme_data(:,i,1:200));
    da{i} = flipud(dt)*255/65535; 
end
[mxi,mxj] = size(da{1});

% prepare 
ppm = 13.5; % pixel/mm
lns = 1;  % pixels for z spacing
pit = 2;  % points of interpolation
pzz = 3;  % times after zz interpolation

dy = 1/ppm/(pit+1);
dz = lns/ppm;
yy = 0:dy:(mxj-1)*dy;
zz = 0:dz:(mxi-1)*dz;
zi = 0:dz/pzz:(mxi-1)*dz;
[CY,CZ] = meshgrid(yy,zz);
[YI,ZI] = meshgrid(yy,zi);
   
mn = [2 3 0.012 0.01]; % m n lg cg   
% plot      
close all
hf = figure(1);
ha = gca;
set(gca,'XColor',[1 1 1],'XTick',[]) 
set(gca,'YColor',[1 1 1],'YTick',[]) 
% set(gca,'visible','off')
set(gcf,'Color',[1 1 1]);
% label, legend
xlabel 'Lateral position (mm)'
ylabel 'Distance from cathode (mm)'
set(gcf,'DefaultLineLineWidth',1.0) % 0.8
% colorbar;
pa = get(gca,'Position'); % [l b w h]
aps = zeros(1,4);
aps(3) = (pa(3)+mn(4))/mn(2)-mn(4); % aw;
aps(4) = (pa(4)+mn(3))/mn(1)-mn(3); % ah;
cp = [pa(1)+pa(3)+0.015 pa(2) 0.03 pa(4)]; % cl cb cw ch
colr = jet(64);
colormap(colr(15:64,:)); 
% colormap(gray(256))
hc = colorbar;
set(hc,'position',cp);
cy = get(hc,'ylim');
tvv = 0:50:255;
cmax = 250;
ytk = (cy(2)-cy(1))*(tvv+10)/(cmax+10);
set(hc,'YTick',ytk);
set(hc,'YTickLabel',tvv);
set(hf,'CurrentAxes',hc); % for ylabel
set(gca,'fontsize',7);    % tick labels 8 point
ylabel 'Intensity (a.u.)'
tx{1} = '(1)';
tx{2} = '(2)';
tx{3} = '(3)';
tx{4} = '(4)';
tx{5} = '(5)';
tx{6} = '(6)';
for i = 1:mn(1)
    for j = 1:mn(2)
        aps(1) = pa(1)+(j-1)*(aps(3)+mn(4));     % al
        aps(2) = pa(2)+(mn(1)-i)*(aps(4)+mn(3)); % ab
        k = (i-1)*mn(2)+j;
        axes('position',aps);
        vmax = max(max(da{k}));
        valu = 5:(vmax-5)/20:vmax; % 20
        valu = [-10 valu];
        tper = interp2(CY,CZ,da{k},YI,ZI,'cubic');
        [C,h,CF]= contourf(YI,ZI,tper,valu);
        set(h,'EdgeColor','none')
        caxis([-5 cmax])
        set(gca,'tickdir','out');
        set(gca,'ticklength',[0.015 0.025]);
    if i==mn(1) & j==1
        hxy = gca;
        xlblv = get(get(ha,'xlabel'),'string');
        ylblv = get(get(ha,'ylabel'),'string');
        delete(get(ha,'xlabel'));
        delete(get(ha,'ylabel'));
        pxl = get(get(gca,'xlabel'),'position'); % xlabel
        pxl(1) = pxl(1)+5.8;
        pxl(2) = pxl(2)+0.1;
        set(get(gca,'xlabel'),'position',pxl);
        pyl = get(get(gca,'ylabel'),'position'); % ylabel
        pyl(1) = pyl(1)+0.3;
        pyl(2) = pyl(2)+6.5;
        set(get(gca,'ylabel'),'position',pyl);
        XLabel(xlblv);
        YLabel(ylblv);
    end
    % set 
    set(gca,'fontsize',7);                 % tick labels 8 point
    set(get(gca,'XLabel'),'FontSize',8);   % xaxis label 9 point
    set(get(gca,'YLabel'),'FontSize',8);   % yaxis label 9 point 
   
    if i~=mn(1)
%         set(gca,'XColor',[1 1 1],'XTick',[]) 
        set(gca,'XTick',0:5);
        set(gca,'XTickLabel',[]);
    else
        set(gca,'XTick',0:5);
    end
    if j~=1
%         set(gca,'YColor',[1 1 1],'YTick',[]) 
        set(gca,'YTick',-0.5:2:15);
        set(gca,'YTickLabel',[]);
    else
        set(gca,'YTick',-0.55:2:15);
        set(gca,'YTickLabel',8:-2:-4);
    end
    text(4.3,10.3,tx{k},'FontSize',8)
    box off
end
end

% print figure
set(hf,'CurrentAxes',ha)
opts = struct('lrdu',[3 4 0 0],'figsize',[12.5 8.2]);
printfig(hf,'fig14',opts);