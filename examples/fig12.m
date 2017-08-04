% fig12

ppm = 12; % pixel/mm
lns = 6;  % pixels for z spacing
pit = 2;  % points of interpolation
pzz = 3;  % times after zz interpolation

file = 'E:\experProc\comp\imlc\';
% load tpr_l
name = 'TMI200Ar794';
tper_file = strcat(file,name);  
tper_temp = load(tper_file);                     
tper_cell = struct2cell(tper_temp);
[tper_data] = deal(tper_cell{:});   
clear tper_file tper_temp tper_cell
tper_data = double(tper_data)/1000;
tper_data = flipud(tper_data);
[mxi,mxj] = size(tper_data);
% load tpr_c
name = 'TAI200Ar794';
tper_file = strcat(file,name);  
tper_temp = load(tper_file);                     
tper_cell = struct2cell(tper_temp);
[tper_datc] = deal(tper_cell{:});   
clear tper_file tper_temp tper_cell
tper_datc = double(tper_datc)/1000;
tper_datc = flipud(tper_datc);

% prepare 
dy = 1/ppm/(pit+1);
dz = lns/ppm;
yy = 0:dy:(mxj-1)*dy;
zz = 0:dz:(mxi-1)*dz;
zi = 0:dz/pzz:(mxi-1)*dz;
[CY,CZ] = meshgrid(yy,zz);
[YI,ZI] = meshgrid(yy,zi);
tpel = interp2(CY,CZ,tper_data,YI,ZI,'cubic');
tpec = interp2(CY,CZ,tper_datc,YI,ZI,'cubic');

txta = '{\itI} = 200 A';
txtb = '\theta = 60 deg';
valu = [10:18 20];
% plot
close all
fgn(1)
hold on
plot([-1 -1], [-1 -1],'-b')
plot([-1 -1], [-1 -1],'--r')
[C,h]= contour(YI,ZI,tpel,valu,'-b');
[C,h]= contour(YI,ZI,tpec,valu,'--r');
hold off
box on
set(h,'linewidth',0.8)
text(5.6,1.4,'10','fontsize',7,'HorizontalAlignment','center',...
        'BackgroundColor',[1 1 1],'Margin',0.01);
text(3.1,1.7,'12','fontsize',7,'HorizontalAlignment','center',...
        'BackgroundColor',[1 1 1],'Margin',0.01);
text(1.4,2.0,'14','fontsize',7,'HorizontalAlignment','center',...
        'BackgroundColor',[1 1 1],'Margin',0.01);
text(0.58,2.8,'16','fontsize',7,'HorizontalAlignment','center',...
        'BackgroundColor',[1 1 1],'Margin',0.01);
text(0.3,3.7,'18','fontsize',7,'HorizontalAlignment','center',...
        'BackgroundColor',[1 1 1],'Margin',0.01);
text(0.33,4.4,'20','fontsize',7,'HorizontalAlignment','center',...
        'BackgroundColor',[1 1 1],'Margin',0.01);   
% label, legend
xlabel 'Radius (mm)'
ylabel 'Distance from anode (mm)'
text(4.46,3.8,txta,'fontsize',7);
text(4.41,3.5,txtb,'fontsize',7);
hl = legend('Line','Line+Continuum',0); 
% print figure
opts = struct('lrdu',[3 0 12 0],'figsize',[6.5 6.5],'ticksize',[0.02 0.025 1], ...
              'xlbl',[0 1 6],'ylbl',[0 1 5],'axis',[0 6 0 5],...
              'legend',[hl,0.62,0.9,8],'relegend',[2 1]);
printfig(gcf,'fig12',opts);