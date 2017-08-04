% ???
file = 'H:\DESKTOP\w2\z\';

fsg = strcat(file,'s');

fa1 = strcat(file,'a1');
fa2 = strcat(file,'a2');
fa3 = strcat(file,'a3');

fd1 = strcat(file,'d1');
fd2 = strcat(file,'d2');
fd3 = strcat(file,'d3');

% load tpr
sg = load(fsg); 

a1 = load(fa1); 
a2 = load(fa2); 
a3 = load(fa3); 

d1 = load(fd1); 
d2 = load(fd2); 
d3 = load(fd3); 

% plot
close all
% plot fig s [1,1]
% ---------------------------------------------------------
figure(1) 
plot(sg,'-r','linewidth',0.5)
% label, legend
%ylabel 's'
text(0.45,0.85,'s','units','normalized','fontsize',8);
text(300,0.21,'Signal and Approximation (s)','fontsize',6);;
% print figure
opts = struct('lbrt',[0 -158 -125 15],'figsize',[8.8 8.0],'ticksize',[0.03 0 1],'tickxylblfs',[6 7],...
              'xlbl',[0 500 1750],'ylbl',[0 0.05 0.2],'axis',[0 1700 0 0.19]);
printfig('fig15',opts);
set(gca,'XTickLabel',[]);

% plot fig a3 [1,2]
% ---------------------------------------------------------
aps = fgx('apos');
aps1 = aps(1,:);
aps1(4) = 1.0*aps1(4);
aps1(2) = aps(1,2)-aps1(4)-0.013;
axes('position',aps1);
plot(a3,'-b','linewidth',0.5);
% label, legend
%ylabel 'a_3'
text(0.45,0.85,'a_3','units','normalized','fontsize',8);
opts = struct('ticksize',[0.02 0 1],'tickxylblfs',[6 7],...
              'xlbl',[0 500 1750],'ylbl',[0 0.05 0.14],'axis',[0 1700 0 0.14]);          
printfig('set',opts); 
set(gca,'XTickLabel',[]);

%% save pos for right subfig
aps0 = fgx('apos');

% plot fig a2 [1,3]
% ---------------------------------------------------------
aps = fgx('apos');
aps1 = aps(1,:);
aps1(2) = aps(1,2)-aps1(4)-0.013;
axes('position',aps1);
plot(a2,'-b','linewidth',0.5);
% label, legend
%ylabel 'a_2'
text(0.45,0.85,'a_2','units','normalized','fontsize',8);
opts = struct('ticksize',[0.02 0 1],'tickxylblfs',[6 7],...
              'xlbl',[0 500 1750],'ylbl',[0 0.05 0.15],'axis',[0 1700 0 0.16]);          
printfig('set',opts); 
set(gca,'XTickLabel',[]);

% plot fig a1 [1,4]
% ---------------------------------------------------------
aps = fgx('apos');
aps1 = aps(1,:);
aps1(2) = aps(1,2)-aps1(4)-0.013;
axes('position',aps1);
plot(a1,'-b','linewidth',0.5);
% label, legend
%ylabel 'a_1'
text(0.45,0.85,'a_1','units','normalized','fontsize',8);
opts = struct('ticksize',[0.02 0 1],'tickxylblfs',[6 7],...
              'xlbl',[0 500 1750],'ylbl',[0 0.06 0.2],'axis',[0 1700 0 0.22]);          
printfig('set',opts); 

%% right figs
% ==========================================================
% ----------------------------------------------------------
% plot fig d3 [2,2]
aps1 = aps0(1,:);
aps1(1) = aps1(1)+aps1(3)+0.09;
axes('position',aps1);
plot(d3,'-g','linewidth',0.5);
% label, legend
%ylabel 'd_3'
text(0.45,0.85,'d_3','units','normalized','fontsize',8);
opts = struct('ticksize',[0.02 0 1],'tickxylblfs',[6 7],...
              'xlbl',[0 500 1750],'ylbl',[-0.03 0.03 0.04],'axis',[0 1700 -0.05 0.04]);          
printfig('set',opts); 
set(gca,'XTickLabel',[]);

% plot fig d2 [2,3]
% ---------------------------------------------------------
aps = fgx('apos');
aps1 = aps(1,:);
aps1(2) = aps(1,2)-aps1(4)-0.013;
axes('position',aps1);
plot(d2,'-b','linewidth',0.5);
% label, legend
%ylabel 'd_2'
text(0.45,0.85,'d_2','units','normalized','fontsize',8);
opts = struct('ticksize',[0.02 0 1],'tickxylblfs',[6 7],...
              'xlbl',[0 500 1750],'ylbl',[-0.05 0.05 0.06],'axis',[0 1700 -0.08 0.08]);          
printfig('set',opts); 
set(gca,'XTickLabel',[]);

% plot fig d1 [2,4]
% ---------------------------------------------------------
aps = fgx('apos');
aps1 = aps(1,:);
aps1(2) = aps(1,2)-aps1(4)-0.013;
axes('position',aps1);
plot(d1,'-g','linewidth',0.5);
% label, legend
%ylabel 'd_1'
text(0.45,0.85,'d_1','units','normalized','fontsize',8);
opts = struct('ticksize',[0.02 0 1],'tickxylblfs',[6 7],...
              'xlbl',[0 500 1750],'ylbl',[-0.1 0.05 0.99],'axis',[0 1700 -0.12 0.09]);          
printfig('set',opts); 

fgs('pfg','fig15'); % print fig