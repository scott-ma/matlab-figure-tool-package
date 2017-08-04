% example of one fig
% -------------------

clear all
tt   = 0:3:60;
val(1,:) = tt;
val(2,:) = 4*tt;
val(3,:) = 2*tt;
val(4,:) = 1*tt;

% plot   
close all
fgx(1);
plot(tt,val(1,:),'-k',tt,val(2,:),'--b',tt,val(3,:),':r',tt,val(4,:),'-.g')
% label, legend
fsz = [8 7]; 
xtk = [0 10 60];
ytk = [0 50 250];
axy = [0 62 0 260];
% set
XLabel('Time (s)');
YLabel('Arc current (A)');
hl = legend('No','Just go out','Yes','But this is OK!',0);
% print figure
opts = struct('lbrt',[0 0 3 2],'figsize',fsz,'ticksize',[0.02 0.025 0],...
              'xlbl',xtk,'ylbl',ytk,'axis',axy,...
              'legend',[hl 0.17 0.5 7],'relegend',[3 2 0.2]);
printfig('fig1',opts);
% 
% apos = fgx('apos');
% fgx(2);
% val = val/10;
% plot(tt,val(1,:),'-k',tt,val(2,:),'--b',tt,val(3,:),':r',tt,val(4,:),'-.g')
% % label, legend
% fsz = [8 7]; 
% xtk = [0 10 60 0 62];
% ytk = [0 5  25 0 26];
% % set
% XLabel('Time (s)');
% YLabel('Arc current (A)');
% hl = legend('No','Just go out','Yes','But this is OK!',0);
% % print figure
% opts = struct('lbrt',[0 0 3 2],'figsize',fsz,'ticksize',[0.02 0.025 1],...
%               'xlbl',xtk,'ylbl',ytk,'axis',[xtk(4:5) ytk(4:5)],'aps',apos,...
%               'legend',[hl 0.17 0.5 7],'relegend',[3 2 0.2]);
% printfig('fig2',opts);