% % example of two column figs
% % --------------------------
% clear all
% tt   = 0:60;
% val(1,:) = 4*tt;
% val(2,:) = tt/2;
% 
% % plot   
% close all
% figure(1);
% mn = [2 1 0 0]; % subfigs; m n lingap colgap 
% % label, legend
% fsz = [8 7]; 
% printfig('set','figsize',fsz);
% delete(gca);
% xlb =  'Time (s)';
% ylb = {'Arc current (A)','Arc voltage (V)'};
% xtk = [0 10 60];
% ytk = [0 50 250; 13 2 21];
% axy = [0 60 0 250; 0 60 13 22];
% % apos
% aps = zeros(1,4);
% aps(3) = (0.8+mn(4))/mn(2)-mn(4); % aw;
% aps(4) = (0.8+mn(3))/mn(1)-mn(3); % ah;
% for i = 1:mn(1)
%     for j = 1:mn(2)
%         aps(1) = 0.1+(j-1)*(aps(3)+mn(4));     % al
%         aps(2) = 0.1+(mn(1)-i)*(aps(4)+mn(3)); % ab
%         axes('position',aps); % pos of subfig
%         k = (i-1)*mn(2)+j;
%         plot(tt,val(k,:))
%         % set xylbl
%         if i==mn(1) % set xlbl
%             XLabel(xlb);
%         end
%         YLabel(ylb(k));
%         % set other parameters
%         opts = struct('ticksize',[0.02 0.025 1],...
%                       'xlbl',xtk,'ylbl',ytk(k,:),'axis',axy(k,:));
%         printfig('set',opts);
%         if i~=mn(1) % delete xticklbl
%             set(gca,'XTickLabel',[]);
%         end
%     end
% end
% % print figure
% printfig('fig2','lbrt',[0 0 0 0],'figsize',fsz);

% fgx1: one x-axis for two axes
% clear all
% xx   = 0:60;
% val(1,:) = xx;
% val(2,:) = 3*xx;
% val(3,:) = xx/3;
% 
% % plot 
% close all
% fgx(2);
% fsz = [8 6]; 
% xyl = {'Time (s)','Arc current (A)','Arc voltage (V)'};
% xyt = [0 10 60 0 60; 0 50 250 0 250; 13 2 21 13 22];
% tks = [0.02 0.025 1];
% fgx('fgx1',val,fsz,xyl,xyt,tks);
% % print fig
% printfig('figx1','lbrt',[0 0 0 0],'figsize',fsz);

% -----------------------------
% fgy2: one x-axis & two y-axes
clear all
tt   = 0:3:60;
val(1,:) = tt;
val(2,:) = 4*tt;
val(3,:) = tt.^0.5;
val(4,:) = tt.^0.8;
% plot 
close all
fgx(1);
fsz = [8 6]; 
xyl = {'Time (s)','Arc current (A)','Arc voltage (V)'};
xyt = [0 10 60 0 60; 0 50 250 0 250; 0 5 30 0 30];
tks = [0.02 0.025 1];
fgx('fgy2','a',fsz,xyl,xyt,tks);
plot(tt,val(2,:))
fgx('fgy2','ab',fsz,xyl,xyt,tks);
plot(tt,val(3:4,:))
fgx('fgy2','b',fsz,xyl,xyt,tks);
% print fig
printfig('figy2','lbrt',[0 0 2 1],'figsize',fsz);

% % -----------------------------
% % fgbx: break x-axis
% clear all
% xx1 = 0:10;
% yy1 = xx1;
% xx2 = 10:0.01:11;
% yy2 = 10+2*sin(50*xx2);
% xx3 = 11:20;
% yy3 = 21-xx3;
% % plot 
% close all
% fgx(1);
% fsz = [9 8]; 
% rxy = [0.2 0.4 0.2];
% xyl = {'Time (s)','Arc current (A)'};
% xyt = [0 2 9 0 10; 10 0.2 10.9 10 11; 11 3 20 11 20; 0 2 12 0 12];
% tks = [0.02 0.025 1];
% fgx('fgbx','a',fsz,rxy);
% plot(xx1,yy1)
% fgx('fgbx','ab',xyt,tks,rxy,1);
% plot(xx2,yy2)
% fgx('fgbx','ab',xyt,tks,rxy,2);
% plot(xx3,yy3)
% fgx('fgbx','b' ,xyt,tks,xyl);
% % print fig
% printfig('fgbx','lbrt',[0 0 0 2],'figsize',fsz,'bkxy','x');

% % -----------------------------
% % fgby: break y-axis
% clear all
% xx1 = 0:0.1:10;
% yy1 = xx1;
% yy2 = 10+2*sin(50*[10:0.01:11]);
% yy3 = 12:0.1:22;
% % plot 
% close all
% fgx(1);
% fsz = [8 10]; 
% rxy = [0.3 0.2 0.3];
% xyl = {'Time (s)','Arc current (A)'};
% xyt = [0 2 10 0 10; 0 2 10 0 10; 10 1 12 10 12; 12 2 22 12 22];
% tks = [0.02 0.025 1];
% fgx('fgby','a',fsz,rxy);
% plot(xx1,yy1)
% fgx('fgby','ab',xyt,tks,rxy,1);
% plot(xx1,yy2)
% fgx('fgby','ab',xyt,tks,rxy,2);
% plot(xx1,yy3)
% fgx('fgby','b' ,xyt,tks,xyl);
% % print fig
% printfig('fgby','lbrt',[0 4 0 0],'figsize',fsz,'bkxy','y');