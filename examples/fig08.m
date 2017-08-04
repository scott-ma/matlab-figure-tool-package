% fig08

file = 'G:\eProc\timevary\';
ftpr = strcat(file,'TMI200upAr794');
fps  = 955;
sers = round([1 2 4.4]*fps+196);
li1 = 2;
li2 = 6;
lll = [2,4,6];

% load tpr
tema_temp = load(ftpr);                     
tema_cell = struct2cell(tema_temp);
[tema_data] = deal(tema_cell{:});  
tema_data = double(tema_data);
da1s = squeeze(tema_data(lll(1),:,1));
da2s = squeeze(tema_data(lll(2),:,1));
da3s = squeeze(tema_data(lll(3),:,1));
ttt  = -196/fps:1/fps:(length(da1s)-1-196)/fps;
for i = 1:5
    if i==1
        tpr = squeeze(tema_data(li1,sers+i,:));
        tpr2= squeeze(tema_data(li2,sers+i,:));
        tpz = squeeze(tema_data(:,sers+i,1));
    else
        tpr = squeeze(tema_data(li1,sers+i,:))+tpr;
        tpr2= squeeze(tema_data(li2,sers+i,:))+tpr2;
        tpz = squeeze(tema_data(:,sers+i,1))+tpz;
    end
end
tpr = tpr/5;
tpr2= tpr2/5;
tpz = tpz'/5;
[mi,mj] = size(tpr);
rr = (0:mj-1)/33;
[mi,mj] = size(tpz);
zz = (0:mj-1)/2;
zz(1) = 0.1;
clear ftpr tema_temp tema_cell tema_data

% plot
close all
figure(1) 
plot(ttt,da1s/1000,'-k',ttt,da2s/1000,'-b',ttt,da3s/1000,'-r')
% label, legend
xlabel '{\itt} (s)'
ylabel 'Temperature (1000 K)'
text(1.9,24.45,'{\itz} = 0.5 mm','fontsize',7);
text(1.6,20.25,'{\itz} = 1.5 mm','fontsize',7);
text(1.5,17.80,'{\itz} = 2.5 mm','fontsize',7);
text(0.8,0.95,'(a)','parent',gca,'units','normalized','fontsize',9);
% print figure
opts = struct('lbrt',[0 0 -148.15 0],'figsize',[8.5 5.8],'ticksize',[0.02 0 1],'tickxylblfs',[7 8],...
              'xlbl',[0 2 4],'ylbl',[16 2 26],'axis',[0 4.5 15 26]);
printfig('fig08',opts);

% plot fig [1,2]
aps = fgx('apos');
aps1 = aps(1,:);
aps1(1) = aps(1,1)+aps(1,3)+0.01;
axes('position',aps1);
tpz = tpz/1000;
plot(zz(2:end),tpz(1,2:end),'-ok',zz(2:end),tpz(2,2:end),'--^b',zz(2:end),tpz(3,2:end),':dr',...
    'MarkerFaceColor',[1 1 1],'MarkerSize',3,'linewidth',0.7);
% label, legend
xlabel '{\itz} (mm)'
text(0.63,15.2,'Cathode','Rotation',90,'fontsize',7);
text(0.8,0.95,'(b)','parent',gca,'units','normalized','fontsize',9);
hl = legend('{\itt} = 1.0 s','{\itt} = 2.0 s','{\itt} = 4.4 s',0); 
opts = struct('ticksize',[0.02 0 1],'tickxylblfs',[7 8],...
              'xlbl',[1 1 3],'ylbl',[16 2 26],'axis',[0.4 3.1 15 26],...
              'legend',[hl 0.44 0.65 7],'relegend',[3 1]);
printfig('set',opts); 
set(gca,'YTickLabel',[]);

% plot fig [1,3]
aps2 = aps1;
aps2(1) = aps1(1)+aps1(3)+0.05;
axes('position',aps2);
tpr = tpr/1000;
tpr2=tpr2/1000;
plot(rr,tpr(1,:),'-k',rr,tpr(2,:),'--b',rr,tpr(3,:),':r','linewidth',0.7);
hold on
plot(rr,tpr2(1,:),'-k',rr,tpr2(2,:),'--b',rr,tpr2(3,:),':r','linewidth',0.7);
hold off
% label, legend
xlabel '{\itr} (mm) '
text(0.8,17.9,'{\itz} = 0.5 mm','fontsize',7);
text(1.8,14.3,'2.5 mm','fontsize',7);
text(0.8,0.95,'(c)','parent',gca,'units','normalized','fontsize',9); 
hl = legend('{\itt} = 1.0 s','{\itt} = 2.0 s','{\itt} = 4.4 s',0); 
% print figure
opts = struct('ticksize',[0.02 0 1],'tickxylblfs',[7 8],...
              'xlbl',[0 1 3],'ylbl',[10 4 26],'axis',[0 3 10 26],...
              'legend',[hl 0.76 0.65 7],'relegend',[3 1]);
printfig('set',opts);
fgs('pfg','fig08'); % print fig
