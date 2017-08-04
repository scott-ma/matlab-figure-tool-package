% fig06

file = 'E:\experProc\pulsdc\081216\';
af = strcat(file,'As3f2Tg794');
tf = strcat(file,'TMs3f2Ar794');
ll = [1:7];

af = load(af);                     
af_cell = struct2cell(af);
[af] = deal(af_cell{:});   
af = double(af);
int = af(9,:)/65535*2550;
tf = load(tf);                     
tf_cell = struct2cell(tf);
[tf] = deal(tf_cell{:});   
tf = double(tf);
tpr = tf(3,:,1);
t = (0:(size(int,2)-1))/955;

% plot
close all
figure(1)
[ax,h1,h2]=plotyy(t,int,t,tpr);
grid on
set(h1,'lineStyle','-');
set(h2,'lineStyle','-');
set(h1,'color','b');
set(h2,'color','r');
axes(ax(1));
hold on
plot(t,int,'-b',t,tpr,'--w');
hold off
set(ax(1),'ycolor','k');
set(ax(2),'ycolor','k');
% label, legend
xlabel('Time (s)');
set(get(ax(1),'Ylabel'),'String','Intensity (a.u.)');
set(get(ax(2),'Ylabel'),'String','Temperature (K)');
text(0.9,230,'{\itI}_p = 150 A, {\itI}_b = 120 A; {\itf} = 2 Hz','fontsize',7);
text(2.55,180,'\rightarrow','fontsize',15);
text(2.45,200,'{\itz} = 1.0 mm','fontsize',7);
text(0.3,60 ,'\leftarrow','fontsize',15);
text(0.23,75,'{\itz} = -2.0 mm','fontsize',7);
% print figure
opts = struct('lrdu',[4 1 -4 -1],'figsize',[10.3 6],'ticksize',[0.02 0.025 1], ...
              'xlbl',[0 0.5 3],'ylbl',[0 50 300],'ylb2',[15000 1000 20000],'axis',[0 3 0 250 15000 20000]);
msizefig(gcf,'fig06',opts);