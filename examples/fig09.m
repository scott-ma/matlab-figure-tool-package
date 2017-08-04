%fig09

file = 'G:\eProc\pulsdc\081216\';
tf = strcat(file,'TMs3f2Ar794');
af = 'intp_dat';
ll = [1:7];

af = proc('load',af);                     
int= double(af)*255;
tf = proc('load',tf); 
tpr= double(tf([2,4],:,1));
t  = (0:(size(int,2)-1))/955+0.5;

% move aver
s1 = [1   117 478 594  957 1072 1434 1550 1911 2027 2388 2504 2867 2982 3343 3457 3822 3936 4299 4414];
s2 = [118 477 596 953 1076 1431 1553 1908 2031 2386 2506 2863 2984 3341 3465 3819 3941 4297 4416 4636];
stpr = tpr;
int(4,:) = int(2,:)-int(3,:)/2;
sint = int;
span = 3;
for i = 1:length(s1)/2
    for j = 1:size(int,1)        
%         sint(j,s1(2*i-1):s1(2*i)) = smooth(int(j,s1(2*i-1):s1(2*i)),span,'moving')';
%         sint(j,s2(2*i-1):s2(2*i)) = smooth(int(j,s2(2*i-1):s2(2*i)),span,'moving')';
        sint(j,s1(2*i-1):s1(2*i)) = medfilt1(int(j,s1(2*i-1):s1(2*i)),span);
        sint(j,s2(2*i-1):s2(2*i)) = medfilt1(int(j,s2(2*i-1):s2(2*i)),span);              
    end   
    for j = 1:size(tpr,1)     
%         stpr(j,s1(2*i-1):s1(2*i)) = smooth(tpr(j,s1(2*i-1):s1(2*i)),span,'moving')';
%         stpr(j,s2(2*i-1):s2(2*i)) = smooth(tpr(j,s2(2*i-1):s2(2*i)),span,'moving')';
        stpr(j,s1(2*i-1):s1(2*i)) = medfilt1(tpr(j,s1(2*i-1):s1(2*i)),span);
        stpr(j,s2(2*i-1):s2(2*i)) = medfilt1(tpr(j,s2(2*i-1):s2(2*i)),span);
    end
end

% plot
close all
figure(1) 
plot(t,stpr(1,:),'-k','linewidth',0.5)
hold on
plot(t,stpr(2,:),'-k','linewidth',0.5)
hold off
grid
% label, legend
ylabel 'Temperature (K)'
text(0.04,0.917,'(a)','parent',gca,'units','normalized','fontsize',10);
text(0.172,0.87,'{\itz} = 0.5 mm','parent',gca,'units','normalized','fontsize',7);
text(0.172,0.51,'{\itz} = 1.5 mm','parent',gca,'units','normalized','fontsize',7);
text(0.412,0.89,'{\itI}_p = 150 A, {\itI}_b = 120 A','parent',gca,'units','normalized','fontsize',7);
% print figure
opts = struct('lbrt',[0.5 -149 0 0],'figsize',[15 9],'ticksize',[0.01 0 1],...
              'xlbl',[0.5 0.5 5.5],'ylbl',[16000 1000 21000],'axis',[.5 5.5 15500 21000]);
printfig('fig09',opts);
set(gca,'XTickLabel',[]);

% plot fig [1,2]
aps = fgx('apos');
aps1 = aps(1,:);
aps1(4) = 0.745*aps1(4);
aps1(2) = aps(1,2)-aps1(4)-0.01;
axes('position',aps1);
plot(t,sint(4,:),'-b','linewidth',0.5)
grid
% label, legend
ylabel 'Intensity (a.u.)'
text(0.04,0.84,'(b)','parent',gca,'units','normalized','fontsize',10); 
text(0.172,0.66,'{\itz} = -0.08 mm','parent',gca,'units','normalized','fontsize',7); 
opts = struct('ticksize',[0.01 0 1],...
              'xlbl',[0.5 0.5 5.5],'ylbl',[25 5 50],'axis',[.5 5.5 26 52]);
printfig('set',opts);
set(gca,'XTickLabel',[]);

% plot fig [1,3]
aps = fgx('apos');
aps1 = aps(1,:);
aps1(2) = aps(1,2)-aps1(4)-0.01;
axes('position',aps1);
plot(t,sint(1,:),'-r','linewidth',0.5)
grid
% label, legend
xlabel 'Time (s)'
ylabel 'Intensity (a.u.)'
text(0.04,0.86,'(c)','parent',gca,'units','normalized','fontsize',10); 
text(0.172,0.82,'{\itz} = -1.5 mm','parent',gca,'units','normalized','fontsize',7); 
opts = struct('ticksize',[0.01 0 1],...
              'xlbl',[0.5 0.5 5.5],'ylbl',[4 4 16],'axis',[.5 5.5 4 16]);
printfig('set',opts);
fgs('pfg','fig09'); % print fig
