% fig07
file = 'G:\eProc\timevary\surface\'; % 50-18:50
fint = strcat(file,'CI200upAr783');
fps  = 955;
colc = 84;
line = [43 48]-31;
colm = [87 85];

% load ctuum int
orig_temp = load(fint);                     
orig_cell = struct2cell(orig_temp);
[orig_data] = deal(orig_cell{:});  
int = double(orig_data);  % uint8-->double
clear orig_temp orig_cell orig_data
ttt = (0:size(int,2)-1)/fps+1;

int1c = squeeze(int(line(1),:,colc));
int1a = squeeze(int(line(1),:,colm(1)));
int1b = squeeze(int(line(1),:,colm(1)+1));
int2c = squeeze(int(line(2),:,colc));
int2a = squeeze(int(line(2),:,colm(2)));
int2b = squeeze(int(line(2),:,colm(2)+1));

% plot
close all
figure(1) 
plot(ttt,int2c,'-k',ttt,int2a,'-r',ttt,int2b,'-b',ttt,int2a-0.5*int2b,'-g')
% label, legend
xlabel '{\itt} (s)'
ylabel 'Intensity (a.u.)'
text(1.9,152,'{\itz} = -0.09 mm','fontsize',7);
text(1.9,138,'p_1','fontsize',7);
text(3.8,115,'p_2','fontsize',7);
text(1.9,100,'p_3','fontsize',7);
text(2.7,53,'{\itI}_2-{\itI}_3/2','fontsize',7);
text(0.03,0.96,'(a)','parent',gca,'units','normalized','fontsize',9);
% print figure
opts = struct('lbrt',[0 0 -113 0.4],'figsize',[8.5 5.8],'ticksize',[0.02 0 1],'tickxylblfs',[7 8],...
              'xlbl',[1 1 4.5],'ylbl',[40 20 160],'axis',[1 4.5 40 160]);
printfig('fig07',opts);

% plot fig [1,2]
aps = fgx('apos');
aps1 = aps(1,:);
aps1(1) = aps(1,1)+aps(1,3)+0.07;
axes('position',aps1);
plot(ttt,int1c,'-k',ttt,int1a,'-r',ttt,int1b,'-b',ttt,int1a-0.5*int1b,'-g')
% label, legend
xlabel '{\itt} (s)'
text(1.9,191,'{\itz} = -0.55 mm','fontsize',7);
text(2.2,170,'p_4','fontsize',7);
text(2.2,117,'p_5','fontsize',7);
text(2.2,85,'p_6','fontsize',7);
text(2.7,57,'{\itI}_5-{\itI}_6/2','fontsize',7);
text(0.03,0.96,'(b)','parent',gca,'units','normalized','fontsize',9); 
opts = struct('ticksize',[0.02 0 1],'tickxylblfs',[7 8],...
              'xlbl',[1 1 4.5],'ylbl',[50 25 200],'axis',[1 4.5 50 200]);
printfig('set',opts);
fgs('pfg','fig07'); % print fig
