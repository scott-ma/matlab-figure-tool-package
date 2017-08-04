% fig05

sp = [145 525 910];  % peak series
% sb = [ 50 435 810];  % base series
sb = [240 620 994];  % base series
cc = 88;             % column
fg_name = 'G:\eData1\pulsdc\080320\f005c\f005d5c_';

% load image & get
lp = length(sp);
lb = length(sb);
av = 4; % half average value
da = zeros(lp,150); % peak
for k = 1:lp
    for i = sp(k)-av:sp(k)+av
        if i < 10
            fg_file = strcat(fg_name,'00',num2str(i),'.tif');  
        elseif i <100
            fg_file = strcat(fg_name,'0', num2str(i),'.tif');  ;  
        else
            fg_file = strcat(fg_name,     num2str(i),'.tif');  
        end
        infig = imread(fg_file);
        infig = double(infig);
        % drop bkgrd
%         bakx = infig(1:5,1:5);
%         [bxi,bxj] = size(bakx);
%         back = sum(sum(bakx))/(bxi*bxj);        
        da(k,:) = da(k,:)+infig(:,cc)';
    end
    da(k,:) = da(k,:)/(2*av+1)-5; % drop bkgrd
end
db = zeros(lb,150); % base 
for k = 1:lb
    for i = sb(k)-av:sb(k)+av
        if i < 10
            fg_file = strcat(fg_name,'00',num2str(i),'.tif');  
        elseif i <100
            fg_file = strcat(fg_name,'0', num2str(i),'.tif'); 
        else
            fg_file = strcat(fg_name,     num2str(i),'.tif');  
        end
        infig = imread(fg_file);
        infig = double(infig);
        db(k,:) = db(k,:)+infig(:,cc)';
    end
    db(k,:) = db(k,:)/(2*av+1)-5;  % drop bkgrd
end
zz = (1-50:150-50)/12;

% plot
close all
fgx(1)
plot(zz,da(1,:),'-k',zz,da(2,:),'--b',zz,da(3,:),':r',...
     zz,db(1,:),'-k',zz,db(2,:),'--b',zz,db(3,:),':r','linewidth',0.7)
% label, legend
xlabel 'Distance from cathode tip (mm)'
ylabel 'Intensity (a.u.)'
text(0.56,0.88,'Peak','parent',gca,'units','normalized','fontsize',7);
text(0.65,0.61,'Base','parent',gca,'units','normalized','fontsize',7);
text(0.74,0.27,'Peak  Base','parent',gca,'units','normalized','fontsize',7);
text(0.63,0.36,'{\itI}_p = 200 A, {\itI}_b = 100 A','parent',gca,'units','normalized','fontsize',7);
text(0.9,0.94,'(b)','parent',gca,'units','normalized','fontsize',10);
hl = legend('0.55, 0.65 s','0.95, 1.05 s','1.35, 1.45 s',0);
% print figure
opts = struct('lbrt',[0 0 1 2],'figsize',[8.3 6.3],'ticksize',[0.02 0.4 1], ...
              'xlbl',[-2 .5 .5],'ylbl',[0 50 200],'axis',[-2 .5 0 210], ...
              'legend',[hl,0.6,0.1,8],'relegend',[3 1]);
printfig('fig05',opts);

% plot subfig
dat = load('sfcz_dat.dat','-ascii');
tt = (0:size(dat,2)-1)/955+0.002+0.4;
aps = [0.22 0.685 0.335 0.27];
axes('position',aps);
plot(tt,dat(2,:)-5,'-b','linewidth',0.5); grid
% label
xlabel 'Time (s)'
ylabel 'Intensity (a.u.)'
text(0.06,0.85,'{\itz} = -2.0 mm','parent',gca,'units','normalized','fontsize',7);
% set
opts = struct('ticksize',[0.03 0 1],'tickxylblfs',[7 7.5],...
              'xlbl',[0.5 0.2 1.5],'ylbl',[10 5 25],'axis',[0.5 1.5 10 25]);
printfig('set',opts);
fgs('lxy',[0 2 4 0]); % xy label position
fgs('pfg','fig05')     % print fig
 