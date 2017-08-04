% fig03

l1 = 3;
l2 = 5;
r1 = 1;  % central pixel
r2 = 1;  % mm
pm = 12;   % pixels/mm
% load data
temp_file = 'G:\eProc\timevary\TAI200s1Ar794';  % start
temp_temp = load(temp_file);                     
temp_cell = struct2cell(temp_temp);
[temp_data] = deal(temp_cell{:});  
temp_data = double(temp_data)/1000;
da1s = squeeze(temp_data(l1,:,r1));
da2s = squeeze(temp_data(l2,:,r1));
da3s = squeeze(temp_data(l2,:,r2*3*pm+1));
clear temp_temp temp_cell temp_data
temp_file = 'G:\eProc\timevary\TAI200e1Ar794';  % end
temp_temp = load(temp_file);                     
temp_cell = struct2cell(temp_temp);
[temp_data] = deal(temp_cell{:});  
temp_data = double(temp_data)/1000;
da1e = squeeze(temp_data(l1,:,r1));
da2e = squeeze(temp_data(l2,:,r1));
da3e = squeeze(temp_data(l2,:,r2*3*pm+1));
clear temp_temp temp_cell temp_data

% fig
xs = (1:length(da1s))/955;
xe = 6.5+(1:length(da1e))/955;
% plot 
close all
figure(1)
fsz = [8.3 5.8]; 
rxy = [0.45 1-0.45];
xyl = {'Time (s)','Temperature (1000 K)'};
xyt = [0 2 4 -0.5 5.5; 30 2 36 29.5 36.1; 14 1 20 13.5 20.5];
tks = [0.02 0 1];
fgx('fgbx','a',fsz,rxy);
plot(xs,da1s,'-k',xs,da2s,'--b',xs,da3s,'-.r')
text(0.98,0.61,'a','fontsize',8,'units','normalized');
text(0.98,0.40,'c','fontsize',8,'units','normalized');
text(0.98,0.15,'e','fontsize',8,'units','normalized');
fgx('fgbx','ab',xyt,tks,rxy,1);
plot(xe+24,da1e,'-k',xe+24,da2e,'--b',xe+24,da3e,'-.r')
text(0.09,0.50,'b','fontsize',8,'units','normalized');
text(0.09,0.33,'d','fontsize',8,'units','normalized');
text(0.09,0.12,'f','fontsize',8,'units','normalized');
hl = legend('a,b: z = 1 mm, r = 0 mm','c,d: z = 2 mm, r = 0 mm','e, f: z = 2 mm, r = 1 mm',0); 
fgx('fgbx','b',xyt,tks,xyl);
% print fig
printfig('fig03','lbrt',[0 0 4.5 1],'figsize',fsz,'bkxy','x','legend',[hl,0.57,0.63,8]);
