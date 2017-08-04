% fig04

file = 'E:\experProc\pulsdc\080916\';        
feqs = {'ef02','ef05','ef10','ef20','ef50','ef1h'};
lf = length(feqs);
ll = [1:7 100];
lt = length(ll)-1;
dat = zeros(lf,lt,999);
for i = 1:lf
    orig_file = strcat(file,'A',feqs{i},'Ar794');
    % load matrix
    orig_temp = load(orig_file);                     
    orig_cell = struct2cell(orig_temp);
    [orig_data] = deal(orig_cell{:});  
    orig_data = double(orig_data);   % uint8-->double
    clear orig_file orig_temp orig_cell
    dat(i,:,:) = orig_data(ll(1:end-1),:,ll(end));
end

% plot
tt = (1:999)/955-0.025;
clor = {'k','r','b','m','r','b'};
close all
figure(1)
k = 2;
sep = [0 24 48 62 91 115];
da = squeeze(dat(:,k,:));
hold on
for i = 1:lf
    dd = squeeze(da(i,:));
    plot(tt,dd+sep(i),clor{i})
end
hold off
box on
% label, legend
xlabel 'Time (s)'
ylabel 'Intensity (a.u.)'
text(.605,157,'{\itf} = 2 Hz','fontsize',7);
text(.675,180,'5 Hz','fontsize',7);
text(.665,202,'10 Hz','fontsize',7);
text(.665,226,'20 Hz','fontsize',7);
text(.665,250,'50 Hz','fontsize',7);
text(.645,270,'100 Hz','fontsize',7);
text(.1,228,'{\itI}_p = 200 A, {\itI}_b = 150 A','fontsize',7);
text(.2,251,'{\itz} = 0.5 mm','fontsize',7);
% print figure
opts = struct('lrdu',[7 1 9 3],'figsize',[7 5],'ticksize',[0.02 0.025 1], ...
              'xlbl',[0 .2 1],'ylbl',[100 350 800],'axis',[-0.04 1.04 130 270]);
printfig(gcf,'fig04',opts);