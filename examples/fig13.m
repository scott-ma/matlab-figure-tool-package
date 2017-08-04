% fig13

dt1 = 7;
dt2 = 10;
x0 = 200;
y0 = 205;
yd = 20;
ym = 300;

% load
syme_temp = load('fig2a');                     
syme_cell = struct2cell(syme_temp);
[temp_data] = deal(syme_cell{:});  

[j_d,i_d] = size(temp_data);
syme_data = zeros(j_d,i_d);
syme_data(1:end-15,:) = temp_data(16:end,:);
syme_data(end-19:end,:) = syme_data(end-19:end,:);
syme_data = syme_data(:,1:i_d/2+1);
syme_data(200:300,1:40) = 0;
% syme_data(200:300,end-40:end) = 0;

[j_temp, i_temp]=size(syme_data);
j_t = 0:1:j_temp-1;
i_t = -i_temp+1:0;
[X,Y]=meshgrid(i_t,j_t);
% close all
figure(1)
contourf(X,Y,syme_data,10);
h = findobj('Type','patch');
set(h,'LineWidth',1);
hold on

% plot emissivity
% load file
emis_temp = load('fig2b');                     
emis_cell = struct2cell(emis_temp);
[temp_data] = deal(emis_cell{:});
% [j_d,i_d] = size(temp_data);
% emis_data = zeros(j_d+2,i_d);
% emis_data(1:end-2,:) = temp_data(1:end,:);
emis_data = temp_data;

% emis_data = [fliplr(temp_data) temp_data(:,2:end)];
[j_temp, i_temp] = size(emis_data);
i_t = 0:i_temp-1;
% i_t(end) = 300;
j_t = 20:9.4:20+(j_temp-1)*9.4;
[X,Y]=meshgrid(i_t,j_t);
contourf(X,Y,emis_data,10);
h = findobj('Type','patch');
set(h,'LineWidth',1);

% plot cathode
xc = [-48 -48 0 48 48 -48];
yc = [ym y0+83 y0 y0+83 ym ym];
fill(xc,yc,'r');
plot(xc,yc);
b = y0;
xt = b-ym;
while xt<48
    if xt>-48  % P1
        x1 = xt;
        y1 = ym;
    else
        yt = b+48;
        if yt>=y0+83
            x1 = -48;
            y1 = yt;
        else
            x1 = 48/35*(y0-b);
            y1 = 48/35*(83/48*b-y0);
        end
    end
    yt = b-48; % P2
    if yt>=y0+83
        x2 = 48;
        y2 = yt;
    else
        x2 = 48/131*(b-y0);
        y2 = 48/131*(y0+83/48*b);
    end
% plot
 plot([x1 x2],[y1 y2]);
 b = b+dt1;
 xt = b-ym;
end

% plot anode
xa = [-x0 -x0 x0 x0 -x0];
ya = [yd 0 0 yd yd];
fill(xa,ya,[0.7 0.6 0.5]);
plot(xa,ya);
b = -x0;
xt = b-yd;
while xt<x0
    yt = b+x0; % P1
    if yt<yd
        x1 = -x0;
        y1 = yt;
    else
        x1 = b-yd;
        y1 = yd;
    end
    if b<x0 % P2
        x2 = b;
        y2 = 0;
    else
        x2 = x0;
        y2 = b-x0;
    end
    % plot
    plot([x1 x2],[y1 y2]);
    b = b+dt2; 
    xt = b-yd;
end
plot([0 0],[yd y0]);
hold off
axis([-200 200 0 300]);

% print figure 
opts = struct('lbrt',[-1 0 0 0],'figsize',[8.9 6/8*8.9],'tickxylblfs',[8 9 9],'ticksize',[0.015 0.025 0], ...
              'xlbl',[-200 400 200],'ylbl',[-50 400 400],  'figcolor','c');
outmyfig('fig13',opts);