% fig10

clear all
% load
infig = imread('fig1.bmp');
infig = im2double(infig);
infig = infig-60/256;
maxf = zeros(1,3);
minf = zeros(1,3);
for i = 1:3
    maxf(i) = max(max(infig(:,:,i)));
    minf(i) = min(min(infig(:,:,i)));
end
for i = 1:3
    ttfig(:,:,i) = (infig(:,:,i)-minf(i)+5/256)/(maxf(i)-minf(i)+5/256);
    otfig(:,:,i) = ttfig(:,:,i);
    otfig(11:end,:,i) = ttfig(1:end-10,:,i);
end
otfig = im2uint8(otfig);
for i = 1:3
    oofig(:,:,i) = otfig(21:end-10,:,i);
end

ttfig = sum(oofig,3)/3;
ttfig = uint8(ttfig);
% ttfig = medfilt2(ttfig,[5 5]);
ttfig = wiener2(ttfig,[3 3]); 

close all
figure(1);
imshow(ttfig);

% print figure
opts = struct('lbrt',[0 0 0 0.5],'figsize',[8.9 5.4/8*8.9],'ticksize',[0.01 0.025 0], ...
              'xlbl',[-10 500 400],'ylbl',[-10 400 300],'figcolor','c');
outmyfig('fig10',opts);