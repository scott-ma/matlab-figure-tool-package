%fig11

i_s = 31;  % w_lef
j_s = 10;  % h_ups
i_n = 160; % width
j_n = 120; % hight
gap = 2;
ss = [167 156 159 163 582 585 591 593]-1;
fg_name = 'E:\HighSpeed\ac081011\sa300\sa300_';
% load image 
dat = ones(2*j_n+gap,4*i_n+3*gap);
dat = im2uint8(dat);

for k = 1:8
    i = ss(k);
    if i < 10
        fg_file = strcat(fg_name,'00',num2str(i),'.tif');  
    elseif i <100
        fg_file = strcat(fg_name,'0', num2str(i),'.tif');   
    else
        fg_file = strcat(fg_name,     num2str(i),'.tif');  
    end
    infig = imread(fg_file);
    infig = 1.13*double(infig)-4;
    otfig = zeros(j_n,i_n);
    otfig = infig(j_s:j_s+j_n-1,i_s:i_s+i_n-1);
    if k<5
        m = 1;
        n = (k-1)*(i_n+gap)+1;
    else
        m = j_n+gap+1;
        n = (k-5)*(i_n+gap)+1;
    end
    if k==2 | k==6
        temp = uint8(round((otfig-4)*245/141+4));
    else
        temp = uint8(otfig);
    end
    dat(m:m+j_n-1,n:n+i_n-1) = temp;
end
[mi mj] = size(dat);

% plot
close all
fgn(1)
imshow(dat)
dx = -5;
text(dx+15,15,'SP','fontsize',7,'color',[1 1 1]); 
text(dx+175,15,'SP\rightarrowRP','fontsize',7,'color',[1 1 1]); 
text(dx+340,15,'RP','fontsize',7,'color',[1 1 1]); 
text(dx+500,15,'RP\rightarrowSP','fontsize',7,'color',[1 1 1]); 
dy = 120;
text(dx+15,15+dy,'SP','fontsize',7,'color',[1 1 1])
text(dx+175,15+dy,'SP\rightarrowRP','fontsize',7,'color',[1 1 1])
text(dx+340,15+dy,'RP','fontsize',7,'color',[1 1 1])
text(dx+500,15+dy,'RP\rightarrowSP','fontsize',7,'color',[1 1 1])
% print figure
opts = struct('lrdu',[0 0 0 0],'figsize',[12 12*mi/mj],...
              'xlbl',[0 mj+10 mj],'ylbl',[0 mi+10 mi]);
outmyfig(gcf,'fig11',opts);