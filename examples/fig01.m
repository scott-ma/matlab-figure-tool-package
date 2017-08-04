% fig 01
close all

n = 20;
dw= pi;
% calculate values of matrix H & F
i=0:n;
j=0:n;
[k_vh,j_v]=meshgrid(i,j);   % for matrix j,k,i
i_v=k_vh;
k_vf=j_v;
H = besselj(0,j_v.*k_vh*dw/n); % calculate H
F = k_vf.*cos(k_vf.*i_v*dw/n); % calculate F
F(:,1) = F(:,1)/2;

idat = pfm(n,20); 
Gk = F*idat(1,:)'/n; 
dg = H-gjkm(n,dw);

k = 1:n;
een = dg(n,:).*Gk';
figure(1); 
plot(k,Gk(2:end)','-k',k,dg(1,2:end),'-.b',k,een(2:end)','-or','MarkerFaceColor',[1 1 1],'MarkerSize',3);

% label, legend
xlabel 'k'
ylabel 'kG(\alphak)/n, f, kG(\alphak)f/n'
text(18,0.215,'(b)','fontsize',10);
text(10,0.21,'n = 20, \alpha = 1','fontsize',8);
hl = legend('kG(\alphak)/n','f(k\pi, \alphar_n)','kG(\alphak)f(k\pi, \alphar_n)/n',0);

% print figure 
opts = struct('lbrt',[-1 2 1 2],'figsize',[8.5 6],'tickxylblfs',[8 9],'ticksize',[0.015 0.025 0], ...
              'xlbl',[0 5 20],'ylbl',[-0.2 0.1 0.2],'axis',[0 20 -0.25 0.25], ...
              'figcolor','k','legend',[hl,0.49,0.56,8],'relegend',[3 1]);
printfig('fig01',opts);