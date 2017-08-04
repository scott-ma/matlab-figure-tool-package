% fig02

k  = 1.3806580*10^(-23);      % J/K
h  = 6.62606896*10^(-34);     % J s
e  = 1.6021892*10^(-19);      % C
c  = 299792458;               % m s^-1
C1 = 1.632*10^(-43);          % J m^4 K^1/2  s^-1 sr^-1
C2 = 1.026*10^(-34);          % J m^2 K^-3/2 s^-1 sr^-1

lam= 100:20:1000;             % nm

ff1= 1.23;
ff2= 1.0;

% load
Arnu = proc('load','Ar_T_n0e_u02');
Te = Arnu(:,1); u0 = Arnu(:,6); u1 = Arnu(:,7); u2 = Arnu(:,8);
n0 = Arnu(:,2); n1 = Arnu(:,3); n2 = Arnu(:,4); ne = Arnu(:,5);

QT = load('ArQT','-ascii');

BII = load('Arfb1.mat','-ascii');     
BIII= load('Arfb2.mat','-ascii'); 

fb1a = interp1(BII(1,:),BII(102,:),lam,'pchip');
fb1b = interp1(BII(1,:),BII(202,:),lam,'pchip');
fb1c = interp1(BII(1,:),BII(302,:),lam,'pchip');
fb1d = interp1(BII(1,:),BII(402,:),lam,'pchip');

fb2a = interp1(BIII(1,:),BIII(102,:),lam,'pchip');
fb2b = interp1(BIII(1,:),BIII(202,:),lam,'pchip');
fb2c = interp1(BIII(1,:),BIII(302,:),lam,'pchip');
fb2d = interp1(BIII(1,:),BIII(402,:),lam,'pchip');

% Biberman
lam  = lam*1e-9;  % nm->m
xi1a = (1-exp(-h*c./(lam*k*Te(101)))).*fb1a + exp(-h*c./(lam*k*Te(101)))*ff1;
xi1b = (1-exp(-h*c./(lam*k*Te(201)))).*fb1b + exp(-h*c./(lam*k*Te(201)))*ff1;
xi1c = (1-exp(-h*c./(lam*k*Te(301)))).*fb1c + exp(-h*c./(lam*k*Te(301)))*ff1;
xi1d = (1-exp(-h*c./(lam*k*Te(401)))).*fb1d + exp(-h*c./(lam*k*Te(401)))*ff1;
xi2a = (1-exp(-h*c./(lam*k*Te(101)))).*fb2a + exp(-h*c./(lam*k*Te(101)))*ff2;
xi2b = (1-exp(-h*c./(lam*k*Te(201)))).*fb2b + exp(-h*c./(lam*k*Te(201)))*ff2;
xi2c = (1-exp(-h*c./(lam*k*Te(301)))).*fb2c + exp(-h*c./(lam*k*Te(301)))*ff2;
xi2d = (1-exp(-h*c./(lam*k*Te(401)))).*fb2d + exp(-h*c./(lam*k*Te(401)))*ff2;

factor = 2.0;
xi2af = factor*(1-exp(-h*c./(lam*k*Te(101)))).*fb2a + exp(-h*c./(lam*k*Te(101)))*ff2;
xi2bf = factor*(1-exp(-h*c./(lam*k*Te(201)))).*fb2b + exp(-h*c./(lam*k*Te(201)))*ff2;
xi2cf = factor*(1-exp(-h*c./(lam*k*Te(301)))).*fb2c + exp(-h*c./(lam*k*Te(301)))*ff2;
xi2df = factor*(1-exp(-h*c./(lam*k*Te(401)))).*fb2d + exp(-h*c./(lam*k*Te(401)))*ff2;

% load data
[d1,d2] = textread('f6_1998_jqsrt.txt','%f%f','headerlines',2);
l1 = d1(1:round(d1(1))-2);
l2 = d1(round(d1(1))-1:round(d2(1))-3);
x1 = d2(1:round(d1(1))-2);
x2 = d2(round(d1(1))-1:round(d2(1))-3);

[d1,d2] = textread('f4_1999_jpd.txt','%f%f','headerlines',2);
d1 = 468.8*ones(length(d1));

% plot
lam = lam*1e9; 
close all
fgx(1)
plot(lam,xi1a,':r', 'LineWidth',0.5)
hold on
plot(lam,xi1b,'-.b','LineWidth',0.5)
plot(lam,xi1c,'--m','LineWidth',0.5)
plot(lam,xi1d,'-g', 'LineWidth',0.5)
plot(l1,x1,'-k','LineWidth',0.8)
plot([0,1],[-1,-1],'w')
plot(l2,x2,'--k','LineWidth',0.8)
plot([0,1],[-1,-1],'w')
hold off
% label, legend
xlabel 'Wavelength (nm)'     
ylabel 'Total Biberman factor \xi_1'
hl = legend('{\fontsize {7.5} 10000 K}','{\fontsize {7.5} 15000 K}',...
    '{\fontsize {7.5} 20000 K}','{\fontsize {7.5} 25000 K}',...
    '{\fontsize {7.5} D''yachkov et al}','{\fontsize {7.5} 12000 K, n_e = 10^{16} cm^{-3}}',...
    '{\fontsize {7.5} D''yachkov et al}','{\fontsize {7.5} 12000 K, n_e = 10^{17} cm^{-3}}',0); 
tl = text(0.90,0.93,'(a)','parent',gca,'units','normalized','fontsize',12);
tl = text(0.2,0.34,'Hofsaess','parent',gca,'units','normalized','fontsize',7);
% print figure
opts = struct('lbrt',[0 0 -4 0],'figsize',[8.3 6],'ticksize',[0.02 0.025 1], ...
              'xlbl',[200 100 1000], 'ylbl',[0 0.4 2.4],'axis',[200 1000 0 2.4], ...
              'legend',[hl 0.24 0.13],'relegend',[4 2 0.3]);
printfig('fig02',opts)