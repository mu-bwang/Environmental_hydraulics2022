clear all
close all
clc

lw = 1.5;
ms = 10;

dxy = 0.0023;
v = xlsread('jet_dataset.xlsx');
[m,n] = size(v);

% define physical scale 
xPIV = [1:n] .* dxy;
zPIV = [1:m] .* dxy;

% define jet half-width at all heights, give zeros fow now.
bg = zeros(1,m);

%% Plot velocity map
figure
imagesc(xPIV,zPIV,v);
title('Vertical velocity (m/s)')
colorbar
xlabel('xPIV (m)')
ylabel('zPIV (m)')
set(gca,'fontsize',18)


%% Fit Gaussian distribution to velocity profile at each height 
for i=1:m
    vz = v(i,:);
    I = find(vz == max(vz));
    r = xPIV - xPIV(I);
    figure(9)
    plot(r, vz);
    hold on
    f = @(x, p) p(1) * exp(-(x/p(2)).^2);
    p = lsqnonlin(@(p) f(r, p) - vz, [max(vz) 0.2]);
    bg(i) = abs(p(2));
    plot(r, f(r,p),'r')
    hold off    
end


%% Plot fitted b_g as a function of relative height
figure(2)
plot(bg,zPIV,'o','LineWidth',lw,'MarkerSize',ms)
set(gca,'fontsize',18)
xlabel('b_g (m)')
ylabel('zPIV (m)')

%% Determine the origin of the jet, so that we can convert relative height to absolute height above the nozzle
f = fittype('x+z0'); 
fitresult = fit(zPIV',bg'/0.11,f);
z_origin = fitresult.z0;

jet_height = linspace(0,0.2);
figure(3)
hold on
box on
grid on
plot(bg,zPIV+z_origin,'o','LineWidth',lw,'MarkerSize',ms)
plot(0.11*jet_height,jet_height, 'r','LineWidth',lw);
xlabel('b_g (m)')
ylabel('z (m)')
tt = sprintf('z_{origin} = %3.2f m', z_origin);
text(0.008,0.06,tt,'FontSize',18);
axis([0 0.02 0 0.2])
h = legend('Data','Theory');
set(h,'location','northwest')
set(gca,'fontsize',18)

