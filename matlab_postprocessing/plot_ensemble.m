% read and plot the 'new' ensembles from THB within HC.
clear;
close all;

% ensemble_file = '../ensemble_semucb_l2-7.txt';
ensemble_file = '../test_mpi.txt';
result = load_ensemble(ensemble_file);
disp('done loading ensemble');
%%
figure;
subplot(2,2,1);
hist(result.residual,50);
title('residuals');

subplot(2,2,2);
hist(log10(result.var),20);
title('log of var');

subplot(2,2,3);
hist(result.nlayer,[0:1:17]);
title('number of layers');

subplot(2,2,4);
plot(result.residual);
hold on
plot(result.var);
title('residual')

%%
nr = 100;
rmin = result.rad(1,1);
rmax = 1.0;
r = linspace(rmin,rmax,nr);
visc = zeros(nr,result.n);
rr = zeros(size(visc));

for i=1:result.n
   visc(:,i) = interp1( result.rad(1:result.nlayer(i),i), result.visc(1:result.nlayer(i),i),r);
   rr(:,i) = r;
end

%%
[N,c] = hist3([visc(:),rr(:)],'Nbins',[200 100]);
figure;
surf(c{1},c{2},N'); shading interp;
colorbar();
title('non-aligned solutions');

%%
visca = zeros(size(visc));
for i=1:result.n
    visca(:,i) = visc(:,i) - mean(visc(:,i));
    rr(:,i) = r;
end

% figure, imagesc(visc); colorbar;

[N,c] = hist3([visca(:),rr(:)],'Nbins',[110 100]);
figure;
pcolor(c{1},c{2},N'); shading interp;
colorbar();
title('aligned solutions');

lmean = mean(visca,2);
hold on
plot(lmean,r,'r');