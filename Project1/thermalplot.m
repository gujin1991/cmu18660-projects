function thermalplot(x)
pcolor(x');
shading interp;
caxis([27 50]);
colorbar;
xlabel('X Axis');
ylabel('Y Axis');