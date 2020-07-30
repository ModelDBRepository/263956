raster = load('raster.dat');

figure; 
plot(raster(:,2:end),raster(:,1),'k.','MarkerSize',0.5)
xlabel('Time [ms]')
ylabel('Neuron #')
saveas(gcf,'raster.png')