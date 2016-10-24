file = fopen('parallel.txt','r');
format = '%f %f %f';
size = [16 781 2];
Pdata = fscanf(file, format, size);
fclose(file);

thread = 1;
while thread <= 16
    gridSize = Pdata(thread,:,1);
    time = Pdata(thread,:,2);
    
    figure('position', [0, 0, 800, 800]);
    plot(gridSize, time, 'r');
    grid on;
    xlabel('Grid size', 'FontSize', 12);
    ylabel('Time', 'FontSize', 12);
    legend('Time for grid size with ' + thread + ' threads', 'h^2');
    saveas(gcf, 'parallel' + threads + '.png');
    thread = thread + 1;
end