file = fopen('timetable.txt','r');
format = '%d %d %f';
size = [3 Inf];
pdata = fscanf(file, format, size);
fclose(file);

tablematrix = zeros(12, 3, 781);

startIdx = 1;
endIdx = 781;

for i = 1:16
  tablematrix(i, :, :) = pdata(:, startIdx:endIdx);
  startIdx = startIdx + 781;
  endIdx = endIdx + 781;
end

thread = 1;
while thread <= 16
  gridsize = reshape(tablematrix(thread,2,:), [1, 781]);
  time = reshape(tablematrix(thread,3,:), [1, 781]);
  scaledtime = time ./ (gridsize .^ 2);
  
  figure('position', [0, 0, 800, 800]);
  plot(gridsize, scaledtime, 'r');
  grid on;
  xlabel('Grid Size', 'FontSize', 12);
  ylabel('Scaled Wall Clock Time', 'FontSize', 12);
  filename = strcat('parallel', num2str(thread), '.png');
  saveas(gcf, filename);
  thread = thread + 1;
end