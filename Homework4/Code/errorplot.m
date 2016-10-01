fid = ('errorvals.txt', 'r');
formatSpec = '%f %f';
sizeData = [2 Inf];
data = fscanf(fid, formatSpec, sizeData);
fclose(fid);

effectiveh = data(1,:);
error = data(2,:);

figure('position', [0, 0, 800, 800]);
loglog(effectiveh, error, 'r', effectiveh, effectiveh.^2, '--b');
grid on;
xlabel('effective h', 'FontSize', 12);
ylabel('max error', 'FontSize', 12);
legend('max error', 'h^2');
saveas(gcf, 'errorgraph.png');
