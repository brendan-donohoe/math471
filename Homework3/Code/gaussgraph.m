fid = fopen('gaussQuad.txt','r');
formatSpec = '%d %d %f';
sizeData = [3 Inf];
data = fscanf(fid, formatSpec, sizeData);
fclose(fid);

data1 = data(:, 1:length(data) / 2);
data2 = data(:, length(data) / 2 + 1:end);

ns = data1(2, :);

ts1 = data1(3, :);
ts2 = data2(3, :);

best1 = ts1(length(ts1));
best2 = ts2(length(ts2));

e1 = abs(ts1 - best1);
e2 = abs(ts2 - best2);
e_approx = (1.2).^(-(4).*(ns));

figure('position', [0, 0, 800, 800]);
loglog(ns, e1, 'r', ns, e2, 'b', ns, e_approx, 'g');
%loglog(ns, e1, 'r', ns, e2, 'b');
ylim([10.^(-16),100]);
grid on;
legend('e^{cos(\pi*x)}','e^{cos(\pi^2*x)}','(1.2)^{(-4*x)}');
%legend('e^{cos(\pi*x)}','e^{cos(\pi^2*x)}');
xlabel('n', 'FontSize', 12);
ylabel('ERROR', 'FontSize', 12);
saveas(gcf,'gausserrgraph2.png');