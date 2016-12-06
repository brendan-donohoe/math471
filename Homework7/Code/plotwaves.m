function plotwaves(t, numProcs)
file = fopen(strcat('t',int2str(t),'up0.txt'), 'r');
formatSpec = '%f %f %f';
sizeData = [3 Inf];
tempdata = fscanf(file, formatSpec, sizeData);
data = tempdata;
fclose(file);
i = 1;
while i < numProcs
    file = fopen(strcat('t',int2str(t),'up',int2str(i),'.txt'), 'r');
    tempdata = fscanf(file, formatSpec, sizeData);
    data = [data,tempdata];
    fclose(file);
    i = i + 1;
end

x = data(1,:);
y = data(2,:);
u = data(3,:);


figure('visible','off');
plot3(x,y,u);
print(['frame' num2str(t,'%5.5i') '.jpg'],'-djpeg100');
end
