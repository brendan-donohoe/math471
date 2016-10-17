function frm = getbirdframe(xmin, xmax, ymin, ymax, birdxs, birdys, feedx, feedy)
  colors = zeros(length(birdxs) + 1, 3);
  colors(1,:) = [0, 255, 0];
  colors(2,:) = [255, 0, 0];
  for i = 3:(length(colors) - 1);
    colors(i,:) = [0, 0, 255];
  end
  colors(length(colors),:) = [0, 0, 0];
  figure('visible', 'off');
  scatter([feedx birdxs], [feedy birdys], 50, colors);
  axis([xmin, xmax, ymin, ymax]);
  frm = getframe;