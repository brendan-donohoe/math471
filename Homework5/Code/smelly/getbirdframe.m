function frm = getbirdframe(xmin, xmax, ymin, ymax, birdxs, birdys, feedx, feedy)
  colors = zeros(length(birdxs) + 1, 3);
  colors(1,:) = [0, 255, 0];
  colors(2,:) = [255, 0, 0];
  colors(3,:) = [0, 0, 0];
  for i = 4:length(colors);
    colors(i,:) = [0, 0, 255];
  end
  figure('visible', 'off');
  scatter([feedx birdxs], [feedy birdys], 50, colors);
  axis([xmin, xmax, ymin, ymax]);
  frm = getframe;