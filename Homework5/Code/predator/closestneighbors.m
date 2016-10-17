function neighborindices = closestneighbors(dead, location, neighbors, n)
  dist = @(x1, y1, x2, y2) (sqrt((x2 - x1)^2 + (y2 - y1)^2));
  disttable = zeros(length(neighbors), 2);
  x1 = location(1);
  y1 = location(2);
  for i = 1:length(neighbors(:,1))
    if neighbors (i, 1) == dead
        disttable(i, :) = [dist(x1, y1, 100, 100), i];
    else
        x2 = neighbors(i, 1);
        y2 = neighbors(i, 2);
        disttable(i, :) = [dist(x1, y1, x2, y2), i];
    end
  end
  disttable = sortrows(disttable, 1);
  neighborindices = disttable(1:n, 2);
end