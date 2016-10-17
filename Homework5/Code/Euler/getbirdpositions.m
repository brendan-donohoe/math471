function positions = getbirdpositions(xmin, xmax, ymin, ymax, nbirds, cx, cy, gamma1, gamma2, kappa, ro, lambda, delta, h, endt)
  numiterations = round(endt / h);
  positions = zeros(nbirds, 2, numiterations);

  % First, for each bird, generate a random position (x, y) within the given bounds.
  for i = 1:nbirds
    positions(i,:,1) = [((xmax - xmin) * rand) + xmin, ((ymax - ymin) * rand) + ymin];
  end
  
  t = 0;
  
  % Next, we begin getting the x and y value of each bird at each time t.
  for i = 1:numiterations
    birdcur = positions(:,:,i);
    birdnext = zeros(nbirds, 2);
    
    % First, compute the center of the current flock.
    centerx = sum(birdcur(:,1)) / nbirds;
    centery = sum(birdcur(:,2)) / nbirds;
    
    % Next, approximate the bird leader's next position at time t + h using
    % Runge Kutta.
    birdnext(1, 1) = rungekuttanext(@(t, y) leadereqn(t, y, cx, gamma1), t, birdcur(1, 1), h);
    birdnext(1, 2) = rungekuttanext(@(t, y) leadereqn(t, y, cy, gamma1), t, birdcur(1, 2), h);
    % Forward Euler.
    % birdnext(1, 1) = myforwardeuler(@(t, y) leadereqn(t, y, cx, gamma1), t, birdcur(1, 1), h);
    % birdnext(1, 2) = myforwardeuler(@(t, y) leadereqn(t, y, cy, gamma1), t, birdcur(1, 2), h);
    
    % Find the x and y coordinates at time t + h for the rest of the flock
    % as well.
    
    for k = 2:nbirds
      % First, for the current bird, find its lambda closest neighbors.
      nlist = birdcur(1:nbirds,:);
      nlist(k,:) = [];
      nindices = closestneighbors(birdcur(k,:), nlist, lambda);
    
      % Using Runge Kutta with the above parameters, approximate the bird's
      % next position.
      birdnext(k, 1) = rungekuttanext(@(t, y) birdeqn(t, y, gamma2, birdcur(1, 1), kappa, centerx, ro, delta, nlist(nindices, 1)), t, birdcur(k, 1), h);
      birdnext(k, 2) = rungekuttanext(@(t, y) birdeqn(t, y, gamma2, birdcur(1, 2), kappa, centery, ro, delta, nlist(nindices, 2)), t, birdcur(k, 2), h);
      % Forward Euler:
      % birdnext(k, 1) = myforwardeuler(@(t, y) birdeqn(t, y, gamma2, birdcur(1, 1), kappa, centerx, ro, delta, nlist(nindices, 1)), t, birdcur(k, 1), h);
      % birdnext(k, 2) = myforwardeuler(@(t, y) birdeqn(t, y, gamma2, birdcur(1, 2), kappa, centery, ro, delta, nlist(nindices, 2)), t, birdcur(k, 2), h);
    end
    
    % Store the information of birdnext in the positions array and
    % increment t.
    positions(:,:,i + 1) = birdnext;
    t = t + h;
  end
end