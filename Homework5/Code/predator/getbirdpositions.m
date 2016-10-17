function positions = getbirdpositions(xmin, xmax, ymin, ymax, nbirds, cx, cy, gamma1, gamma2, kappa, ro, lambda, delta, h, endt, pursuit)
  numiterations = round(endt / h);
  positions = zeros(nbirds + 1, 2, numiterations);
  eaten = 0;
  lambda2 = lambda;
  if xmin < ymin
      dead = -100;
  else
      dead = -100;
  end
  
  % First, for each bird, generate a random position (x, y) within the given bounds.
  for i = 1:(nbirds + 1)
    positions(i,:,1) = [((xmax - xmin) * rand) + xmin, ((ymax - ymin) * rand) + ymin];
  end
  
  t = 0;
  
  % Next, we begin getting the x and y value of each bird at each time t.
  for i = 1:numiterations
    birdcur = positions(:,:,i);
    birdnext = zeros(nbirds + 1, 2);
    
    % Decrease lambda, if necessary.
    if lambda2 > (nbirds - eaten)
        lambda2 = nbirds - eaten;
    end
        
    % First, compute the center of the current flock.
    % centerx = sum(birdcur(:,1)) / nbirds;
    % centery = sum(birdcur(:,2)) / nbirds;
    centerx = livingSum(dead, nbirds, birdcur, 1);
    centery = livingSum(dead, nbirds, birdcur, 2);
    
    % Next, approximate the bird leader's next position at time t + h using
    % Runge Kutta.
    birdnext(1, 1) = rungekuttanext(@(t, y) leadereqn(t, y, cx, gamma1), t, birdcur(1, 1), h);
    birdnext(1, 2) = rungekuttanext(@(t, y) leadereqn(t, y, cy, gamma1), t, birdcur(1, 2), h);
    
    % Find the x and y coordinates at time t + h for the rest of the flock
    % as well.
    predlist = birdcur(1:nbirds,:);
    for k = 2:nbirds
      if birdcur(k, 1) == dead
          birdnext(k, 1) = dead;
          birdnext(k, 2) = dead;
          %predlist(k,:) = [];
      else
          % First, for the current bird, find its lambda closest neighbors.
          nlist = birdcur(1:nbirds,:);
          nlist(k,:) = [];
          nindices = closestneighbors(dead, birdcur(k,:), nlist, lambda2);
    
          % Using Runge Kutta with the above parameters, approximate the bird's
          % next position.
          birdnext(k, 1) = rungekuttanext(@(t, y) birdeqn(t, y, gamma2, birdcur(1, 1), kappa, centerx, ro, delta, nlist(nindices, 1)), t, birdcur(k, 1), h);
          birdnext(k, 2) = rungekuttanext(@(t, y) birdeqn(t, y, gamma2, birdcur(1, 2), kappa, centery, ro, delta, nlist(nindices, 2)), t, birdcur(k, 2), h);
      end
    end
    
    % Update the predator.
    closestindices = closestneighbors(dead, birdcur(nbirds + 1,:), predlist, 1);
    birdnext(nbirds + 1, 1) = rungekuttanext(@(t, y) predeqn(t, y, pursuit, birdcur(closestindices, 1), eaten), t, birdcur(nbirds + 1, 1), h);
    birdnext(nbirds + 1, 2) = rungekuttanext(@(t, y) predeqn(t, y, pursuit, birdcur(closestindices, 2), eaten), t, birdcur(nbirds + 1, 2), h);
    
    % Check for dead birdies.
    for j = 2:nbirds
        if abs(birdcur(j, 1) - birdcur(nbirds + 1, 1)) < 0.05 && abs(birdcur(j, 2) - birdcur(nbirds + 1, 2)) < 0.05
            birdnext(j, 1) = dead;
            birdnext(j, 2) = dead;
            eaten = eaten + 1;
            %kappa = kappa - 0.25;
        end
    end
    
    % Store the information of birdnext in the positions array and
    % increment t.
    positions(:,:,i + 1) = birdnext;
    t = t + h;
  end
end