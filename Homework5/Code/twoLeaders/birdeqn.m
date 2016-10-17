function slope = birdeqn(t, b, gamma, othergamma, leader, otherleader, kappa, center, ro, delta, neighbors)
  leaderterm = gamma * (leader - b);
  flockterm = kappa * (center - b);
  otherleaderterm = othergamma * (otherleader - b); 
  repulsionterm = 0;
  for bm = neighbors'
    s = (b - bm) / ((b - bm)^2 + delta);
    repulsionterm = repulsionterm + s;
  end
  repulsionterm = repulsionterm * ro;
  slope = leaderterm + flockterm + repulsionterm + otherleaderterm;
end