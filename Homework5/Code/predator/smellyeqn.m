function slope = smellyeqn(t, b, gamma, leader, kappa, center)
  leaderterm = gamma * (leader - b);
  flockterm = kappa * (center - b);
  slope = leaderterm + flockterm;
end