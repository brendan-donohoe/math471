function slope = predeqn(t, b, pursuit, closest, eaten)
  slope = (pursuit - (0.7 * eaten)) * (closest - b);
end