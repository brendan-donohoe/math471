function savemovie(xmin, xmax, ymin, ymax, nbirds, cx, cy, gamma1, gamma2, kappa, ro, lambda, delta, h, endt, filename, pursuit)
  fprintf('Running simulation...\n');
  positions = getbirdpositions(xmin, xmax, ymin, ymax, nbirds, cx, cy, gamma1, gamma2, kappa, ro, lambda, delta, h, endt, pursuit);
  
  fprintf('Creating movie...\n');
  numframes = length(positions(1,1,:));
  mov(numframes) = struct('cdata',[],'colormap',[]);
  
  t = 0;
  for i=1:numframes
    curpositions = positions(:,:,i);
    mov(i) = getbirdframe(xmin, xmax, ymin, ymax, curpositions(:, 1)', curpositions(:, 2)', cx(t), cy(t));
    t = t + h;
  end
  
  fprintf('Saving movie...\n');
  name = strcat(filename, '.avi');
  video = VideoWriter(name, 'Uncompressed AVI');
  video.FrameRate = 12;
  open(video);
  writeVideo(video, mov);
  close(video);
  fprintf('Complete!\n');
end