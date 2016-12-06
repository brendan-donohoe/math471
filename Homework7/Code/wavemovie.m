function wavemovie(numprocs, t_max) 
  fprintf('Creating movie...\n');
  %mov(t_max) = struct('cdata',[],'colormap',[]);
  
  for i=0:t_max
    plotwaves((i), numprocs);
  end
  
  fprintf('Saving movie...\n');
  %name = strcat(filename, '.avi');
  %video = VideoWriter(name, 'Uncompressed AVI');
  %video.FrameRate = 12;
  %open(video);
  %writeVideo(video, mov);
  %close(video);
  %fprintf('Complete!\n');
end