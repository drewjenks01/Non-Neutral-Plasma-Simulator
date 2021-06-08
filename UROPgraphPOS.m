%Specify # of electrons and timesteps (should be same as julia code)
 num_e = 5;
 num_steps=1000.0;


v = VideoWriter("C:\Users\drewj\Documents\atom\uropjulia\posPlots1",'MPEG-4');
open(v);


Array = readtable('positions.csv');

for k = progress(1:num_steps)
        
    x = Array{((k-1)*num_e)+1:((k-1)*num_e)+num_e,1};
    y = Array{((k-1)*num_e)+1:((k-1)*num_e)+num_e,2};
        
    scatter(x,y,3,"filled")
    xlim([-50 50])
    ylim([-50 50])
    
    
  %  drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
   
  
    
end
    
 close(v);
 
 





    