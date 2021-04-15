
%Specify # of electrons and timesteps (should be same as julia code)
num_e = 10;
num_steps=1000;

%create axis, animated lines, and video file
axis([-80,80,-80,80])
g = arrayfun(@(x) animatedline(), 1:num_e);
v = VideoWriter('C:\Users\drewj\Documents\atom\uropjulia\uropPlot','MPEG-4');
open(v);

%change width of lines
for i=1:length(g)
    g(i).LineWidth=1;
end
%plot each timestep
for k = 1:num_steps-1
    count=0;
    for e =1:num_e
        xvec = posx(count+k:count+k+1);
        yvec = posy(count+k:count+k+1);
        addpoints(g(e),xvec,yvec)
        drawnow
        count=count+num_steps;
    end
    
    % Capture the plot as an image, write to video
     frame = getframe(gcf);
     writeVideo(v,frame);
end
close(v);
