v=VideoReader('autumn_in_central_park_nyc.mp4')
load kobe32_cacti.mat

Y=zeros(270,275,100);
currAxes = axes;
i=1;
while hasFrame(v) && i <= 600
    vidFrame = readFrame(v);
    image(vidFrame, 'Parent', currAxes);
    currAxes.Visible = 'off';
    pause(1/v.FrameRate);
    
   
    Y(:,:,i) = double(rgb2gray(vidFrame(1:4:1080, 1:7:1920, :)));
    i = i+1;
   
end


%% Measurements
measure = zeros(256,256);
orig = zeros(256, 256, 8);
for i= 537:544
    orig(:,:,i-536) = Y(1:256,1:256,i);
    measure = measure + Y(1:256,1:256,i) .* mask(:,:,i-536);
end

meas = measure;
save meas meas
save orig orig;
save mask mask;





   
    