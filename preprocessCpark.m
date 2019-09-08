v=VideoReader('autumn_in_central_park_nyc.mp4')
load kobe32_cacti.mat

X=zeros(270,275,100);
currAxes = axes;
i=1;
while hasFrame(v) && i <= 100
    vidFrame = readFrame(v);
    image(vidFrame, 'Parent', currAxes);
    currAxes.Visible = 'off';
    pause(1/v.FrameRate);
    
   
    X(:,:,i) = double(rgb2gray(vidFrame(1:4:1080, 1:7:1920, :)));
    i = i+1;
   
end


%% Measurements
measure = zeros(256,256);
orig = zeros(256, 256, 8);
for i= 90:97
    orig(:,:,i-89) = X(1:256,1:256,i);
    measure = measure + X(1:256,1:256,i) .* mask(:,:,i-89);
end

meas = measure;
save meas meas

save orig orig;
save mask mask;



   
    