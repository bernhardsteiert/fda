%% Color circle HSV
% Adapted from: http://stackoverflow.com/questions/3339692/modeling-hsv-color-space-in-matlab

radius = 0:0.2:1;
% radius = sqrt(radius);
angle = 0:0.05:1;

X = radius'*cos(2*pi*angle);
Y = radius'*sin(2*pi*angle);

H = repmat(linspace(0,1,length(angle)),length(radius),1);
S = ones(length(radius),length(angle));
V = repmat(linspace(0,1,length(radius))',1,length(angle));

hsvImage = cat(3,H,S,V);
C = hsv2rgb(hsvImage);

close all

figure
hold on

for r = 1:length(radius)-1
    for a = 1:length(angle)-1
        x = X(r:r+1,a:a+1);
        y = Y(r:r+1,a:a+1);
        fill(x([1 2 4 3]),y([1 2 4 3]),C(r,a,:))
        % patch works just as good ...
    end
end

axis equal