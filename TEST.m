clear all;close all;

fs = (1/(1*60));
t = 0:1/fs:3000/fs-1;
x = cos(2*pi*(1/(20*60))*t)+sin(2*pi*(1/(50*60))*t);

[pxx,f] = periodogram(x,hamming(length(x)),length(x),fs);
plot(f,10*log10(pxx)); hold on;
xlabel('Hz'); ylabel('dB');
title('Periodogram with 95%-Confidence Bounds');

figure(5);

y = [0,diff(x)];
z = [0,diff(y)];
for i=1:length(x)
    plot(x(1),y(1),'ro');hold on
    plot(x(1:i),y(1:i),'k'); 
    plot(x(i),y(i),'sb');hold off;
    pause(0.05);drawnow;
end

