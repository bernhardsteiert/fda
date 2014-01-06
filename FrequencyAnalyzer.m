clear all;
%sites = [4	5	6	7	8	9];
sites = [37	36	35	34	33	32];


for i=1:length(sites)
    [radial_dist c_signal_woNharm] = radial_dist(sites(i));
    myrad{i} = radial_dist;
    mysignal{i} = c_signal_woNharm;
    clear radial_dist c_signal_woNharm;
end
%%
for i=1:length(sites)
    figure(2);subplot(1,6,i),plot(mysignal{i});
    title(num2str(sites(i)));
    
    figure(3);myp(i) = subplot(1,6,i),hist(myrad{i},linspace(0,0.1,32));
    xlim([0 0.1]);
    title(num2str(sites(i)));
end

old_pxx = zeros(32,1);
fs = (1/(5*60));
for i=1:length(sites)
    figure(4);s(i) = subplot(1,6,i);
end
linkaxes(s);

for i=1:length(sites)
    figure(4);s(i) = subplot(1,6,i);
    for c = 1:size(mysignal{i},2)
        x = mysignal{i}(:,c);
        [pxx,f] = periodogram(x,hamming(length(x)),length(x),fs,'onesided','power');
        old_pxx = (old_pxx*(c-1)+pxx)/c;
            
        plot(f,10*log10(old_pxx));
        title(num2str(sites(i)));
        drawnow;
    end
end
figure(5);
set(gcf,'Position',[10,200,1000,200]);
mycolor = hsv(20);
for i=1:length(sites)
    s(i) = subplot(1,6,i);
    p = randperm(size(mysignal{i},2));
    for c =  p(1:30)
        x = mysignal{i}(:,c);
        y = [0;diff(x)];
        z = [0;diff(y)];
        plot(y,z,'color',mycolor(randi(20,1),:)); hold on; axis image;
        pause(1);
        title(num2str(sites(i)));
        drawnow;
    end
    hold off;
end



