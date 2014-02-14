% Figure box: Explanation of fPCA
addpath('./Functions/')

figure

nrows = 3;
ncols = 2;

t = [linspace(-2,0,100) linspace(0,2,100)];
step1 = [ones(1,100) zeros(1,100)]*sqrt(pi)/2;
step2 = [zeros(1,100) ones(1,100)]*sqrt(pi)/2;
strans = [zeros(1,100) sin(t(101:end)*pi)];

weights = [[-1 1.2 .1];[-1.2 .7 -.05];[-.4 .3 -.15]];

hold on
colmap = lines(size(weights,1));
for iw = 1:size(weights,1)
    plot(t,weights(iw,1)*step1+weights(iw,2)*step2+weights(iw,3)*strans,'Color',colmap(iw,:))
end

figure

subplot(nrows,ncols,1)
plot(t,step1)
hold on
plot([-2 2],[0 0],':')
set(gca,'YLim',[-1.1 1.1])

subplot(nrows,ncols,3)
plot(t,step2)
hold on
plot([-2 2],[0 0],':')
set(gca,'YLim',[-1.1 1.1])

subplot(nrows,ncols,5)
plot(t,strans)
hold on
plot([-2 2],[0 0],':')
set(gca,'YLim',[-1.1 1.1])

subplot(nrows,ncols,2)
bar(1:3,weights(:,1))

subplot(nrows,ncols,4)
bar(1:3,weights(:,2))

subplot(nrows,ncols,6)
bar(1:3,weights(:,3))