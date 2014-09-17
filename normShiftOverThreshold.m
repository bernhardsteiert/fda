x = linspace(-4,8,5001);
dx = x(2)-x(1);

thres = linspace(0,3,201);

y = normpdf(x);
normfac = 1-normcdf(thres);

expect = nan(size(thres));
for i = 1:length(thres)
    
    myind = x > thres(i);
    expect(i) = x(myind)*y(myind)'*dx / normfac(i);
    
end

figure
plot(thres,expect)
set(gca,'XDir','reverse')


figure

h = nan(1,6);
h(1) = subplot(2,3,1);
mypos = get(gca,'Position');

plot(x,y,'k')
hold on
plot([1 1],[0 .5],'k--')
set(gca,'XLim',[-4 6])
annotation('arrow',mypos(1)+mypos(3)*[.4 .8],mypos(2)+[.5 .5]*mypos(4))
ylabel('PDF')
xlabel('pulsatile score')

h(2) = subplot(2,3,2);
mypos = get(gca,'Position');

y2 = normpdf(x-4);
plot(x,.7*y+.3*y2,'k')
hold on
plot([2 2],get(gca,'YLim'),'k--')
set(gca,'XLim',[-4 8])
annotation('arrow',mypos(1)+mypos(3)*[1/3 1/3],mypos(2)+[.6 .4]*mypos(4))
annotation('arrow',mypos(1)+mypos(3)*[2/3 2/3],mypos(2)+[0.05 .25]*mypos(4))

h(3) = subplot(2,3,3);
mypos = get(gca,'Position');

plot(x,.7*y+.3*y2,'k')
hold on
plot([2 2],get(gca,'YLim'),'k--')
set(gca,'XLim',[-4 8])
annotation('arrow',mypos(1)+mypos(3)*[1/3 1/3],mypos(2)+[.6 .4]*mypos(4))
annotation('arrow',mypos(1)+mypos(3)*[2/3 2/3],mypos(2)+[.05 .25]*mypos(4))
annotation('arrow',mypos(1)+mypos(3)*[2/3 .9],mypos(2)+[.15 .15]*mypos(4))

h(4) = subplot(2,3,4);
zo = [0 1];
plot(zo,.1+.8*zo,'b')
hold on
plot(zo,.95-.9*zo,'r')
ylabel('Mean')
xlabel('Ligand dose')

h(5) = subplot(2,3,5);
plot(zo,.1+.8*zo,'b')
hold on
plot(zo,.5+0*zo,'r')

h(6) = subplot(2,3,6);
plot(zo,.1+.8*zo,'b')
hold on
plot(zo,.05+.65*zo,'r')

legend('fraction','features')

for i = 1:length(h)
    set(h(i),'XTick',[],'YTick',[])
end