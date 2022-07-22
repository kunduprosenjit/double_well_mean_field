r1=1;
r2=2;
r3=5;

x=[0:.0001:1,1:.01:6];

y1=-(x-r1).*(x-r2).*(x-r3);
%plot(x,y1,x1(TFmin),y1(TFmin),'r*')
%plot(x,-y1)
TFmin = islocalmin(y1);
TFmax = islocalmax(y1);
plot(x,-y1,x(TFmin),-y1(TFmin),'r*',x(TFmax),-y1(TFmax),'r*')
hold on

plot([0,6],[-y1(TFmin)+1,-y1(TFmin)+1])
plot([0,6],[-y1(TFmin)-2,-y1(TFmin)-2])
plot([0,6],[-y1(TFmax)-1,-y1(TFmax)-1])