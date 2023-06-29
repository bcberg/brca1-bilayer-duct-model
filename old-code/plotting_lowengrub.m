
figure(1)
plot(tA1/1.2,solA1(:,2),'k')
hold on
plot((tA2+100)/1.2,solA2(:,2),'k')
plot((tB+200)/1.2,solB(:,2),'k')

plot(tA1/1.2,solA1(:,1),'--k')
hold on
plot((tA2+100)/1.2,solA2(:,1),'--k')
plot((tB+200)/1.2,solB(:,1),'--k')

figure(2)
plot(tA1/1.2,solA1(:,3),'k')
hold on
plot((tA2+100)/1.2,solA2(:,3),'k')
plot((tB+200)/1.2,solB(:,3),'k')