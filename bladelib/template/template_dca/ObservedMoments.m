m1=load('m1_a.obs.dat');
m2=load('m2_a.obs.dat');

m1p=load('m1_a.LM.dat');
m2p=load('m2_a.LM.dat');

fig=1
figure(fig)
loglog(m1,m1p,'.')
set(fig,'PaperSize',[3.25,2.5])
set(fig,'PaperPosition',[0,0,3.25,2.5])
title('First Moment')
xlabel('Observed p(A)')
ylabel('DCA p(A)')
saveas(fig,'m1.pdf')

fig=2
figure(fig)
loglog(m2,m2p,'.')
set(fig,'PaperSize',[3.25,2.5])
set(fig,'PaperPosition',[0,0,3.25,2.5])
title('Second Moment')
xlabel('Observed p(A,B)')
ylabel('DCA p(A,B)')
saveas(fig,'m2.pdf')

