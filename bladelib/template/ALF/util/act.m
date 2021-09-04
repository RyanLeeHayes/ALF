ex=load('132/ex.dat');
nreps=load('../nreps');

for irep=0:(nreps-1)
  s=(ex==irep)*((0:(nreps-1))');
  ss=s(1:(end-50));
  for i=1:50
    sr=s((1:(end-50))+i);
    ac(irep+1,i)=mean(ss==sr);
  end
end

figure(1)
colors=hsv(nreps);
hold off
for i=1:nreps
  plot(0.001:0.001:0.05,ac(i,:),'Color',colors(i,:))
  hold on
end
xlabel('time ns')
ylabel('Autocorrelation')
