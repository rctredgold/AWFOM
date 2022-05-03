a = linspace(0,1,10001);

psi = 0.7545;

f = 0.9;
ct1d = 4.*a.*(1-a);
cp = 4.*a.*((1-a).^2);
cpnorm = (4.*a.*((1-a).^2)).*psi;
ctBuhl(1:4000) = 4.*a(1:4000).*f.*(1-a(1:4000));
ctBuhl(4001:10001) =  0.8889 +  ((4.*f -4.4444).*a(4001:10001)) + ((5.5556 - 4*f).*(a(4001:10001).^2));

figure(1)
plot(a,ct1d);
hold on
plot(a,ctBuhl);
plot(a,cp);
plot(a,cpnorm)
xlim([0,1]);
hold off



