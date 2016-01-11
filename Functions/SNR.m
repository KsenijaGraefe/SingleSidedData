function [spaltex,spaltey,w,standartabweichung,SNR] = SNR_Energie_Standardabweichung(SM_particle,SM_empty_1,minSNRx,minSNRy,p)

SystemMatrix_P = ((SM_particle(1:2:end,:,:,:) + 1i*SM_particle(2:2:end,:,:,:)));
SystemMatrix_E1 = ((SM_empty_1(1:2:end,:,:,:) + 1i*SM_empty_1(2:2:end,:,:,:)));

EminusE = SystemMatrix_E1;
EminusE = squeeze(EminusE);
EminusE = reshape(EminusE,[(size(EminusE,1)*size(EminusE,2)) size(EminusE,3)]);

PminusE = SystemMatrix_P-SystemMatrix_E1;
PminusE = squeeze(PminusE);
PminusE = reshape(PminusE,[(size(PminusE,1)*size(PminusE,2)) size(PminusE,3)]);

standartabweichung = sqrt(var(EminusE));

for j= 1:size(PminusE,1)
  w(j,:) = sqrt(((abs(PminusE(j,:))).^2));
  SNR(j,:) = squeeze(w(j,:))./standartabweichung;
end
 
figure
A=minSNRx.*ones(1000,1);
semilogy(SNR([p],1:1000).')
hold all;
semilogy(A,'--','LineWidth',1.5);
xlabel('frequency components','FontSize',21);
ylabel('SNR of the x-receive channel','FontSize',21);
set(gca,'FontSize',21)
str=sprintf('snrpunkt %d', p);
axis([1 1000 0 inf])

figure
B=minSNRy.*ones(1000,1);
semilogy(SNR([p],size(SM_empty_1,4)/2+1:size(SM_empty_1,4)/2+1+999).')
hold all
semilogy(B,'--','LineWidth',1.5)
xlabel('frequency components','FontSize',21);
ylabel('SNR of the y-receive channel','FontSize',21);
set(gca,'FontSize',21) 
axis([1 1000 0 inf])

[reihex,spaltex]=find(SNR([p],1:end/2)>minSNRx);

[reihey,spaltey]=find(SNR([p],end/2+1:end)>minSNRy);