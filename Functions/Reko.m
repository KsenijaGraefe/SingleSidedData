function [S,u,c] = Reko(SM_particle,SM_empty_1,particle2,particle3,e,p)
%% SNR -> sort auf the frequency components
[gutfreqX_alle,gutfreqY_alle] = SNR(SM_particle,SM_empty_1(:,:,:,1:(size(SM_particle,4))),3.65,3.65,120);
gutfreqX=gutfreqX_alle(~(gutfreqX_alle>635));
gutfreqY=gutfreqY_alle(~(gutfreqY_alle>635));
%% system matrix
SystemMatrix_2 = ((SM_particle(1:2:end,:,:,:) + 1i*SM_particle(2:2:end,:,:,:)));
SystemMatrixX_1 = SystemMatrix_2(:,:,1:size(SystemMatrix_2,4)/2);
SystemMatrixY_1 = SystemMatrix_2(:,:,size(SystemMatrix_2,4)/2+1:size(SystemMatrix_2,4));

for k = 1:size(SystemMatrix_2,1)
    for j = 1:size(SystemMatrix_2,2)
    SystemMatrixX_1(k,j,:) = reshape(SystemMatrix_2(k,j,1:size(SystemMatrix_2,4)/2),1,size(SystemMatrix_2,4)/2);
    SystemMatrixY_1(k,j,:) = reshape(SystemMatrix_2(k,j,size(SystemMatrix_2,4)/2+1:size(SystemMatrix_2,4)),1,size(SystemMatrix_2,4)/2);
    end
end

for a = 1:size(gutfreqX,2)
    SystemMatrixX(:,:,(a)) = SystemMatrixX_1(:,:,gutfreqX(a));
end
for a = 1:size(gutfreqY,2)
    SystemMatrixY(:,:,(a)) = SystemMatrixY_1(:,:,gutfreqY(a));
end

SystemMatrix_1 = cat(3,SystemMatrixX,SystemMatrixY);

%% phantom measurement
empty = particle3(:,e,1:2:end) + 1i*particle3(:,e,2:2:end);
phantom_2_partikel = (particle2(:,p,1:2:end) + 1i*particle2(:,p,2:2:end))-empty;
phantomX_1 = phantom_2_partikel(:,1:size(phantom_2_partikel,3)/2);
phantomY_1 = phantom_2_partikel(:,size(phantom_2_partikel,3)/2+1:size(phantom_2_partikel,3));

for a = 1:size(gutfreqX,2)
       phantomX(:,a) = phantomX_1(:,gutfreqX(a));
end
for a = 1:size(gutfreqY,2)
    phantomY(:,a) = phantomY_1(:,gutfreqY(a));
end
u = [phantomX,phantomY].';
%% reconstruction

    SM = reshape(SystemMatrix_1,size(SystemMatrix_1,1)*size(SystemMatrix_1,2),size(SystemMatrix_1,3));%SystemMatrix_1(:,m,:);%aussortiertSM(m,:,:); % Warum ist m hier an der y Position und im folgenden an der x position?????
    S = permute(SM,[2 1]);
   
   
%% weighting matrix and regularization
w = [];
    SH = ctranspose(S);

    for j= 1:size(S,1)
        w(j) = sqrt(sum((abs(S(j,:))).^2));
    end

    W = zeros(size(S,1),size(S,1));
    W = bsxfun(@times,(1./w.^2),diag(ones(1,size(S,1))));

[U,s,V] = csvd((W*S));
[x_lambda] = tikhonov(U,s,V,u,s(end));
I = diag(ones(1,size(SystemMatrix_1,1)*size(SystemMatrix_1,2)));
lambda = mean(abs(x_lambda));

%% ART
    k = 5;
    ARegularisiert = (SH*W*S+lambda*I);
    bRegularisiert = (SH*W*u);
    [XRegularisiert(:,:,:),XbisRegularisiert(:,:,:),rho,eta] = art_Ksenija(ARegularisiert,bRegularisiert,k);
%% images
figure
B=real(XbisRegularisiert(:,k));
c=reshape(B,size(SystemMatrix_1,1),size(SystemMatrix_1,2),1);
imagesc(c)
colormap('gray')
colorbar

xlabel('y-direction / mm','FontSize',20);
ylabel('x-direction / mm','FontSize',20);

