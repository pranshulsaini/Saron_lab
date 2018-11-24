% SNR_CORR is a program to compute the corollation map from mixing matrix
% outputed from BSS/ICA/SOBI for estimation of SNR and showing the efficacy
% of denoising procedure
% Date Modified: May 29 2012
% By: Iman M.Rezazadeh irezazadeh@ucdavis.edu 





load('C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI 003\GOTO SOBI\SOBI_Output\ASDMSI_003_forSOBI\ASDMSI_003_forSOBI_W.mat')
no_source=size(W,1);
Corr_Mat=zeros(no_source,no_source);
for i=1:no_source
    for j=1:no_source
        temp=xcorr(W(:,i),W(:,j));
        Corr_Mat(i,j)=temp(129);
    end
end
Corr_Mat=abs(Corr_Mat);
max_corr=max(max(Corr_Mat));
min_corr=min(min(Corr_Mat));

scale=255/(max_corr);
scaled_Corr_Mat=(Corr_Mat)*scale;
scaled_Corr_Mat=log10(scaled_Corr_Mat);
subplot(311);imshow(scaled_Corr_Mat); title('Pre-DENOISING After SOBI Normalized Correlation Matrices of Sources')
xlabel('Source No 1-128');ylabel('Source No 1-128');
before_corr=scaled_Corr_Mat;

load('C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI 003\GOTO SOBI\SOBI_Output\ASDMSI_003_forSOBI\003_VOTES_IMAN05022012forRecon\with 112\again SOBI\output\recon_ASDMSI_003_stim_forSOBI_ACCA\W_scaled.mat')
W=W_scaled;
no_source=size(W,1);
Corr_Mat=zeros(no_source,no_source);
for i=1:no_source
    for j=1:no_source
        temp=xcorr(W(:,i),W(:,j));
        Corr_Mat(i,j)=temp(129);
    end
end
Corr_Mat=abs(Corr_Mat);
max_corr=max(max(Corr_Mat));
min_corr=min(min(Corr_Mat));
scale=255/(max_corr);
scaled_Corr_Mat=(Corr_Mat)*scale;
scaled_Corr_Mat=log10(scaled_Corr_Mat);
subplot(312);imshow(scaled_Corr_Mat);title('Post-DENOISING After SOBI Normalized Correlation Matrices of Sources')
after_corr=scaled_Corr_Mat;
xlabel('Source No 1-128');ylabel('Source No 1-128');



diff=after_corr-before_corr;
%max_corr=max(max(diff));
%min_corr=min(min(diff));
%scale=255/(max_corr-min_corr)/255;
scaled_Corr_Mat=(abs(diff));
%scaled_Corr_Mat=log10(scaled_Corr_Mat);
subplot(313);imshow(scaled_Corr_Mat);title('Difference Between Pre and Post denoising on  Normalized Correlation Matrices of Sources')
xlabel('Source No 1-128');ylabel('Source No 1-128');