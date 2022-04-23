% MU-MISO

clc
clear
close all

M = 4;
N = 1;
K = 4;

Hnum = 500;
snr = -10:5:30;
sigma2 = 1;
sumratedataMRT = zeros(1,length(snr));
sumratedataZF = zeros(1,length(snr));
sumratedataRZF = zeros(1,length(snr));
sumratedataWMMSE = zeros(1,length(snr));
for idx1 = 1:1:length(snr)
    disp(['SNR: ' num2str(snr(idx1))])
    Pow = sigma2*10^(snr(idx1)/10);
    SRrzf = 0;
    SRzf = 0;
    SRmrt = 0;
    SRwmmse = 0;
    for idx2 = 1:1:Hnum
        H = sqrt(1/2)*(randn(K,M) + 1i*randn(K,M));
        % MRT
        PMRT = MRT(H,Pow);
        % ZF
        PZF = ZF(H,Pow);
        % MMSE
        PRZF = RZF(H,Pow);
        [Bwmmse,SR] = WMMSEpreandSR(H,H,Pow,ones(1,K),M,N,K,20);
        SRzf = SRzf + SumRate(H,PZF,sigma2);
        SRrzf = SRrzf + SumRate(H,PRZF,sigma2);
        SRmrt = SRmrt + SumRate(H,PMRT,sigma2);
        SRwmmse = SRwmmse + SR(end);
     end
    sumratedataMRT(idx1) = SRmrt/Hnum;
    sumratedataZF(idx1) = SRzf/Hnum;
    sumratedataRZF(idx1) = SRrzf/Hnum;
    sumratedataWMMSE(idx1) = SRwmmse/Hnum;
end

sumratedataMRT
sumratedataZF
sumratedataRZF
sumratedataWMMSE

figure
plot(snr,sumratedataZF,'r-*','LineWidth',2)
hold on
plot(snr,sumratedataRZF,'k-+','LineWidth',2)
plot(snr,sumratedataMRT,'g->','LineWidth',2)
plot(snr,sumratedataWMMSE,'b-^','LineWidth',2)
xlabel('SNR (dB)')
ylabel('Sum Rate (bps/Hz)')
legend('ZF','MMSE','MRT','WMMSE','Location','best')
grid on