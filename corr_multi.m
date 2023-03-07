clear
close all

nSeq=149; %Longitud de la secuencia
root=primes(233); %semilla de la secuencia
bitRound=10; %número de bits del conversor
nSim=5;      %número de ciclos por símbolo
signalNoiseR=100;   %Relacion señal ruido para el ruido blanco gaussiano
nTransmisores=5;    %Numero de transmisores

fc=110e3; %frecuencia de la señal portadora
tc=1/fc; %periodo de la señal portadora
fs=1e6; %Frecuencia de muestreo
ts=0:1/fs:tc; %Vector temporal de un ciclo de portadora
t=0:1/fs:nSim*length(nSeq*nTransmisores)*(ts(length(ts))+1/fs)-1/fs; %Vector temporal de la secuencia completa
sI=kron(ones(1,nSim),square(2*pi*fc*ts)); %onda cuadrada en fase
sQ=kron(ones(1,nSim),square(2*pi*fc*ts+pi/2)); %onda cuadrada en cuadratura

ZCseq=zeros(nTransmisores,nSeq);
for i=1:nTransmisores
    ZCseq(i,:)=zadoffChuSeq(root(i),nSeq)';
end

modZCseq=zeros(nTransmisores,nSeq*length(sI));
for i=1:nTransmisores
    modZCseq(i,:)=modularSecuencia(ZCseq(i,:),sI,sQ,0,signalNoiseR);
end

modTotalSeq=reshape(modZCseq',1,[]);

corr=zeros(nTransmisores,2*length(modTotalSeq)-1);
lag=corr;

for i=1:nTransmisores
    [corr(i,:),lag(i,:)]=xcorr(modTotalSeq,modZCseq(i,:));
end

figure;
for i=1:nTransmisores
    subplot(nTransmisores,1,i);plot(lag(i,:),corr(i,:));
end

figure;
fourier_transform(modTotalSeq,fs,'frec');

function modSeq = modularSecuencia(seq,sampleI,sampleQ,bR,noise)
    xI=reshape(kron(imag(seq),sampleI)',1,[]);
    xQ=reshape(kron(real(seq),sampleQ)',1,[]);

    modSeq=xI-xQ; %Secuencia modulada con onda cuadrada
    
    if bR>0
        modSeq=round((modSeq-min(modSeq))/(max(modSeq)-min(modSeq))*(2^bR-1))-(2^(bR-1)-1/2);  %Se simula el efecto de cuantizar la secuencia
        modSeq=modSeq/(2^bR-1); %Se devuelve a la secuencia a valores aproximadamente de la misma magnitud original
    end
    modSeq=modSeq/(abs(max(modSeq)));

    modSeq=awgn(modSeq,noise,'measured');
end