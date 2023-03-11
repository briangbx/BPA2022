clear
close all

nSeq=349; %Longitud de la secuencia
root=primes(nSeq); %semilla de la secuencia
bitRound=10; %número de bits del conversor
nSim=1;      %número de ciclos por símbolo
signalNoiseR=100;   %Relacion señal ruido para el ruido blanco gaussiano
nTransmisores=3;    %Numero de transmisores
gap=200;            %Desfase entre señales recibidas

fc=110e3; %frecuencia de la señal portadora
tc=1/fc; %periodo de la señal portadora
fs=1e6; %Frecuencia de muestreo
ts=0:1/fs:tc; %Vector temporal de un ciclo de portadora
t=0:1/fs:nSim*length(nSeq*nTransmisores)*(ts(length(ts))+1/fs)-1/fs; %Vector temporal de la secuencia completa
sI=kron(ones(1,nSim),square(2*pi*fc*ts)); %onda cuadrada en fase
sQ=kron(ones(1,nSim),square(2*pi*fc*ts+pi/2)); %onda cuadrada en cuadratura

gap=gap*length(sI); %Corrección del desfasaje

% Generación de las N secuencias
ZCseq=zeros(nTransmisores,nSeq);
for i=1:nTransmisores
    ZCseq(i,:)=zadoffChuSeq(root(i),nSeq)';
end

% Modulación de las N secuencias patron
modZCseq=zeros(nTransmisores,nSeq*length(sI));
for i=1:nTransmisores
    modZCseq(i,:)=modularSecuencia(ZCseq(i,:),sI,sQ,bitRound,0);
end

%modTotalSeq=reshape(modZCseq',1,[]);

%Generación de una única secuencia a correlar. Superponiendo las secuencias
%patrón

modTotalSeq=zeros(1,length(modZCseq)+gap*(nTransmisores-1));

for i=1:nTransmisores
    desfase1=zeros(1,(i-1)*gap);
    desfase2=zeros(1,(nTransmisores-i)*gap);
    modTotalSeq=modTotalSeq+[desfase1 modZCseq(i,:) desfase2];
end

modTotalSeq=awgn(modTotalSeq,signalNoiseR,'measured');

corr=zeros(nTransmisores,2*length(modTotalSeq)-1);
lag=corr;

for i=1:nTransmisores
    [corr(i,:),lag(i,:)]=xcorr(modTotalSeq,modZCseq(i,:));
end

figure;
for i=1:nTransmisores
    subplot(nTransmisores,1,i);plot(lag(i,:),corr(i,:));
end

% figure;
% fourier_transform(modTotalSeq,fs,'frec');

function modSeq = modularSecuencia(seq,sampleI,sampleQ,bR,noise)
    xI=reshape(kron(imag(seq),sampleI)',1,[]);
    xQ=reshape(kron(real(seq),sampleQ)',1,[]);

    modSeq=xI-xQ; %Secuencia modulada con onda cuadrada
    
    if bR>0
        modSeq=round((modSeq-min(modSeq))/(max(modSeq)-min(modSeq))*(2^bR-1))-(2^(bR-1)-1/2);  %Se simula el efecto de cuantizar la secuencia
        modSeq=modSeq/(2^bR-1); %Se devuelve a la secuencia a valores aproximadamente de la misma magnitud original
    end
    modSeq=modSeq/(abs(max(modSeq)));
    
    if noise~=0
    modSeq=awgn(modSeq,noise,'measured');
    end
end