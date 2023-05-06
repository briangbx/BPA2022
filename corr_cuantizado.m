clear
%close all

nSeq=349; %Longitud de la secuencia
root=primes(nSeq); %semilla de la secuencia
bitRound=10; %número de bits del conversor
bitFracc=bitRound-1; %número de bits que representan la parte fraccionaria---Se define de esta forma para normalizar la entrada y salida
bitMax=18; %longitud máxima del vector de entrada de los multiplicadores
nMultiplicadores=300; %Indica el número de multiplicadores de la FPGA que se pueden usar en simultáneo.
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
    modZCseq(i,:)=modularSecuencia(ZCseq(i,:),sI,sQ,bitRound,1);
end

%Generación de una única secuencia a correlar. Superponiendo las secuencias
%patrón

modTotalSeq=zeros(1,length(modZCseq)+gap*(nTransmisores-1));

for i=1:nTransmisores
    desfase1=zeros(1,(i-1)*gap);
    desfase2=zeros(1,(nTransmisores-i)*gap);
    modTotalSeq=modTotalSeq+[desfase1 modZCseq(i,:) desfase2];
    clear desfase1 desfase2
end

%modTotalSeq=awgn(modTotalSeq,signalNoiseR,'measured');

%Correlación entre la secuencia total con cada secuencia patrón

corr=zeros(nTransmisores,2*length(modTotalSeq)-1);
%Se calcula la cantidad de etapas en las que se deben dividir las
%multiplicaciones
nEtapas=fix(length(modZCseq)/nMultiplicadores);
if nEtapas~=length(modZCseq)/nMultiplicadores
    nEtapas=nEtapas+1;
end

%Se crea un vector con la secuencia completa rodeada por ceros de donde se
%tomaran los vectores cola
seqPad=[zeros(1,length(modTotalSeq)-1) modTotalSeq zeros(1,length(modZCseq)-1)];
cola=zeros(1,length(modTotalSeq));

for i=1:nTransmisores
    SecPat=modZCseq(i,:);
    for j=1:2*length(modTotalSeq)-1
        cola=seqPad(j:length(modZCseq)-1+j);
        for k=1:nEtapas-1
            corr(i,j)=corr(i,j)+sum(cola(1+(k-1)*nMultiplicadores:k*nMultiplicadores).*SecPat(1+(k-1)*nMultiplicadores:k*nMultiplicadores),'all');
        end
        corr(i,j)=corr(i,j)+sum(cola(1+(nEtapas-1)*nMultiplicadores:end).*SecPat(1+(nEtapas-1)*nMultiplicadores:end),'all');
    end
    
    %Se devuelven los valores binarios a decimal.
    for m=1:length(corr(i,:))
        neg=1;
        if corr(i,m)<0
            neg=-1;
        end
        binN=dec2bin(abs(corr(i,m)),2*bitFracc);
        decN=bin2dec(binN(length(binN)-2*bitFracc+1:end));
        intN=0;
        if length(binN)>2*bitFracc
            intN=neg*bin2dec(binN(1:length(binN)-2*bitFracc));
        end
        corr(i,m)=intN+decN/(2^(2*bitFracc));
    end
end
    
    
figure;
for i=1:nTransmisores
    subplot(nTransmisores,1,i);plot(corr(i,:));
end

% figure;
% fourier_transform(modTotalSeq,fs,'frec');

function modSeq = modularSecuencia(seq,sampleI,sampleQ,bR,int_output)
    xI=reshape(kron(imag(seq),sampleI)',1,[]);
    xQ=reshape(kron(real(seq),sampleQ)',1,[]);

    modSeq=xI-xQ; %Secuencia modulada con onda cuadrada
    
    ppSeq=max(modSeq)-min(modSeq);
    
    if bR>0
        modSeq=round((modSeq-min(modSeq))/ppSeq*(2^bR))-(2^(bR-1));  %Se simula el efecto de cuantizar la secuencia
        modSeq(modSeq==2^(bR-1))=2^(bR-1)-1;
        if int_output==0
            modSeq=modSeq/(2^bR-1)*ppSeq; %Se devuelve a la secuencia a valores aproximadamente de la misma magnitud original
        elseif int_output~=1
            error("int_output debe ser 0 o 1")
        end
    end
end