clear
%close all

bitADC=10; %número de bits del conversor
bitFracc_entrada=bitADC-1; %número de bits que representan la parte fraccionaria---Se define de esta forma para normalizar la entrada y salida
bitFracc_salida=4;  %número de bits que representan la parte fraccionaria luego de la correlación
bitMax_salida=18; %número máximo de bits que se permite durante el proceso de correlación

nSeq=349; %Longitud de la secuencia
root=primes(nSeq); %semilla de la secuencia
nMultiplicadores=300; %Indica el número de multiplicadores de la FPGA que se pueden usar en simultáneo.
nSim=1;      %número de ciclos por símbolo
signalNoiseR=100;   %Relacion señal ruido para el ruido blanco gaussiano
nTransmisores=3;    %Numero de transmisores
gap=100;            %Desfase entre señales recibidas

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
    modZCseq(i,:)=modularSecuencia(ZCseq(i,:),sI,sQ,bitADC);
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
seqPadded=[zeros(1,length(modTotalSeq)-1) modTotalSeq zeros(1,length(modZCseq)-1)];
cola=zeros(1,length(modTotalSeq));
bitsTruncados=2*bitFracc_entrada-bitFracc_salida; %Se definen los bits que se deben despreciar luego de la multiplicación
if bitsTruncados<0
    bitsTruncados=0;
end

for i=1:nTransmisores
    secPat=modZCseq(i,:); %Se toma una de las secuencias patrón para correlar
    for j=1:2*length(modTotalSeq)-1
        cola=seqPadded(j:length(modZCseq)-1+j); %Se arma un vector que se correlará con secPat de igual longitud 
        
        %Se realiza la multiplicación elemento a elemento por etapas
        for k=1:nEtapas-1
            prodParcial=sum(cola(1+(k-1)*nMultiplicadores:k*nMultiplicadores).*secPat(1+(k-1)*nMultiplicadores:k*nMultiplicadores),'all');
            prodParcial=round(prodParcial/2^bitsTruncados); %Se desplaza hacia la izquierda la secuencia de bits, truncando su valor.
            corr(i,j)=corr(i,j)+prodParcial;
            
            %Si el resultado parcial de la correlación es mayor que el
            %mayor valor representable con los bits Max, se trunca.
            if corr(i,j)>=2^(bitMax_salida-1)
                corr(i,j)=2^(bitMax_salida-1)-1;
            elseif corr(i,j)<-2^(bitMax_salida-1)
                corr(i,j)=-2^(bitMax_salida-1);
            end
        end
        prodParcial=sum(cola(1+(nEtapas-1)*nMultiplicadores:end).*secPat(1+(nEtapas-1)*nMultiplicadores:end),'all');
        prodParcial=round(prodParcial/2^bitsTruncados);
        corr(i,j)=corr(i,j)+prodParcial;
        if corr(i,j)>=2^(bitMax_salida-1)
            corr(i,j)=2^(bitMax_salida-1)-1;
        elseif corr(i,j)<-2^(bitMax_salida-1)
            corr(i,j)=-2^(bitMax_salida-1);
        end
    end
    
    %Se devuelven los valores binarios a decimal.
    for j=1:length(corr(i,:))
        neg=1;
        if corr(i,j)<0
            neg=-1;
        end
        binN=dec2bin(abs(corr(i,j)),bitFracc_salida);
        decN=bin2dec(binN(length(binN)-bitFracc_salida+1:end));
        intN=0;
        if length(binN)>bitFracc_salida
            intN=neg*bin2dec(binN(1:length(binN)-bitFracc_salida));
        end
        corr(i,j)=intN+decN/(2^(bitFracc_salida));
    end
end
    
    
figure;
for i=1:nTransmisores
    subplot(nTransmisores,1,i);plot(corr(i,:));
end

% figure;
% fourier_transform(modTotalSeq,fs,'frec');

function modSeq = modularSecuencia(seq,sampleI,sampleQ,bR)
    xI=reshape(kron(imag(seq),sampleI)',1,[]);
    xQ=reshape(kron(real(seq),sampleQ)',1,[]);

    modSeq=xI-xQ; %Secuencia modulada con onda cuadrada
    
    ppSeq=max(modSeq)-min(modSeq);
    
    modSeq=round((modSeq-min(modSeq))/ppSeq*(2^bR))-(2^(bR-1));  %Se simula el efecto de cuantizar la secuencia
    modSeq(modSeq==2^(bR-1))=2^(bR-1)-1;
end