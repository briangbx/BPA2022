Nseq=79; %Longitud de la secuencia
root=31; %semilla de la secuencia
ZCseq=zadoffChuSeq(root,Nseq)';
bitRound=12; %número de bits del conversor

A=randomMatrix(1,21);
B=randomMatrix(1,30);

totalSeq=[A ZCseq B]; %Con ruido

fc=110e3; %frecuencia de la señal portadora
tc=1/fc; %periodo de la señal portadora
fs=1e6; %Frecuencia de muestreo
ts=0:1/fs:tc; %Vector temporal de un ciclo de portadora
t=0:1/fs:length(totalSeq)*(ts(length(ts))+1/fs)-1/fs; %Vector temporal de la secuencia completa

sI1=square(2*pi*fc*ts); %onda cuadrada en fase
sQ1=square(2*pi*fc*ts+pi/2); %onda cuadrada en cuadratura

sI2=sin(2*pi*fc*ts); %onda senoidal en fase
sQ2=sin(2*pi*fc*ts+pi/2); %onda senoidal en cuadratura

xI1=reshape(kron(imag(totalSeq),sI1)',1,[]);
xQ1=reshape(kron(real(totalSeq),sQ1)',1,[]);

xI2=reshape(kron(imag(totalSeq),sI2)',1,[]);
xQ2=reshape(kron(real(totalSeq),sQ2)',1,[]);

modSeq1=xI1-xQ1; %Secuencia modulada con onda cuadrada
modSeq2=xI2-xQ2; %Secuencia modulada con onda senoidal

modSeq1=round(((modSeq1+2)/4)*2^bitRound);
modSeq2=round(((modSeq2+1)/2)*2^bitRound);

modSeq1=modSeq1*4/(2^bitRound)-2;
modSeq2=modSeq2*2/(2^bitRound)-1;

modZCseq1=reshape(kron(imag(ZCseq),sI1)',1,[])-reshape(kron(real(ZCseq),sQ1)',1,[]);
modZCseq2=reshape(kron(imag(ZCseq),sI2)',1,[])-reshape(kron(real(ZCseq),sQ2)',1,[]);

modZCseq1=round(((modZCseq1+2)/4)*2^bitRound);
modZCseq2=round(((modZCseq2+1)/2)*2^bitRound);

modZCseq1=modZCseq1*4/(2^bitRound)-2;
modZCseq2=modZCseq2*2/(2^bitRound)-1;

[c1,lag1]=xcorr(modSeq1,modZCseq1);
[c2,lag2]=xcorr(modSeq1,modZCseq2);

figure;
subplot(2,1,1);plot(lag1,c1); 
subplot(2,1,2);plot(lag2,c2);

figure;
subplot(2,1,1);plot(t,modSeq1); 
subplot(2,1,2);plot(t,modSeq2);

%figure;
%fourier_transform(modSeq1,fs,'frec');title("TF de la secuencia modulada con onda cuadrada")
%figure;
%fourier_transform(modSeq2,fs,'frec');title("TF de la secuencia modulada con onda senoidal")

function R=randomMatrix(m,n)
    R=zeros(m,n);
    for i=1:m
        for j=1:n
            R(i,j)=exp(1i*2*pi*rand);
        end
    end
end