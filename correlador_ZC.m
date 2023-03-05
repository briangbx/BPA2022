clear
close all

Nseq=149; %Longitud de la secuencia
root=17; %semilla de la secuencia
ZCseq=zadoffChuSeq(root,Nseq)';
bitRound=10; %número de bits del conversor
nSim=1;      %número de ciclos por símbolo

totalSeq=ZCseq;

fc=110e3; %frecuencia de la señal portadora
tc=1/fc; %periodo de la señal portadora
fs=1e6; %Frecuencia de muestreo
ts=0:1/fs:tc; %Vector temporal de un ciclo de portadora
t=0:1/fs:nSim*length(totalSeq)*(ts(length(ts))+1/fs)-1/fs; %Vector temporal de la secuencia completa

sI1=kron(ones(1,nSim),square(2*pi*fc*ts)); %onda cuadrada en fase
sQ1=kron(ones(1,nSim),square(2*pi*fc*ts+pi/2)); %onda cuadrada en cuadratura

sI2=kron(ones(1,nSim),sin(2*pi*fc*ts)); %onda senoidal en fase
sQ2=kron(ones(1,nSim),sin(2*pi*fc*ts+pi/2)); %onda senoidal en cuadratura

xI1=reshape(kron(imag(totalSeq),sI1)',1,[]);
xQ1=reshape(kron(real(totalSeq),sQ1)',1,[]);

xI2=reshape(kron(imag(totalSeq),sI2)',1,[]);
xQ2=reshape(kron(real(totalSeq),sQ2)',1,[]);

modSeq1=xI1-xQ1; %Secuencia modulada con onda cuadrada
modSeq2=xI2-xQ2; %Secuencia modulada con onda senoidal

min1=min(modSeq1);min2=min(modSeq2);
max1=max(modSeq1);max2=max(modSeq2);

figure
subplot(1,4,1);plot(modSeq1);axis([140 170 -inf inf]) 
modSeq1=round((modSeq1-min1)/(max1-min1)*(2^bitRound-1))-(2^(bitRound-1)-1/2);  %Se simula el efecto de cuantizar la secuencia
modSeq2=round((modSeq2-min2)/(max2-min2)*(2^bitRound-1))-(2^(bitRound-1)-1/2);
hold on
subplot(1,4,2);
plot(modSeq1);axis([140 170 -inf inf]) 

modSeq1=modSeq1/(2^bitRound-1); %Se devuelve a la secuencia a valores aproximadamente de la misma magnitud original
modSeq2=modSeq2/(2^bitRound-1);

subplot(1,4,3);plot(modSeq1);axis([140 170 -inf inf])

modSeq1=awgn(modSeq1,100,'measured');
modSeq2=awgn(modSeq2,100,'measured');

subplot(1,4,4);plot(modSeq1);axis([140 170 -inf inf])

modZCseq1=reshape(kron(imag(ZCseq),sI1)',1,[])-reshape(kron(real(ZCseq),sQ1)',1,[]);
modZCseq2=reshape(kron(imag(ZCseq),sI2)',1,[])-reshape(kron(real(ZCseq),sQ2)',1,[]);

[c1,lag1]=xcorr(modSeq1,modZCseq1); %correlación de las señales
[c2,lag2]=xcorr(modSeq2,modZCseq2);

%c1=c1/(abs(max(c1)));
%c2=c2/(abs(max(c2)));

figure;
subplot(2,1,1);plot(lag1,c1);title("Correlacion de la secuencia modulada con onda cuadrada");
subplot(2,1,2);plot(lag2,c2);title("Correlacion de la secuencia modulada con onda senoidal");

%figure;
%subplot(2,1,1);plot(t,modSeq1); 
%subplot(2,1,2);plot(t,modSeq2);

figure;
fourier_transform(modSeq1,fs,'frec');title("TF de la secuencia modulada con onda cuadrada")
figure;
fourier_transform(modSeq2,fs,'frec');title("TF de la secuencia modulada con onda senoidal")