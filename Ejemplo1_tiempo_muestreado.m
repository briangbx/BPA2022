Nseq=47; %Longitud de la secuencia
root=5; %semilla de la secuencia
ZCseq=zadoffChuSeq(root,Nseq)';
figure(4);
plot(real(ZCseq),imag(ZCseq),'x');
grid on;
nRound=6;

%A=randomMatrix(1,21);
%B=randomMatrix(1,30);

totalSeq=ZCseq; %Con ruido
%%
fc=110e3; %frecuencia de la señal portadora
tc=1/fc; %periodo de la señal portadora
N=10; %Muestras por ciclo
nTemp= 100; %Numero de puntos por ciclo
t1=0:tc/nTemp:tc-tc/nTemp; %paso temporal en un ciclo de la portadora
t=0:tc/nTemp:length(totalSeq)*tc-tc/nTemp; %vector temporal de secuencia completa
sI=sin(2*pi*fc*t1); %onda cuadrada en fase
sQ=sin(2*pi*fc*t1+pi/2); %onda cuadrada en cuadratura

% figure;
% plot(t1,sQ);
%%
xI=reshape(kron(imag(totalSeq),sI)',1,[]);
xQ=reshape(kron(real(totalSeq),sQ)',1,[]);

modZCseq=reshape(kron(imag(ZCseq),sI)',1,[])-reshape(kron(real(ZCseq),sQ)',1,[]);
modSeq=xI-xQ;

[c1,lag]=xcorr(modSeq,modZCseq);

figure(1);
plot(t,modSeq);

figure(2);
plot(lag,c1/max(abs(c1)));

figure;
fourier_transform(modSeq,fc*nTemp,'frec')

ts=0:tc/N:length(totalSeq)*tc-tc/nTemp;
sampledZCseq=zeros(1,length(modZCseq)*N/nTemp);
for i=0:length(sampledZCseq)-1
    sampledZCseq(1+i)=modZCseq(1+i*nTemp/N);
end
figure;
subplot(1,2,1);plot(sampledZCseq)
subplot(1,2,2);plot(modZCseq)
cola=zeros(1,length(sampledZCseq));
c=zeros(1,length(t)*N/nTemp);

for i=0:length(t)*N/nTemp-1
    cola=[round(modSeq(1+i*nTemp/N),nRound) cola(1,1:length(sampledZCseq)-1)];
    c(1+i)=xcorr(cola,flip(sampledZCseq),0);
end

figure(3);
plot(ts,c);

function R=randomMatrix(m,n)
    R=zeros(m,n);
    for i=1:m
        for j=1:n
            R(i,j)=exp(1i*2*pi*rand);
        end
    end
end
        