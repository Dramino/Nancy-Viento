clear,clc
clear all
format short
NancyT=cputime;

%disp('                                     Programa para simulaci?n de viento')
%disp('                                   Universidad Nacional Aut?noma de M?xico')
%disp('                                     Nancy Alejandra Mosqueda Santacruz')
%
%fprintf('\n\n')
%disp('Programa para calcular viento sint?tico mediante m?todo de Davenport')
%
%
%

%Vo=input('Introduzca la velocidad caracter?stica del viento en m/s:   ');
Vo=40;
%Cambio a velocidad promediada a 10 min
Uo=0.69*Vo;
w0=0.0017;
dw=0.0017;
noFrecuencias=1470;
wf=w0*noFrecuencias;
x=zeros(1,noFrecuencias);
w=linspace(w0 ,wf ,noFrecuencias);

%Función de densidad de turbulencia
for i=1:noFrecuencias
    x(i)=1220*dw*i/Uo;
end
Sr=4*((x.^2)./((1+x.^2).^(4/3)))./w;
%plot(w,Sr)
%grid on

%r=input('Introduzca la frecuencia natural de la estructura en ciclos/s:   ');
r=0.85;
% Esta parte define los coeficientes de arrastre
%coef=input('Introduzca el n?mero de coeficientes de arrastre de la estructura:   ');
numCD=37;  
%  
%  
% arrastre=zeros(coef,1);
%     for i=1:coef
%         CD(i,:)=input('Ingrese los coeficientes de arrastre de la estructura de abajo hacia arriba:   ')
%     end

%Esta parte es ?nicamente de apoyo
CD=[2.5;3.15;3.15;3.05;3.05;3.05;2.9;2.9;2.38;2.78;2.78;2.78;2.78;2.86;2.86;2.9;2.9;3.02;3.02;3.1;3.1;3.15;3.15;3.11;3.11;3.05;3.05;3.1;3.1;3.15;3.15;3.15;3.15;3.2;3.2;3.2;3.2];
    
% Esta parte define las ?reas
% numAreas=numCD;
% areas=zeros(numCD,1);
%     for i=1:numCD
%         areas(i,:)=input('Ingrese las ?reas frontales de secci?n media en m2 de abajo hacia arriba:   ')       
%     end
%Esta parte es ?nicamente de apoyo
areas=[.531;.531;.531;.602;.602;.602;.708;.708;1.239;1.4145;1.4145;1.78;1.78;1.967;1.967;2.2;2.2;2.23;2.23;2.28;2.28;2.3775;2.3775;2.765;2.765;3.247;3.247;3.341;3.341;3.353;3.353;3.59;3.59;3.59;3.58;3.808;3.808];       
%h=input('Introduzca la altura de la torre en metros:   ');
h=100.3;
%div=input('?En cu?ntas partes desea subdividir la torre?:   ');
div=17;

%por lo tanto lo que sigue es dividir h/hi para conocer la altura de cada tramo
hTramo=h/div;

hTorre=zeros(numCD,1);
hTorre(1,:)=h;
     
for i=1:10
    hTorre(i+1,:)=  hTorre(i,:)-(hTramo/3);
end
     %Esta parte es para ajustar el n?mero de coeficientes de arrastre 
     %con el m?mero de secciones en las que se desea dividir la toore, de
     %talmanera que el vector se distribuye de manera uniforme hasta llegar
     %a cero.
    
aux1=hTorre(10)/(numCD-9);   
for i=11:numCD
   hTorre(i,:)=hTorre(i-1,:)-aux1;
end
    
%Descomposici?n de presiones fluctuantes
%armonicos=input('?En cu?ntos arm?nicos desea descomponer la respuesta? (al menos 12):   ');
armonicos=12;
k=zeros(armonicos,1);
for i=1
    k(i,:)=1;
end
for i=1:armonicos-1
    k(i+1,:)=k(i,:)+1;
end

R=3; %Arm?nico resonante
fk=zeros(armonicos,1);
fk=r./(2.^(k-R));

Tk=zeros(armonicos,1);
Tk=1./fk;

fak=zeros(armonicos,1);
fak=r./(2.^(k-0.5-R));

fpk=zeros(armonicos,1);
fpk=r./(2.^(k+0.5-R));


%Evaluando x respecto a los armónicos
xf=zeros(armonicos,1);
for i=1:armonicos
xf(i,:)=(1220*fk(i,:)/Uo);
end

xf=xf.^2;
syms fk
SrF=5953600*fk/((1+1488400*fk^2/Uo^2)^(4/3)*Uo^2)
%***********************************************************
%Integral punto por punto
for i=1: armonicos
    ckSim(i)=sqrt(2*int(SrF,fpk(i),fak(i)));
    Ck(i)=double(ckSim(i));
end

aux3=zeros(armonicos,1);
for i=1:armonicos
    aux3=aux3+Ck(i);
end
    
ck1=zeros(armonicos,1);
for i=1:armonicos
    ck1(i)=Ck(i)/aux3(i);
end

%?ngulos de fase generados pseudoaleatoriamente...
theta=(2*pi).*rand(armonicos,1);
%theta para calibrar
theta=[5.417;4.899;6.263;3.842;1.673;5.279;2.362;4.255;0.055;1.733;3.694;5.263];


%En esta parte falta corregir la manera en que recorre las "i" la función
%P, ya que primero recorre los 12 valores de los armónicos y después de
%esto deja en ceros, pues se agotan los 12 armínicos utilizados para los
%72012 valores de tiempo requeridos

%**************************************************************************
% Tiempo total de la simulación
T=6001;
P=[];
p2=zeros(T*armonicos,1);
for i=1:armonicos
    for j=1:T
        Paux(j)=cos(((pi/(Tk(R)*2.^(k(i)-R))*t(j)-theta(i))));
        p2(T*(i-1)+j)=cos(((pi/(Tk(R)*2.^(k(i)-R))*t(j)-theta(i))));
    end
P=[P Paux];
end
P=P';
% xlswrite('DatosP.xlsx',P,'Hoja1','A1');

xk=zeros(armonicos,1); 
for i=1:armonicos
xk(i)=Tk(i)/Tk(R);
end

ccR=ck1(R)/2;
ccR1=ck1(R-1)+(ck1(R)/4);
ccR11=ck1(R+1)+(ck1(R)/4);


% x=zeros(1470,1);
% for i=1
%   ff(:,i)=0.0017:0.0017:2.499;
%   x(:,i)=(1220.*ff)./Uo;
% end   
% syms ff
% Sr=zeros(1470,1);
% Sr=4.*((x.^2)./((1+x.^2).^(4/3))./ff); 

f2=zeros(armonicos,1);
for i=1:armonicos
   f2(i)=1/Tk(i);
end
sum(ck1);

%Armónicos corregidos conforme al espectro del viento (ok)
Pp=zeros(armonicos*6001,1);
for i=1:6001
    Pp(i,:)=P(i,:)*ck1(1);
end

for i=2:armonicos
for j=6001*(i-1)+1:6001*i
    Pp(j,:)=P(j,:)*ck1(i);
end
end

% xlswrite('DatosPp.xlsx',Pp,'Hoja1','A1');

%Correlación Espacial
%Tamaño de ráfaga equivalente
rafagaeq=zeros(armonicos,1);
for i=1:armonicos
    rafagaeq(i,:)=Uo/(7*f2(i,:));
end

Lk=zeros(armonicos,1);
for i=1:armonicos
    Lk(i,:)=2*rafagaeq(i,:);
end

%Posición del centro de ráfaga para más info consultar tesis Celio fontao
%Carril pp 80
Gc=h-(Lk(R)/2);

%Borrar esto posteriormente
Gc=82.6;
 
%coeficiente de reducción de presiones fluctuantes
Cred=[];
Credaux=[];
Cred1=zeros(armonicos*numCD,1);
%**************************************************************************
for j=1:armonicos
   for k=1:numCD
            if Gc<=hTorre(k) && hTorre(k)<=Gc+rafagaeq(j)
           %if round(Gc,4)<=round(hTorre(k),4) && round(hTorre(k),4)<=round(Gc+rafagaeq(j),4);
            Credaux=(1/rafagaeq(j)).*(Gc-hTorre(k))+1;
            Cred1(armonicos*(j-1)+k)=(1/rafagaeq(j)).*(Gc-hTorre(k))+1;
      elseif Gc-rafagaeq(j)<=hTorre(k) && hTorre(k)<=Gc
     %elseif round(Gc-rafagaeq(j),4)<=hTorre(k) && hTorre(k)<=Gc
                 Credaux=(-1/rafagaeq(j))*(Gc-hTorre(k))+1;
                 Cred1(armonicos*(j-1)+k)=(-1/rafagaeq(j))*(Gc-hTorre(k))+1;
            else
                Credaux=0;  
                Cred1(armonicos*(j-1)+k)=0;
            end
     Cred=[Cred Credaux];
      end
end
    Cred=Cred';
% xlswrite('Cred.xlsx',Cred,'Hoja1','A1');

%Velocidad media
Vm=[];
for i=1:numCD
    Vm(:,:)=(0.5934*Vo)*(hTorre(:,:)/10).^0.185;
end

%Velocidad Pico
V=[];
for i=1:numCD
    V(:,:)=0.94*Vo*(hTorre(:,:)/10).^0.1;
end

%Presión Pico
q=[];
for i=1:numCD
    q=0.613*(V).^2;
end

%*****Carga estática *****
qest=[];
for i=1:numCD
    qest=0.613*Vm.^2;
end
Fest=[];
Fest=CD.*areas.*qest;
%*****Carga estática *****

%Presión fluctuante
qf=[];
qf=q-qest;

size(qf)
size(Cred)
size(Pp)
%Presión armónica variable en tiempo, altura y armónicos
datos=numCD*armonicos*6001;
Q=zeros(datos,1);
% for k=1:armonicos
%     for j=1:numCD
%         for t=1:6001
%             sub1=numCD*(t-1)+armonicos*(j-1)+k;
%             sub2=armonicos*(j-1)+k;
%             sub3=armonicos*(t-1)+k;
%             Q(sub1)=qf(j)*Cred(sub2)*Pp(sub3);
%         end
%     end
% end


for t=1:6001
    for j=1:numCD
        for k=1:armonicos
            sub1=numCD*(t-1)+armonicos*(j-1)+k;
            sub2=armonicos*(j-1)+k;
            sub3=armonicos*(t-1)+k;
            Q(sub1)=qf(j)*Cred(sub2)*Pp(sub3);
        end
    end
end
NancyE=cputime-NancyT
 
    size(Q)

dlmwrite('DatosQ',Q);



NancyE=cputime-NancyT







    