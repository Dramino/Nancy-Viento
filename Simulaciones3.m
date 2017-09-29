clear,clc
clear all
format short

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
w0=0.0017;  %frecuencia inicial
dw=0.0017;  %distancia entre frecuencias
noFrecuencias=1470;
wf=w0*noFrecuencias;
x=zeros(1,noFrecuencias);
w=linspace(w0 ,wf ,noFrecuencias);
%Funci?on de densidad de turbulencia Davenport
for i=1:noFrecuencias
    x(i)=1220*dw*i/Uo;
end
Sr=4*((x.^2)./((1+x.^2).^(4/3)))./w;
  %En este punto verificar la gr?fica (espectro de Potencia) de i,Sr
%plot(w,Sr)
%grid on

%r=input('Introduzca la frecuencia natural de la estructura en ciclos/s:   ');
r=0.85;  %frecuencia de la estructura Hz
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
hTorre(1,:)=h%-(hTramo/3)
     
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
xf=(1220*fk/Uo)^2;
SrF=4*(xf/(((1+xf)^(4/3))/(fk)));


%***********************************************************
%Integral punto por punto
for i=1: armonicos
    ckSim(i)=sqrt(2*int(SrF,fpk(i),fak(i)));
    ck(i)=double(ckSim(i));
end


% %
% syms ff
% aux333=int(Sr,ff)
% %
% syms ff
% aux2=zeros(armonicos,1);
% for i=1:armonicos
% aux2(i,:)=int(Sr(i,:),ff,fpk(i,:),fak(i,:))
% end
% % for i=1:armonicos
% % aux2(i,:)=int(Sr,ff)
% % end
% Ck=zeros(armonicos,1);
% for i=1:armonicos
% Ck(1,:)=sqrt.(2*aux2(1,:))
% end
%************************************************************

Ck=[0.444865;0.560439;0.705824;0.887853;1.111508;1.366177;1.575335;1.535191;1.142286;0.672208;0.353216;0.178951];

aux3=zeros(armonicos,1);
for i=1:armonicos
    aux3=aux3+Ck(i);
end
    
ck1=zeros(armonicos,1);
for i=1:armonicos
    ck1(i)=Ck(i)/aux3(i);
end

%?ngulos de fase generados pseudoaleatoriamente...
%theta=0+(2*pi-0).*rand(armonicos,1)
theta=[5.417;4.899;6.263;3.842;1.673;5.279;2.362;4.255;0.055;1.733;3.694;5.263];

for i=1
    t(:,1)=0:0.1:600; %Tiempo
end

%En esta parte falta corregir la manera en que recorre las "i" la función
%P, ya que primero recorre los 12 valores de los armónicos y después de
%esto deja en ceros, pues se agotan los 12 armínicos utilizados para los
%72012 valores de tiempo requeridos

%**************************************************************************
%Tiempo total de la simulacion
T=6001;
P=zeros(T,armonicos);
for i=1:armonicos
	for j=1:T
		P(:,:)=cos(((pi/(Tk(R)*2.^(k(j)-R))*t(i)-theta(j))));
	end
end
xlswrite('DatosP.xlsx',P,'Hoja1','A1');
%**************************************************************************

P=xlsread('P');

xk=zeros(armonicos,1); 
for i=1:armonicos
xk(i)=Tk(i)/Tk(R);
end

ccR=ck1(R)/2;
ccR1=ck1(R-1)+(ck1(R)/4);
ccR11=ck1(R+1)+(ck1(R)/4);


%P=cos(((pi/(Tk(3)*2.^(k-R))).*t)-theta)

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

% for    i=1:armonicos
% for    j=1:armonicos*6001-(armonicos-i)*6001
%     Pp(j,:)=P(:,:)*ck1(i);
% end
% end
xlswrite('DatosPp.xlsx',Pp,'Hoja1','A1');

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
Gc=ones(armonicos*numCD)*h-(Lk(R)/2);
%Gc=h-(Lk(R)/2);

%Borrar esto posteriormente
Gc=ones(armonicos*numCD)*82.6;
 
%coeficiente de reducción de presiones fluctuantes
Cred=zeros(numCD*armonicos,1);

for i=1:numCD
    if Gc(i)<=hTorre(i) & hTorre(i)<=Gc(i)+rafagaeq(i)
         Cred(i)=(1/rafagaeq)*(Gc(i)-hTorre(i)+1)
%             else 
%             if Gc(i)-rafagaeq<=hTorre & hTorre<=Gc(i)    
%                 Cred(i)=(-1/refagaeq)*(Gc(i)-hTorre(i)+1)
%             else
%                 0
%             end
    end
end






















    