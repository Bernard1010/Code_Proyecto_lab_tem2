function [] = readTFMatRF()

clear all;
clc;
p=4001;             %Cantidad de puntos de graficas del VNA


%----MATCH1----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/match1.s2p');      %Se carga archivo s2p con los parametros S de la prueba match1(Agregar ruta de archivo)
frecs = data.Freq;
frecs1=frecs';
S1 = extract(data,'S_PARAMETERS');      %parametros S de prueba match1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----OPEN1----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/open1.s2p');      %Se carga archivo s2p con los parametros S de la prueba open1 (Agregar ruta de archivo)
frecs = data.Freq;
frecs2=frecs';
S2 = extract(data,'S_PARAMETERS');      %parametros S de prueba open1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----SHORT1----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/short1.s2p');      %Se carga archivo s2p con los parametros S de la prueba short1(Agregar ruta de archivo)
frecs = data.Freq;
frecs3=frecs';
S3 = extract(data,'S_PARAMETERS');      %parametros S de prueba short1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----MATCH2----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/match2.s2p');      %Se carga archivo s2p con los parametros S de la prueba match2(Agregar ruta de archivo)
frecs = data.Freq;
frecs4=frecs';
S4 = extract(data,'S_PARAMETERS');      %parametros S de prueba match2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----OPEN2----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/open2.s2p');      %Se carga archivo s2p con los parametros S de la prueba open2(Agregar ruta de archivo)
frecs = data.Freq;
frecs5=frecs';
S5 = extract(data,'S_PARAMETERS');      %parametros S de prueba open2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----SHORT2----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/short2.s2p');      %Se carga archivo s2p con los parametros S de la prueba short2(Agregar ruta de archivo)
frecs = data.Freq;
frecs6=frecs';
S6 = extract(data,'S_PARAMETERS');      %parametros S de prueba short2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----THRU----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/thru.s2p');      %Se carga archivo s2p con los parametros S de la prueba thru(Agregar ruta de archivo)
frecs = data.Freq;
frecs7=frecs';
S7 = extract(data,'S_PARAMETERS');      %parametros S de prueba thru

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----PARAMETROS S MEDIDOS SIN ERROR----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/thru.s2p');      %Se carga archivo s2p con los parametros S REALES SIN ERROR(VNA)(Agregar ruta de archivo)
frecs = data.Freq;
frecs=frecs';
S8 = extract(data,'S_PARAMETERS');      %parametros S de prueba thru
SP = sparameters('C:\Users\berna\Documents\Matlabs/thru.s2p');%!!CARGAR MISMA DIRECCION QUE EN LINEA 58!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INICIALIZACION DE VECTORES DE TERMINOS DE ERROR

%-----FORWARD-----%

e00 = zeros(1,p);                         %variable-vector para cargar valores de e00
e11 = zeros(1,p);                         %variable-vector para cargar valores de e11
e10e01 = zeros(1,p);                      %variable-vector para cargar valores de e10e01
e10e32 = zeros(1,p);                      %variable-vector para cargar valores de e10e32
e22 = zeros(1,p);                         %variable-vector para cargar valores de e22
e30 = 0;                                  %e30 es despreciable
%-----REVERSE-----%!!!!!!!!SE UTILIZA EL GUION BAJO PARA REPRESENTAR EL PRIMADO

e_33 = zeros(1,p);                         %variable-vector para cargar valores de e_33
e_11 = zeros(1,p);                         %variable-vector para cargar valores de e_11
e_23e_32 = zeros(1,p);                     %variable-vector para cargar valores de e_23e_32
e_23e_01 = zeros(1,p);                     %variable-vector para cargar valores de e_23e_01
e_22 = zeros(1,p);                         %variable-vector para cargar valores de e_22
e_03 = 0;                                  %e_03 es despreciable
%-----VECTORES AUXULIARES-----%

De = zeros(1,p);                            %Delta e
De_ = zeros(1,p);                           %Delta e prima


%CICLO DE CARGA DE VALORES DE TERMINOS DE ERROR
for i=1:p                                 
e00(i) = S1(1,1,i); 
e11(i) = (S2(1,1,i)+S3(1,1,i)-2*e00(i))/(S2(1,1,i)-S3(1,1,i));
e10e01(i) = -2*((S2(1,1,i)-e00(i))*(S3(1,1,i)-e00(i)))/(S2(1,1,i)-S3(1,1,i));

e_33(i) = S4(2,2,i); 
e_22(i) = (S5(2,2,i)+S6(2,2,i)-2*e_33(i))/(S5(2,2,i)-S6(2,2,i));
e_23e_32(i) = -2*((S5(2,2,i)-e_33(i))*(S6(2,2,i)-e_33(i)))/(S5(2,2,i)-S6(2,2,i));

De(i) = e00(i)*e11(i)-e10e01(i);
e22(i) = (S7(1,1,i)-e00(i))/(S7(1,1,i)*e11(i)-De(i));
e10e32(i) = (S7(2,1,i)-e30)*(1-e11(i)*e22(i)); 

De_(i) = e_33(i)*e_22(i)-e_23e_32(i);
e_11(i) = (S7(2,2,i)-e_33(i))/(S7(2,2,i)-De_(i));
e_23e_01(i) = (S7(1,2,i)-e_03)*(1-e_11(i)*e_22(i));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INICIALIZACION DE VECTORES DE PARAMETROS S SIN ERROR

a11 = zeros(1,p);
a12 = zeros(1,p);
a21 = zeros(1,p);
a22 = zeros(1,p);

S11 = zeros(1,p);
S12 = zeros(1,p);
S21 = zeros(1,p);
S22 = zeros(1,p);
D = zeros(1,p);                              %Delta s

for i=1:p
    
a11(i) = ((S8(1,1,i)-e00(i))/(e10e01(i)));
a12(i) = ((S8(1,2,i)-e_03)/(e_23e_01(i)));
a21(i) = ((S8(2,1,i)-e30)/(e10e32(i)));
a22(i) = ((S8(2,2,i)-e_33(i))/(e_23e_32(i)));

D(i) = (1+a11(i)*e11(i))*(1+a22(i)*e_22(i)-a21(i)*a12(i)*e22(i)*e_11(i));

S11(i) = (a11(i)*(1+a22(i)*e_22(i))-(e22(i)*a21(i)*a12(i)))/D(i);
S21(i) = (a21(i)*(1+a22(i)*(e_22(i)-e22(i))))/D(i);
S22(i) = (a22(i)*(1+a11(i)*e11(i))-(e_11(i)*a21(i)*a12(i)))/D(i);
S12(i) = (a12(i)*(1+a11(i)*(e11(i)-e_11(i))))/D(i);

end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GRAFICOS DE LOS PARAMETROS S

%-----------GRAFICA DE PARAMETRO S11 SIN ERROR--------%
figure(1);
plot(frecs,20*log10(abs(S11)));
grid on;   
title('Parámetro S11 con error');
xlabel('Frequency (Hz)');
ylabel('Atenuación (dB)');



%-----------GRAFICA DE PARAMETRO S21 SIN ERROR--------%
figure(2);
plot(frecs,20*log10(abs(S21)));
grid on;   
title('Parámetro S21 con error');
xlabel('Frequency (Hz)');
ylabel('Atenuación (dB)');
%axis([0 6e9 -10 1]);

%-----------GRAFICA DE PARAMETRO S22 SIN ERROR--------%
figure(3);
plot(frecs,20*log10(abs(S22)));
grid on;   
title('Parámetro S22 con error');
xlabel('Frequency (Hz)');
ylabel('Atenuación (dB)');


%-----------GRAFICA DE PARAMETRO S12 SIN ERROR--------%
figure(4);
plot(frecs,20*log10(abs(S12)));
grid on;   
title('Parámetro S12 con error');
xlabel('Frequency (Hz)');
ylabel('Atenuación (dB)');
%axis([0 6e9 -10 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------GRAFICA DE PARAMETRO S11 SIN ERROR--------%
figure(5);                                  
rfplot(SP,1,1);
title('Parametro S11');

%-----------GRAFICA DE PARAMETRO S21 SIN ERROR--------%
figure(6);                                  
rfplot(SP,2,1);
title('Parametro S21');

%-----------GRAFICA DE PARAMETRO S22 SIN ERROR--------%
figure(7);                                  
rfplot(SP,2,2);
title('Parametro S22');

%-----------GRAFICA DE PARAMETRO S12 SIN ERROR--------%
figure(8);                                  
rfplot(SP,1,2);
title('Parametro S12');



% figure(1)
% A1 = zeros(1,4001);                         
% for i=1:4001                                
% A1(i)=20*log(abs(e_23e_01(i))); 
% end
% plot(frecs,A1);
% 
% grid on;
% title('Término e´23e´01');
% xlabel('Frequency (Hz)');
% ylabel('Ganancia (dB)');
% 
% 
% A0 = zeros(1,4001);                        
% for i=1:4001                               
% A0(i)=angle(e_23e_01(i)); 
% end
% 
% lonWrapped = wrapTo2Pi(A0);
% 
% figure(2)
% plot(frecs,lonWrapped);
% title('Fase de término e´23´01');
% xlabel('Frequency (Hz)');
% ylabel('Fase (rad)');


end




