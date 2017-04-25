function [] = readTFMatRF()

clear all;
clc;
p=4001;             %Cantidad de puntos de graficas del VNA


%----MATCH1----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/snake.s2p');      %Se carga archivo s2p con los parametros S de la prueba match1(Agregar ruta de archivo)
frecs = data.Freq;
frecs1=frecs';
S1 = extract(data,'S_PARAMETERS');      %parametros S de prueba match1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----OPEN1----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/L50_1.s2p');      %Se carga archivo s2p con los parametros S de la prueba open1 (Agregar ruta de archivo)
frecs = data.Freq;
frecs1=frecs';
S2 = extract(data,'S_PARAMETERS');      %parametros S de prueba open1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----SHORT1----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/L100.s2p');      %Se carga archivo s2p con los parametros S de la prueba short1(Agregar ruta de archivo)
frecs = data.Freq;
frecs1=frecs';
S3 = extract(data,'S_PARAMETERS');      %parametros S de prueba short1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----MATCH2----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/L25.s2p');      %Se carga archivo s2p con los parametros S de la prueba match2(Agregar ruta de archivo)
frecs = data.Freq;
frecs1=frecs';
S4 = extract(data,'S_PARAMETERS');      %parametros S de prueba match2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----OPEN2----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/L50_1.s2p');      %Se carga archivo s2p con los parametros S de la prueba open2(Agregar ruta de archivo)
frecs = data.Freq;
frecs1=frecs';
S5 = extract(data,'S_PARAMETERS');      %parametros S de prueba open2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----SHORT2----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/L100.s2p');      %Se carga archivo s2p con los parametros S de la prueba short2(Agregar ruta de archivo)
frecs = data.Freq;
frecs1=frecs';
S6 = extract(data,'S_PARAMETERS');      %parametros S de prueba short2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----THRU----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/L100.s2p');      %Se carga archivo s2p con los parametros S de la prueba thru(Agregar ruta de archivo)
frecs = data.Freq;
frecs1=frecs';
S7 = extract(data,'S_PARAMETERS');      %parametros S de prueba thru

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----PARAMETROS S MEDIDOS SIN ERROR----%
data = read(rfdata.data,'C:\Users\berna\Documents\Matlabs/L100.s2p');      %Se carga archivo s2p con los parametros S REALES SIN ERROR(VNA)(Agregar ruta de archivo)
frecs = data.Freq;
frecs1=frecs';
S8 = extract(data,'S_PARAMETERS');      %parametros S de prueba thru

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

De(i) = e00(i)*e11(i)-e01e10;
e22(i) = (S7(1,1,i)-e00(i))/(S7(1,1,i)*e11(i)-De(i));
e10e32(i) = (S7(2,1,i)-e30)*(1-e11(i)*e22(i)); 

De_(i) = e_33(i)*e_22(i)-e_23e_32(i);
e_11(i) = (S7(2,2,i)-e_33(i))/(S7(2,2,i)-De_(i));
e_23e_01(i) = (S7(1,2,i)-e_03)*(1-e_11*e_22);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INICIALIZACION DE VECTORES DE PARAMETROS S CON ERROR

S11m = zeros(1,p);
S12m = zeros(1,p);
S21m = zeros(1,p);
S22m = zeros(1,p);
Ds = zeros(1,p);                              %Delta s

for i=1:p
    
Ds(i) = S8(1,1,i)*S8(2,2,i)-S8(1,2,i)*S8(2,1,i);
S11m(i) = e00(i)+ (e10e01(i)*(S8(1,1,i)-e22(i)*Ds(i)))/(1-e11(i)*S8(1,1,i)-e22(i)*S8(2,2,i)+e11(i)*e22(2)*Ds(i));
S21m(i) = e30+ (e10e32(i)*S8(2,1,i))/(1-e11(i)*S8(1,1,i)-e22(i)*S8(2,2,i)+e11(i)*e22(2)*Ds(i));     
S22m(i) = e_33(i)+ (e_23e_32(i)*(S8(2,2,i)-e_11(i)*Ds(i)))/(1-e_11(i)*S8(1,1,i)-e_22(i)*S8(2,2,i)+e_11(i)*e_22(2)*Ds(i));
S12m(i) = e_03+ (e_23e_01(i)*S8(2,1,i))/(1-e_11(i)*S8(1,1,i)-e_22(i)*S8(2,2,i)+e_11(i)*e_22(2)*Ds(i));

end 


%plot(frecs1,abs(A0));



end

