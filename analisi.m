%% Anàlisi de dades

set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

%esta en pixels passar a um 0.65 (10x)
%he hecho la estadistica con 25 timepoints en todos, con f0002 que?, falta
%7 que lo tengo en el lab
At=5;
%{
for j=1:1:7
    file=['f000', num2str(j), '.mat'];
    load(file);
    if j==1
        m=(matpos(7,23)-matpos(7,19))/4;
        n=matpos(7,23)-m*23;
        for i=20:1:22
            matpos(7,i)=m*i+n;
        end
    elseif j==2
        m=(matpos(1,10)-matpos(1,7))/3;
        n=matpos(1,10)-m*10;
        for i=8:9
            matpos(1,i)=m*i+n;
        end
        m=(matpos(8,20)-matpos(8,18))/2;
        n=matpos(8,20)-m*20;
        matpos(1,19)=m*19+n;
    elseif j==3
        m=(matpos(3,20)-matpos(3,18))/2;
        n=matpos(3,20)-m*20;
        matpos(3,19)=m*19+n;
        matpos(7,1)=matpos(7,3)+6;
        matpos(7,2)=matpos(7,3)+3;
        m=(matpos(8,20)-matpos(8,18))/2;
        n=matpos(8,20)-m*20;
        matpos(8,19)=m*19+n;
    elseif j==4
        m=(matpos(6,6)-matpos(6,4))/2;
        n=matpos(6,6)-m*6;
        matpos(6,5)=m*5+n;
        m=(matpos(7,8)-matpos(7,6))/2;
        n=matpos(7,6)-m*6;
        matpos(7,7)=m*7+n;
        m=(matpos(7,21)-matpos(7,19))/2;
        n=matpos(7,21)-m*21;
        matpos(7,20)=m*20+n;
        u=zeros(8,1);
        u(:,1)=matpos(:,1)-(matpos(:,2)-matpos(:,1));
        matpos=[u, matpos];
    elseif j==5
        matpos(7,1)=matpos(7,3)+6;
        matpos(7,2)=matpos(7,3)+3;
        m=(matpos(1,16)-matpos(1,10))/6;
        n=matpos(1,10)-m*10;
        for i=11:1:15
            matpos(1,i)=m*i+n;
        end
    elseif j==6
        matpos(6,1)=matpos(6,2)-2;
    elseif j==7
        m=(matpos(5,23)-matpos(5,22));
        n=matpos(5,22)-m*22;
        for i=24:1:25
            matpos(5,i)=m*i+n;
        end
    end
    %col 1
    velizq1=((matpos(1,2:end)-matpos(1,1:end-1))/At)*0.65; meanvizq1=mean(velizq1);  stdvizq1=mean(velizq1);
    velup1=((matpos(2,2:end)-matpos(2,1:end-1))/At)*0.65; meanvup1=mean(velup1); stdvup1=std(velup1);
    velder1=((matpos(3,1:end-1)-matpos(3,2:end))/At)*0.65; meanvder1=mean(velder1); stdvder1=std(velder1);
    veldown1=((matpos(4,1:end-1)-matpos(4,2:end))/At)*0.65; meanvdown1=mean(veldown1);  stdvdown1=mean(veldown1);
    %col 2
    velizq2=((matpos(5,2:end)-matpos(5,1:end-1))/At)*0.65; meanvizq2=mean(velizq2); stdvizq2=mean(velizq2);
    velup2=((matpos(6,2:end)-matpos(6,1:end-1))/At)*0.65; meanvup2=mean(velup2); stdvup2=std(velup2);
    velder2=((matpos(7,1:end-1)-matpos(7,2:end))/At)*0.65; meanvder2=mean(velder2);  stdvder2=std(velder2);
    veldown2=((matpos(8,1:end-1)-matpos(8,2:end))/At)*0.65; meanvdown2=mean(veldown2); stdvdown2=mean(veldown2);
    
    if j==1
        xlen1=matpos(3,1)-matpos(1,1); xfin1=matpos(3,end-1)-matpos(1,end-1);
        ylen1=matpos(4,1)-matpos(2,1); yfin1=matpos(4,end-1)-matpos(2,end-1);
        xlen2=matpos(3,1)-matpos(1,1); xfin2=matpos(3,end-1)-matpos(1,end-1);
        ylen2=matpos(4,1)-matpos(2,1); yfin2=matpos(4,end-1)-matpos(2,end-1);
    else
        xlen1=matpos(3,1)-matpos(1,1); xfin1=matpos(3,end)-matpos(1,end);
        ylen1=matpos(4,1)-matpos(2,1); yfin1=matpos(4,end)-matpos(2,end);
        xlen2=matpos(3,1)-matpos(1,1); xfin2=matpos(3,end)-matpos(1,end);
        ylen2=matpos(4,1)-matpos(2,1); yfin2=matpos(4,end)-matpos(2,end);
    end

    if j==1
        shortx=[meanvder1, meanvizq2, meanvder2];
        shorty=[meanvup1, meanvup2, meanvdown2];
        vxtshort=[velder1(1:end-1); velizq2(1:end-1); velder2(1:end-1)];
        vytshort=[velup1(1:end-1); velup2(1:end-1); veldown2(1:end-1)];
        stdxshort=[stdvder1(1:end-1) stdvizq2(1:end-1) stdvder2(1:end-1)];
        stdyshort=[stdvup1(1:end-1) stdvup2(1:end-1) stdvdown2(1:end-1)];
    elseif j==2
    % solo va hasta 20
    elseif 
%}
                


%% F0001
load('f0001.mat'); %matpos izq up der down
vec=1:1:length(matpos); K=length(matpos);
figure(1)
%plot(vec,matpos(1,:),'-o','DisplayName','izq')
hold on
 plot(vec,matpos(2,:),'-o','DisplayName','up')
% hold on
% plot(vec,matpos(3,:),'-o','DisplayName','der')
% hold on
% plot(vec,matpos(4,:),'-o','DisplayName','down')
legend
title('Colony one')

figure(2)
%plot(vec,matpos(5,:),'-o','DisplayName','izq')
% hold on
% plot(vec,matpos(6,:),'-o','DisplayName','up')
% hold on
 plot(vec,matpos(7,:),'-o','DisplayName','der')
% hold on
% plot(vec,matpos(8,:),'-o','DisplayName','down')
legend
title('Colony two')

%Eliminamos 1down y 1izq, y a 2der, los puntos 19 y 23 estan bien, los del
%medio no, interpolacion lineal

m=(matpos(7,23)-matpos(7,19))/4;
n=matpos(7,23)-m*23;
for i=20:1:22
    matpos(7,i)=m*i+n;
end

%col 1
%velizq1=((matpos(1,2:end)-matpos(1,1:end-1))/At)*0.65; meanvizq1=mean(velizq1);  stdvizq1=mean(velizq1);
%velup1=((matpos(2,2:end)-matpos(2,1:end-1))/At)*0.65; meanvup1=mean(velup1); stdvup1=std(velup1);
velder1=((matpos(3,1:end-1)-matpos(3,2:end))/At)*0.65; meanvder1=mean(velder1); stdvder1=std(velder1);
%veldown1=((matpos(4,1:end-1)-matpos(4,2:end))/At)*0.65; meanvdown1=mean(veldown1);  stdvdown1=mean(veldown1);
%col 2
velizq2=((matpos(5,2:end)-matpos(5,1:end-1))/At)*0.65; meanvizq2=mean(velizq2); stdvizq2=mean(velizq2);
velup2=((matpos(6,2:end)-matpos(6,1:end-1))/At)*0.65; meanvup2=mean(velup2); stdvup2=std(velup2);
velder2=((matpos(7,1:end-1)-matpos(7,2:end))/At)*0.65; meanvder2=mean(velder2);  stdvder2=std(velder2);
veldown2=((matpos(8,1:end-1)-matpos(8,2:end))/At)*0.65; meanvdown2=mean(veldown2); stdvdown2=mean(veldown2);

xlen11=matpos(3,1)-matpos(1,1); xfin11=matpos(3,end-1)-matpos(1,end-1);
ylen11=matpos(4,1)-matpos(2,1); yfin11=matpos(4,end-1)-matpos(2,end-1);
xlen12=matpos(7,1)-matpos(5,1); xfin12=matpos(7,end-1)-matpos(5,end-1); 
ylen12=matpos(8,1)-matpos(6,1); yfin12=matpos(8,end-1)-matpos(6,end-1); 
x1=[xlen11 xfin11; xlen12 xfin12];  y1=[ylen11 yfin11; ylen12 yfin12]; 
 
 shortx=[meanvder1, meanvizq2, meanvder2];
 shorty=[meanvup2, meanvdown2];
 vxtshort=[velder1(1:end-1); velizq2(1:end-1); velder2(1:end-1)];
vytshort=[velup2(1:end-1); veldown2(1:end-1)];
stdxshort=[stdvder1(1:end-1) stdvizq2(1:end-1) stdvder2(1:end-1)];
stdyshort=[stdvup2(1:end-1) stdvdown2(1:end-1)];


% %% F0002
% close all
% load('f0002.mat'); %matpos izq up der down
% vec=1:1:length(matpos); K=length(matpos);
% figure(1)
% %plot(vec,matpos(1,:),'-o','DisplayName','izq')
% hold on
% % plot(vec,matpos(2,:),'-o','DisplayName','up')
% % hold on
% % plot(vec,matpos(3,:),'-o','DisplayName','der')
% % hold on
%  plot(vec,matpos(4,:),'-o','DisplayName','down')
% legend
% title('Colony one')
% 
% figure(2)
% %plot(vec,matpos(5,:),'-o','DisplayName','izq')
% % hold on
%  %plot(vec,matpos(6,:),'-o','DisplayName','up')
% % hold on
% %plot(vec,matpos(7,:),'-o','DisplayName','der')
% % hold on
%  plot(vec,matpos(8,:),'-o','DisplayName','down')
% legend
% title('Colony two')
% 
% % izq 1 interpolacio, down2 punto 19
% 
% m=(matpos(1,10)-matpos(1,7))/3;
% n=matpos(1,10)-m*10;
% for i=8:9
%     matpos(1,i)=m*i+n;
% end
% 
% m=(matpos(8,20)-matpos(8,18))/2;
% n=matpos(8,20)-m*20;
% matpos(1,19)=m*19+n;
% 
% 
% 
% 
% %col 1
% velizq1=((matpos(1,2:end)-matpos(1,1:end-1))/At)*0.65; meanvizq1=mean(velizq1);  stdvizq1=mean(velizq1);
% velup1=((matpos(2,2:end)-matpos(2,1:end-1))/At)*0.65; meanvup1=mean(velup1); stdvup1=std(velup1);
% velder1=((matpos(3,1:end-1)-matpos(3,2:end))/At)*0.65; meanvder1=mean(velder1); stdvder1=std(velder1);
% veldown1=((matpos(4,1:end-1)-matpos(4,2:end))/At)*0.65; meanvdown1=mean(veldown1);  stdvdown1=mean(veldown1);
% %col 2
% velizq2=((matpos(5,2:end)-matpos(5,1:end-1))/At)*0.65; meanvizq2=mean(velizq2); stdvizq2=mean(velizq2);
% velup2=((matpos(6,2:end)-matpos(6,1:end-1))/At)*0.65; meanvup2=mean(velup2); stdvup2=std(velup2);
% velder2=((matpos(7,1:end-1)-matpos(7,2:end))/At)*0.65; meanvder2=mean(velder2);  stdvder2=std(velder2);
% veldown2=((matpos(8,1:end-1)-matpos(8,2:end))/At)*0.65; meanvdown2=mean(veldown2); stdvdown2=mean(veldown2);
% 
% xlen11=matpos(3,1)-matpos(1,1); xfin11=matpos(3,end)-matpos(1,end);
% ylen11=matpos(4,1)-matpos(2,1); yfin11=matpos(4,end)-matpos(2,end);
% xlen12=matpos(7,1)-matpos(5,1); xfin12=matpos(7,end)-matpos(5,end); 
% ylen12=matpos(8,1)-matpos(6,1); yfin12=matpos(8,end)-matpos(6,end); 
%x1=[x1; xlen11 xfin11; xlen12 xfin12];  y1=[y1; ylen11 yfin11; ylen12 yfin12]; 
 
 
%  shortx=[shortx, meanvizq1, meanvder1, meanvizq2, meanvder2];
%  shorty=[shorty, meanvup1, meanvdown1, meanvup2, meanvdown2];
%  vytshort=[vytshort; velup1(1:end-1); veldown1(1:end-1); velup2(1:end-1); veldown2(1:end-1)];
% vxtshort=[vxtshort; velizq1(1:end-1); velder1(1:end-1); velizq2(1:end-1); velder2(1:end-1)];
% stdxshort=[stdxshort stdvizq1(1:end-1 ) stdvder1(1:end-1) stdvizq2(1:end-1) stdvder2(1:end-1)];
% stdyshort=[stdyshort stdvup1(1:end-1) stdvdown1(1:end-1) stdvup2(1:end-1) stdvdown2(1:end-1)];


%% F0003
close all
load('f0003.mat'); %matpos izq up der down
vec=1:1:length(matpos); K=length(matpos);
figure(1)
%plot(vec,matpos(1,:),'-o','DisplayName','izq')
hold on
% plot(vec,matpos(2,:),'-o','DisplayName','up')
% hold on
 plot(vec,matpos(3,:),'-o','DisplayName','der')
% hold on
 %plot(vec,matpos(4,:),'-o','DisplayName','down')
legend
title('Colony one')

figure(2)
%plot(vec,matpos(5,:),'-o','DisplayName','izq')
% hold on
% plot(vec,matpos(6,:),'-o','DisplayName','up')
% hold on
%plot(vec,matpos(7,:),'-o','DisplayName','der')
% hold on
 plot(vec,matpos(8,:),'-o','DisplayName','down')
legend
title('Colony two')

% 1 der interpolacio lineal punt 19, 2 der primers dos punts interpolar,
% 2 izq eliminar, 2 down 19 filtrar


m=(matpos(3,20)-matpos(3,18))/2;
n=matpos(3,20)-m*20;
matpos(3,19)=m*19+n;
matpos(7,1)=matpos(7,3)+6;
matpos(7,2)=matpos(7,3)+3;
m=(matpos(8,20)-matpos(8,18))/2;
n=matpos(8,20)-m*20;
matpos(8,19)=m*19+n;



%col 1
velizq1=((matpos(1,2:end)-matpos(1,1:end-1))/At)*0.65; meanvizq1=mean(velizq1);  stdvizq1=mean(velizq1);
velup1=((matpos(2,2:end)-matpos(2,1:end-1))/At)*0.65; meanvup1=mean(velup1); stdvup1=std(velup1);
velder1=((matpos(3,1:end-1)-matpos(3,2:end))/At)*0.65; meanvder1=mean(velder1); stdvder1=std(velder1);
veldown1=((matpos(4,1:end-1)-matpos(4,2:end))/At)*0.65; meanvdown1=mean(veldown1);  stdvdown1=mean(veldown1);
%col 2
%velizq2=((matpos(5,2:end)-matpos(5,1:end-1))/At)*0.65; meanvizq2=mean(velizq2); stdvizq2=mean(velizq2);
velup2=((matpos(6,2:end)-matpos(6,1:end-1))/At)*0.65; meanvup2=mean(velup2); stdvup2=std(velup2);
velder2=((matpos(7,1:end-1)-matpos(7,2:end))/At)*0.65; meanvder2=mean(velder2);  stdvder2=std(velder2);
veldown2=((matpos(8,1:end-1)-matpos(8,2:end))/At)*0.65; meanvdown2=mean(veldown2); stdvdown2=mean(veldown2);

xlen11=matpos(3,1)-matpos(1,1); xfin11=matpos(3,end)-matpos(1,end);
ylen11=matpos(4,1)-matpos(2,1); yfin11=matpos(4,end)-matpos(2,end);
xlen12=matpos(7,1)-matpos(5,1); xfin12=matpos(7,end)-matpos(5,end); 
ylen12=matpos(8,1)-matpos(6,1); yfin12=matpos(8,end)-matpos(6,end); 
x1=[x1; xlen11 xfin11; xlen12 xfin12];  y1=[y1; ylen11 yfin11; ylen12 yfin12]; 
 

 shortx=[shortx, meanvizq1, meanvder1, meanvder2];
 shorty=[shorty, meanvup1, meanvdown1, meanvup2, meanvdown2];
 vytshort=[vytshort; velup1; veldown1; velup2; veldown2];
vxtshort=[vxtshort; velizq1; velder1; velder2];
stdxshort=[stdxshort stdvizq1 stdvder1 stdvder2];
stdyshort=[stdyshort stdvup1 stdvdown1  stdvup2 stdvdown2];

%% F0004
load('f0004.mat'); %matpos izq up der down
vec=1:1:length(matpos); K=length(matpos);
figure(1)
%plot(vec,matpos(1,:),'-o','DisplayName','izq')
hold on
 %plot(vec,matpos(2,:),'-o','DisplayName','up')
% hold on
%plot(vec,matpos(3,:),'-o','DisplayName','der')
% hold on
 plot(vec,matpos(4,:),'-o','DisplayName','down')
legend
title('Colony one')

figure(2)
%plot(vec,matpos(5,:),'-o','DisplayName','izq')
% hold on
% plot(vec,matpos(6,:),'-o','DisplayName','up')
% hold on
% plot(vec,matpos(7,:),'-o','DisplayName','der')
% hold on
 plot(vec,matpos(8,:),'-o','DisplayName','down')
legend
title('Colony two')

%Añadir primer punto, up2 un punt, der2 dos punts, 

m=(matpos(6,6)-matpos(6,4))/2;
n=matpos(6,6)-m*6;
matpos(6,5)=m*5+n;
m=(matpos(7,8)-matpos(7,6))/2;
n=matpos(7,6)-m*6;
matpos(7,7)=m*7+n;
m=(matpos(7,21)-matpos(7,19))/2;
n=matpos(7,21)-m*21;
matpos(7,20)=m*20+n;
u=zeros(8,1);
u(:,1)=matpos(:,1)-(matpos(:,2)-matpos(:,1));
matpos=[u, matpos];


%col 1
velizq1=((matpos(1,2:end)-matpos(1,1:end-1))/At)*0.65; meanvizq1=mean(velizq1);  stdvizq1=mean(velizq1);
velup1=((matpos(2,2:end)-matpos(2,1:end-1))/At)*0.65; meanvup1=mean(velup1); stdvup1=std(velup1);
velder1=((matpos(3,1:end-1)-matpos(3,2:end))/At)*0.65; meanvder1=mean(velder1); stdvder1=std(velder1);
veldown1=((matpos(4,1:end-1)-matpos(4,2:end))/At)*0.65; meanvdown1=mean(veldown1);  stdvdown1=mean(veldown1);
%col 2
velizq2=((matpos(5,2:end)-matpos(5,1:end-1))/At)*0.65; meanvizq2=mean(velizq2); stdvizq2=mean(velizq2);
velup2=((matpos(6,2:end)-matpos(6,1:end-1))/At)*0.65; meanvup2=mean(velup2); stdvup2=std(velup2);
velder2=((matpos(7,1:end-1)-matpos(7,2:end))/At)*0.65; meanvder2=mean(velder2);  stdvder2=std(velder2);
veldown2=((matpos(8,1:end-1)-matpos(8,2:end))/At)*0.65; meanvdown2=mean(veldown2); stdvdown2=mean(veldown2);

xlen11=matpos(3,1)-matpos(1,1); xfin11=matpos(3,end)-matpos(1,end);
ylen11=matpos(4,1)-matpos(2,1); yfin11=matpos(4,end)-matpos(2,end);
xlen12=matpos(7,1)-matpos(5,1); xfin12=matpos(7,end)-matpos(5,end); 
ylen12=matpos(8,1)-matpos(6,1); yfin12=matpos(8,end)-matpos(6,end); 
x1=[x1; xlen11 xfin11; xlen12 xfin12];  y1=[y1; ylen11 yfin11; ylen12 yfin12]; 
 


 longx=[meanvizq1, meanvder1, meanvizq2, meanvder2];
 longy=[meanvup1, meanvdown1, meanvup2, meanvdown2];
vytlong=[velup1; veldown1; velup2; veldown2];
vxtlong=[velizq1; velizq2; velder2];
stdxlong=[stdvizq1 stdvizq2 stdvder2];
stdylong=[stdvup1 stdvdown1  stdvup2 stdvdown2];


%% F0005
close all
load('f0005.mat'); %matpos izq up der down
vec=1:1:length(matpos); K=length(matpos);
figure(1)
plot(vec,matpos(1,:),'-o','DisplayName','izq')
hold on
% plot(vec,matpos(2,:),'-o','DisplayName','up')
% hold on
 %plot(vec,matpos(3,:),'-o','DisplayName','der')
% hold on
 %plot(vec,matpos(4,:),'-o','DisplayName','down')
legend
title('Colony one')

figure(2)
plot(vec,matpos(5,:),'-o','DisplayName','izq')
% hold on
%plot(vec,matpos(6,:),'-o','DisplayName','up')
% hold on
 %plot(vec,matpos(7,:),'-o','DisplayName','der')
% hold on
% plot(vec,matpos(8,:),'-o','DisplayName','down')
legend
title('Colony two')

%filtrar primer dos punts de der2, 1 izq (10,16 mig interpolacio)

matpos(7,1)=matpos(7,3)+6;
matpos(7,2)=matpos(7,3)+3;
m=(matpos(1,16)-matpos(1,10))/6;
n=matpos(1,10)-m*10;
for i=11:1:15
    matpos(1,i)=m*i+n;
end

m=(matpos(5,23)-matpos(5,21))/2;
n=matpos(5,21)-m*21;
matpos(5,22)=m*22+n;


%col 1
velizq1=((matpos(1,2:end)-matpos(1,1:end-1))/At)*0.65; meanvizq1=mean(velizq1);  stdvizq1=mean(velizq1);
velup1=((matpos(2,2:end)-matpos(2,1:end-1))/At)*0.65; meanvup1=mean(velup1); stdvup1=std(velup1);
velder1=((matpos(3,1:end-1)-matpos(3,2:end))/At)*0.65; meanvder1=mean(velder1); stdvder1=std(velder1);
veldown1=((matpos(4,1:end-1)-matpos(4,2:end))/At)*0.65; meanvdown1=mean(veldown1);  stdvdown1=mean(veldown1);
%col 2
velizq2=((matpos(5,2:end)-matpos(5,1:end-1))/At)*0.65; meanvizq2=mean(velizq2); stdvizq2=mean(velizq2);
velup2=((matpos(6,2:end)-matpos(6,1:end-1))/At)*0.65; meanvup2=mean(velup2); stdvup2=std(velup2);
velder2=((matpos(7,1:end-1)-matpos(7,2:end))/At)*0.65; meanvder2=mean(velder2);  stdvder2=std(velder2);
veldown2=((matpos(8,1:end-1)-matpos(8,2:end))/At)*0.65; meanvdown2=mean(veldown2); stdvdown2=mean(veldown2);

xlen11=matpos(3,1)-matpos(1,1); xfin11=matpos(3,end)-matpos(1,end);
ylen11=matpos(4,1)-matpos(2,1); yfin11=matpos(4,end)-matpos(2,end);
xlen12=matpos(7,1)-matpos(5,1); xfin12=matpos(7,end)-matpos(5,end); 
ylen12=matpos(8,1)-matpos(6,1); yfin12=matpos(8,end)-matpos(6,end); 
x1=[x1; xlen11 xfin11; xlen12 xfin12];  y1=[y1; ylen11 yfin11; ylen12 yfin12]; 
 

 longx=[longx, meanvizq1, meanvder1, meanvizq2, meanvder2];
 longy=[longy, meanvup1, meanvdown1, meanvup2, meanvdown2];
vytlong=[vytlong; velup1; veldown1; velup2; veldown2];
vxtlong=[vxtlong; velizq1; velder1; velder2];
stdxlong=[stdxlong stdvizq1 stdvder1 stdvder2];
stdylong=[stdylong stdvup1 stdvdown1  stdvup2 stdvdown2];


%% F0006
close all
load('f0006.mat'); %matpos izq up der down
vec=1:1:length(matpos); K=length(matpos);
figure(1)
%plot(vec,matpos(1,:),'-o','DisplayName','izq')
hold on
% plot(vec,matpos(2,:),'-o','DisplayName','up')
% hold on
 plot(vec,matpos(3,:),'-o','DisplayName','der')
% hold on
% plot(vec,matpos(4,:),'-o','DisplayName','down')
legend
title('Colony one')

figure(2)
%plot(vec,matpos(5,:),'-o','DisplayName','izq')
% hold on
 plot(vec,matpos(6,:),'-o','DisplayName','up')
% hold on
%plot(vec,matpos(7,:),'-o','DisplayName','der')
% hold on
% plot(vec,matpos(8,:),'-o','DisplayName','down')
legend
title('Colony two')

%2 izq eliminar y 2 up 1er punto
matpos(6,1)=matpos(6,2)-2;

%col 1
velizq1=((matpos(1,2:end)-matpos(1,1:end-1))/At)*0.65; meanvizq1=mean(velizq1);  stdvizq1=mean(velizq1);
velup1=((matpos(2,2:end)-matpos(2,1:end-1))/At)*0.65; meanvup1=mean(velup1); stdvup1=std(velup1);
%velder1=((matpos(3,1:end-1)-matpos(3,2:end))/At)*0.65; meanvder1=mean(velder1); stdvder1=std(velder1);
veldown1=((matpos(4,1:end-1)-matpos(4,2:end))/At)*0.65; meanvdown1=mean(veldown1);  stdvdown1=mean(veldown1);
%col 2
%velizq2=((matpos(5,2:end)-matpos(5,1:end-1))/At)*0.65; meanvizq2=mean(velizq2); stdvizq2=mean(velizq2);
velup2=((matpos(6,2:end)-matpos(6,1:end-1))/At)*0.65; meanvup2=mean(velup2); stdvup2=std(velup2);
velder2=((matpos(7,1:end-1)-matpos(7,2:end))/At)*0.65; meanvder2=mean(velder2);  stdvder2=std(velder2);
veldown2=((matpos(8,1:end-1)-matpos(8,2:end))/At)*0.65; meanvdown2=mean(veldown2); stdvdown2=mean(veldown2);

xlen11=matpos(3,1)-matpos(1,1); xfin11=matpos(3,end)-matpos(1,end);
ylen11=matpos(4,1)-matpos(2,1); yfin11=matpos(4,end)-matpos(2,end);
xlen12=matpos(7,1)-matpos(5,1); xfin12=matpos(7,end)-matpos(5,end); 
ylen12=matpos(8,1)-matpos(6,1); yfin12=matpos(8,end)-matpos(6,end); 
x1=[x1; xlen11 xfin11; xlen12 xfin12];  y1=[y1; ylen11 yfin11; ylen12 yfin12]; 
 

longx=[longx, meanvizq1, meanvder1, meanvder2];
 longy=[longy, meanvup1, meanvdown1, meanvup2, meanvdown2];
vytlong=[vytlong; velup1; veldown1; velup2; veldown2];
vxtlong=[vxtlong; velizq1; velder2];
stdxlong=[stdxlong stdvizq1 stdvder2];
stdylong=[stdylong stdvup1 stdvdown1  stdvup2 stdvdown2];


%% F0007
close all
load('f0007.mat'); %matpos izq up der down
vec=1:1:length(matpos); K=length(matpos);
figure(1)
%plot(vec,matpos(1,:),'-o','DisplayName','izq')
hold on
% plot(vec,matpos(2,:),'-o','DisplayName','up')
% hold on
% plot(vec,matpos(3,:),'-o','DisplayName','der')
% hold on
 plot(vec,matpos(4,:),'-o','DisplayName','down')
legend
title('Colony one')

figure(2)
%plot(vec,matpos(5,:),'-o','DisplayName','izq')
% hold on
% plot(vec,matpos(6,:),'-o','DisplayName','up')
% hold on
%plot(vec,matpos(7,:),'-o','DisplayName','der')
% hold on
 plot(vec,matpos(8,:),'-o','DisplayName','down')
legend
title('Colony two')

%izq2 last two points, 2up fora (grafica no dolenta, pero gif si),
%down2 s'hauria de fer alguna cosa

m=(matpos(5,23)-matpos(5,22));
n=matpos(5,22)-m*22;
for i=24:1:25
    matpos(5,i)=m*i+n;
end

%col 1
velizq1=((matpos(1,2:end)-matpos(1,1:end-1))/At)*0.65; meanvizq1=mean(velizq1);  stdvizq1=mean(velizq1);
velup1=((matpos(2,2:end)-matpos(2,1:end-1))/At)*0.65; meanvup1=mean(velup1); stdvup1=std(velup1);
velder1=((matpos(3,1:end-1)-matpos(3,2:end))/At)*0.65; meanvder1=mean(velder1); stdvder1=std(velder1);
veldown1=((matpos(4,1:end-1)-matpos(4,2:end))/At)*0.65; meanvdown1=mean(veldown1);  stdvdown1=mean(veldown1);
%col 2
velizq2=((matpos(5,2:end)-matpos(5,1:end-1))/At)*0.65; meanvizq2=mean(velizq2); stdvizq2=mean(velizq2);
%velup2=((matpos(6,2:end)-matpos(6,1:end-1))/At)*0.65; meanvup2=mean(velup2); stdvup2=std(velup2);
%velder2=((matpos(7,1:end-1)-matpos(7,2:end))/At)*0.65; meanvder2=mean(velder2);  stdvder2=std(velder2);
%veldown2=((matpos(8,1:end-1)-matpos(8,2:end))/At)*0.65; meanvdown2=mean(veldown2); stdvdown2=mean(veldown2);

xlen11=matpos(3,1)-matpos(1,1); xfin11=matpos(3,end)-matpos(1,end);
ylen11=matpos(4,1)-matpos(2,1); yfin11=matpos(4,end)-matpos(2,end);
xlen12=matpos(7,1)-matpos(5,1); xfin12=matpos(7,end)-matpos(5,end); 
ylen12=matpos(8,1)-matpos(6,1); yfin12=matpos(8,end)-matpos(6,end); 
x1=[x1; xlen11 xfin11; xlen12 xfin12];  y1=[y1; ylen11 yfin11; ylen12 yfin12]; 
 

longx=[longx, meanvizq1, meanvder1, meanvizq2];
 longy=[longy, meanvup1, meanvdown1, meanvdown2];
vytlong=[vytlong; velup1; veldown1];
vxtlong=[vxtlong; velizq1; velder1; velizq2];
stdxlong=[stdxlong stdvizq1 stdvder1 stdvizq2];
stdylong=[stdylong stdvup1 stdvdown1];




%% Boxplot

figure(166)
x = [shortx longx];
g1 = repmat({'1/3'},1,length(shortx));
g2 = repmat({'1/5'},1,length(longx));
g = [g1 g2];
boxplot(x,g)
title('Mean x velocity')
ylabel('$$\langle v_{x} \rangle$$ ($\mu$m/min)')
xlabel('Width to length colony proportion')

figure(266)
x = [shorty longy];
g1 = repmat({'1/3'},1,length(shorty));
g2 = repmat({'1/5'},1,length(longy));
g = [g1 g2];
boxplot(x,g)
title('Mean y velocity')
ylabel('$$\langle v_{y} \rangle$$  ($\mu$m/min)')
xlabel('Width to length colony proportion')
%%
figure(6)
plot(ones(1,length(shortx)),shortx,'bo')
hold on 
plot(2*ones(1,length(longx)),longx,'ro')
title('Mean x velocities')
ylabel('$$\langle v_{x} \rangle$$ ($\mu$m/min)')
xlabel('Width to length colony proportion')
xlim([0 3])
xticks([1 2])
xticklabels({'1/3','1/5'})

figure(7)
plot(ones(1,length(shorty)),shorty,'bo')
hold on 
plot(2*ones(1,length(longy)),longy,'ro')
title('Mean y velocities')
xlim([0 3])
ylabel('$$\langle v_{y} \rangle$$ ($\mu$m/min)')
xlabel('Width to length colony proportion')
xticks([1 2])
xticklabels({'1/3','1/5'})
%%dibuixar-les per separat
meanvytlong=mean(vytlong);
meanvytshort=mean(vytshort);
meanvxtlong=mean(vxtlong);
meanvxtshort=mean(vxtshort);

%% std
stdxlong=std(vxtlong);
stdylong=std(vytlong);
stdxshort=std(vxtshort);
stdyshort=std(vytshort);


vt1=0:5:5*(length(vxtlong)-1); vt2=5:5:5*(length(vxtlong)); v=2.5:5:2.5+5*(length(vxtlong)-1);
figure(8)
for i=1:24
    t1=vt1(i); t2=vt2(i); vt=[t1 t2]; 
    plot(vt,meanvytlong(i)*ones(1,2),'-b')
    hold on   
end
errorbar(v,meanvytlong,stdylong,'-o')
xlabel('t (min)')
ylabel('$$ \langle vy \rangle $$ ($\mu$m/min)')
title('Mean $v_{y}$ in function of time for 1/5 colonies')
hold off

figure(9)
for i=1:24
    t1=vt1(i); t2=vt2(i); vt=[t1 t2]; 
    plot(vt,meanvxtlong(i)*ones(1,2),'-b')
    hold on   
end
errorbar(v,meanvxtlong,stdxlong,'-o')
xlabel('t (min)')
ylabel('$$ \langle vx \rangle $$ ($\mu$m/min)')
title('Mean $v_{x}$ in function of time for 1/5 colonies')
hold off

figure(10)
for i=1:24
    t1=vt1(i); t2=vt2(i); vt=[t1 t2]; 
    plot(vt,meanvytshort(i)*ones(1,2),'-b')
    hold on   
end
errorbar(v,meanvytshort,stdyshort,'-o')
xlabel('t (min)')
ylabel('$$ \langle vy \rangle $$ ($\mu$m/min)')
title('Mean $v_{y}$ in function of time for 1/3 colonies')
hold off

figure(11)
for i=1:24
    t1=vt1(i); t2=vt2(i); vt=[t1 t2]; 
    plot(vt,meanvxtshort(i)*ones(1,2),'-b')
    hold on   
end
errorbar(v,meanvxtshort,stdxshort,'-o')
xlabel('t (min)')
ylabel('$$ \langle vx \rangle $$ ($\mu$m/min)')
title('Mean $v_{x}$ in function of time for 1/3 colonies')
hold off

%% with the shaded part

%{
y = rand(1,10); % your mean vector;
x = 1:numel(y);
std_dev = 1;
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'g');
hold on;
plot(x, y, 'r', 'LineWidth', 2);
%}

curve1=meanvytlong+stdylong;
curve2=meanvytlong-stdylong;
x2 = [v, fliplr(v)];
figure(55)
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'c');
hold on
plot(v,meanvytlong,'-bo')
xlabel('t (min)')
ylabel('$$ \langle v_{y} \rangle $$ ($\mu$m/min)')
title('Mean $v_{y}$ in function of time for 1/5 colonies')
hold off

curve1=meanvxtlong+stdxlong;
curve2=meanvxtlong-stdxlong;
x2 = [v, fliplr(v)];
figure(65)
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'c');
hold on
plot(v,meanvxtlong,'-bo')
xlabel('t (min)')
ylabel('$$ \langle v_{x} \rangle $$ ($\mu$m/min)')
title('Mean $v_{x}$ in function of time for 1/5 colonies')
hold off

curve1=meanvxtshort+stdxshort;
curve2=meanvxtshort-stdxshort;
x2 = [v, fliplr(v)];
figure(45)
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'c');
hold on
plot(v,meanvxtshort,'-bo')
xlabel('t (min)')
ylabel('$$ \langle v_{x} \rangle $$ ($\mu$m/min)')
title('Mean $v_{x}$ in function of time for 1/3 colonies')
hold off

curve1=meanvytshort+stdyshort;
curve2=meanvytshort-stdyshort;
x2 = [v, fliplr(v)];
figure(455)
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'c');
hold on
plot(v,meanvytshort,'-bo')
xlabel('t (min)')
ylabel('$$ \langle v_{y} \rangle $$ ($\mu$m/min)')
title('Mean $v_{y}$ in function of time for 1/3 colonies')
hold off

%%

%full vxtshort
a=size(vxtshort);
figure(12)
for i=1:1:a(1)
    plot(v,vxtshort(i,:),'-o')
    hold on
end
title('vx short')


a=size(vytshort);
figure(13)
for i=1:1:a(1)
    plot(v,vytshort(i,:),'-o')
    hold on
end
title('vy short')


a=size(vxtlong);
figure(14)
for i=1:1:a(1)
    plot(v,vxtlong(i,:),'-o')
    hold on
end
title('vx long')

a=size(vytlong);
figure(15)
for i=1:1:a(1)
    plot(v,vytlong(i,:),'-o')
    hold on
end
title('vy long')

%% Percentatge de length escurçat, esto cambiar al meter las otras

%que es esto esto
    dx=x1(:,1)-x1(:,2); perx=dx./x1(:,1);
    dy=y1(:,1)-y1(:,2); pery=dy./y1(:,1);
    perxshort=perx(1:4); %falta 2
    perxlong=perx(5:end);
    peryshort=pery(1:4); %falta 2
    perylong=pery(5:end);

figure(16)
g1 = repmat({'1/3'},1,4);
g2 = repmat({'1/5'},1,length(dx)-4);
g = [g1 g2];
boxplot(perx,g)
title('Length reduction in x direction')
xlabel('Width to length colony proportion')
ylabel('%')
 
figure(17)
g1 = repmat({'1/3'},1,4);
g2 = repmat({'1/5'},1,length(dx)-4);
g = [g1 g2];
boxplot(pery,g)
title('Length reduction in y direction')
xlabel('Width to length colony proportion')
ylabel('%')

%% v=dx/dt

vshortx=dx(1:4)/120;
vlongx=dx(5:end)/120;
vshorty=dy(1:4)/120;
vlongy=dy(5:end)/120;

figure(777)
x = [vshortx' vlongx'];
g1 = repmat({'1/3'},1,length(vshortx));
g2 = repmat({'1/5'},1,length(vlongx));
g = [g1 g2];
boxplot(x,g)
title('Mean x velocity')
ylabel('$$\langle v_{x} \rangle$$ ($\mu$m/min)')
xlabel('Width to length colony proportion')
ylim([0 4]);

figure(888)
x = [vshorty' vlongy'];
g1 = repmat({'1/3'},1,length(vshorty));
g2 = repmat({'1/5'},1,length(vlongy));
g = [g1 g2];
boxplot(x,g)
title('Mean y velocity')
ylabel('$$\langle v_{y} \rangle$$  ($\mu$m/min)')
xlabel('Width to length colony proportion')
ylim([0 4]);
    %% thickness to length 
load("height.mat");
%falta dos 
vh=[vh(1:2), vh(5:end)];
sratio=x1(:,1)'./vh

figure(1000)
g1 = repmat({'1/3'},1,length(vshorty));
g2 = repmat({'1/5'},1,length(vlongy));
g = [g1 g2];
boxplot(sratio,g)
title('Slenderness')
ylabel('Length to thickness ratio')
xlabel('Width to length colony proportion')
