%% A script to segment cells from BF images

%%https://es.mathworks.com/matlabcentral/answers/76556-improve-ellipse-fit-around-a-binary-image
%%https://es.mathworks.com/matlabcentral/fileexchange/15125-fitellipse-m
%A simple method for fitting of bounding rectangle to closed regions
%https://es.mathworks.com/matlabcentral/fileexchange/71491-largest-inscribed-rectangle-square-or-circle
%https://www.mathworks.com/matlabcentral/fileexchange/28155-inscribed_rectangle#functions_tab

%https://www.unioviedo.es/compnum/labs/lab05_interpol/lab05_interpol.html


%Considering maximum number of colonies of interest=2
clear all
close all

set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

%rootname = 'C:\Users\kande\Desktop\LAB\11112022\20221111_noY27\20221111_BF_10x\20221111_BF_10x_f0006.tif';
rootname = 'C:\Users\kballerini\Desktop\Lab\11112022\20221111_noY27\20221111_BF_10x\20221111_BF_10x_f0007.tif';
%figure 4 peta
%AUMENTA EL TAMAÑO DE LA MATRIZ, Y  EN LA ULTIMA APARECEN LINEAS --> son
%culpa de la delta
%tuning parameters: delta, pixlim and <10



area_min=10000; 
min_object_size=area_min; %smallest object in segmented image that will be kept
min_hole_size=1000;  %800; %holes larger than this in the segmented image will be filled in
treshold_finetune=0.0001; 
delta=8; %multiple of 2048, so that the programme is faster
pixlim=80; %minimum # of valid subsequent pixels, the bigger the more it cuts
limcolonies=5;
num=1;

Fluo_tiff_info = imfinfo(rootname);
K = length(Fluo_tiff_info); %# of images in the stack


meanxder1=zeros(1,length(K)); meanxizq1=zeros(1,length(K));
meanxder2=zeros(1,length(K)); meanxizq2=zeros(1,length(K));
meanydown2=zeros(1,length(K)); meanyup2=zeros(1,length(K));
meanydown1=zeros(1,length(K)); meanyup1=zeros(1,length(K));
At=1;

for q=1:1:K
  q
%segmentation
imBF1 = double(imread(rootname,q));
res_EGT = EGT_Segmentation(imBF1, min_object_size,min_hole_size,treshold_finetune);

%just stay with the biggest islands
CC0 = bwconncomp(res_EGT);
numPixels = cellfun(@numel,CC0.PixelIdxList); %vec
[biggest,idx] = max(numPixels);
def=zeros(size(res_EGT));
def(CC0.PixelIdxList{idx}) = 1;
numsort=sort(numPixels);
if fix(numsort(end-1)/numsort(end))<limcolonies %of the same order of magnitude
      idx2=numsort(end-1); 
      k=find(numPixels==idx2);
      def(CC0.PixelIdxList{k}) = 1;
end
def=logical(def); %not necessary
CC1=bwconncomp(def);


%Start sectioning row-wise
reduceddef=def(1:delta:end,:); %reduced matrix
divisions=length(def(1,:))/delta;
verase=[]; j=0; 
vec=find(reduceddef'); %finds values dif from 0
tic
step=vec(2)-vec(1); i=2;

while i<length(vec)
    verase=[verase vec(i-1)]; 
    if step==1
        j=j+1;
        step=vec(i+1)-vec(i);
        i=i+1;
        if j==pixlim 
            vv=vec(i+1:end-1)-vec(i:end-2);
            i=i+find(vv>1,1); 
            j=0; verase=[]; step=1;
            if isempty(i)
                i=length(vec)+1;
            end
        end
    else
         col=fix(verase/(divisions*delta))+1;  
         col=col*delta-(delta-1);
         row=mod(verase,(divisions*delta));
         n=find(row==0); 
         if isempty(n)
         else  
             row(n)=divisions*delta;
         end
         def(col,row)=0; 
         step=1; j=0; i=i+1; verase=[];
    end
    if i==length(vec)
        verase=[verase vec(end-1) vec(end)];  
        verase=verase(2:end);
        col=fix(verase/(divisions*delta))+1;   
        col=col*delta-(delta-1);
        row=mod(verase,(divisions*delta));  
        n=find(row==0); 
        if isempty(n)
        else  
           row(n)=divisions*delta;
        end
        def(col,row)=0;
    end
end

CC2 = bwconncomp(def);
numPixels = cellfun(@numel,CC2.PixelIdxList); 
[biggest,idx] = max(numPixels);
def1=zeros(size(res_EGT));
def1(CC2.PixelIdxList{idx}) = 1;
numsort=sort(numPixels);
if fix(numsort(end-1)/numsort(end))<limcolonies %of the same order of magnitude
      idx2=numsort(end-1); 
      k=find(numPixels==idx2);
      def1(CC2.PixelIdxList{k}) = 1;
end
def1=logical(def1);
CC3 = bwconncomp(def1);

%Sectioning column-wise
reduceddeff=def1(:,1:delta:end);
divisions=length(def1(1,:))/delta;
verase=[]; j=0; 
vec=find(reduceddeff); 
tic
step=vec(2)-vec(1); i=2;
while i<length(vec)
    verase=[verase vec(i-1)]; 
    if step==1
        j=j+1;
        step=vec(i+1)-vec(i);
        i=i+1;
        if j==pixlim 
            vv=vec(i+1:end-1)-vec(i:end-2);
            i=i+find(vv>1,1); 
            j=0; verase=[]; step=1;
            if isempty(i)
                i=length(vec)+1;
            end
        end
    else
         row=fix(verase/(divisions*delta))+1;  
         row=row*delta-(delta-1);
         col=mod(verase,(divisions*delta));
         n=find(col==0); 
         if isempty(n)
         else  
             col(n)=divisions*delta;
         end
         def1(col,row)=0; 
         step=1; j=0; i=i+1; verase=[];
    end
    if i==length(vec)
        verase=[verase vec(end-1) vec(end)];  
        verase=verase(2:end);
        row=fix(verase/(divisions*delta))+1;   
        row=row*delta-(delta-1);
        col=mod(verase,(divisions*delta));  
        n=find(col==0); 
        if isempty(n)
        else  
           col(n)=divisions*delta;
        end
        def1(col,row)=0;
    end
end

CC4 = bwconncomp(def1);
numPixels = cellfun(@numel,CC4.PixelIdxList); 
[biggest,idx] = max(numPixels);
def2=zeros(size(res_EGT));
def2(CC4.PixelIdxList{idx}) = 1;
numsort=sort(numPixels);
if fix(numsort(end)/numsort(end-1))<limcolonies %of the same order of magnitude
      idx2=numsort(end-1); 
      k=find(numPixels==idx2);
      def2(CC4.PixelIdxList{k}) = 1;
end
def2=logical(def2);
CC5 = bwconncomp(def2);
im1=figure(1);
imshow(imBF1,[]); axis equal; colormap gray;hold on; %this just the image
visboundaries(def2) 
drawnow

pause(2)
title(['EGT'])
imBF1=uint16(imBF1);

% AFTER SECTIONING COLONY SEPARATION



numPixels = cellfun(@numel,CC5.PixelIdxList);
a=CC5.PixelIdxList{1}; amin=min(a); amax=max(a);
%(rowmina,colmina), da las coordenadas del pixel "mas pequeño", el que
%tiene la columna mas txiki, lo usamos solo para esto
rowmina=mod(amin,2048); colmina=fix(amin/2048)+1; 
rowmaxa=mod(amax,2048);
b=CC5.PixelIdxList{2}; bmin=min(b);  bmax=max(b);
rowminb=mod(bmin,2048); rowmaxb=mod(bmax,2048);
if rowmina<rowminb
    up=a;
    down=b;
else
    up=b;
    down=a;
end

% COLONY 1, BORDER
one=zeros(size(res_EGT));
one(up) = 1;
mask1=bwperim(one); maskcoor1=find(mask1);

col1=fix(maskcoor1/2048)+1;  
 row1=mod(maskcoor1,2048); 
 n1=find(row1==0); 
 if isempty(n1)
 else  
     row1(n1)=2048;
 end
rowmin1=min(row1); rowmax1=max(row1);
colmin1=min(col1); colmax1=max(col1);

%con esto, golpeando desde la izquierda, encontramos el pixel de la columna
%mas pequeña para cada row, y el mas grande tambien
out=round((rowmax1-rowmin1)*0.05);
vec1=rowmin1+out:1:rowmax1-out; border1=zeros(1,length(vec1)); bordeizq1=zeros(1,length(vec1));
for jj=1:1:length(vec1)
     cizq1=find(one(vec1(jj),:),1);
     bordeizq1(jj)=cizq1;
     cder1=find(one(vec1(jj),:),1,'last');
     border1(jj)=cder1;
end


%con esto, golpeando desde ARRIBA, encontramos el pixel de la row
%mas pequeña para cada columna, y el mas grande tambien
out=round((colmax1-colmin1)*0.05);
 vec11=colmin1+out:1:colmax1-out; bordup1=zeros(1,length(vec11)); bordown1=zeros(1,length(vec11));
for jjj=1:1:length(vec11)
     cup1=find(one(:,vec11(jjj)),1);
     bordup1(jjj)=cup1;
     cdown1=find(one(:,vec11(jjj)),1,'last');
     bordown1(jjj)=cdown1;
end

 %se deberia coger un valor el mas similar al anterior y dejarlo como ref//
 %o mirar la media de este (probamos el primero)
 
  %coger un valor de en medio, y quedarse con ese como referencia, calcular
  %distancia relativa y eliminar los que no
 if q==1
      meanizq1=mean(bordeizq1);
       meander1=mean(border1); 
        meanup1=mean(bordup1);
        meandown1=mean(bordown1);
        limh=40; limv=70;


%      if mod(jj,2)==0
%       meanizq1=bordeizq1(jj/2);
%        meander1=border1(jj/2);
%   
%      else
%           meanizq1=bordeizq1((jj+1)/2);
%          meander1=border1((jj+1)/2);
%      end
%      if mod(jjj,2)==0
%             meanup1=bordup1(jjj/2);
%        meandown1=bordown1(jjj/2);
%      else
%          meanup1=bordup1((jjj+1)/2);
%          meandown1=bordown1((jjj+1)/2);
%      end
 else
     meanizq1=meanxizq1(q-1); %referencia la media del paso anterior
     meander1=meanxder1(q-1);
     meanup1=meanyup1(q-1);
     meandown1=meanydown1(q-1);
     limh=50; limv=50;
 end

 %acceptaremos aprox 100 pixeles de diferencia, nos tenemos que quedar con
 %los indices
 vecder1=vec1; vecizq1=vec1; vecup1=vec11; vecdown1=vec11;
 iborder1=find(abs(meander1-border1)>limh);
 ibordeizq1=find(abs(meanizq1-bordeizq1)>limh);
 ibordup1=find(abs(meanup1-bordup1)>limv);
ibordown1=find(abs(meandown1-bordown1)>limv);


border1(iborder1) = []; vecder1(iborder1)=[];
bordeizq1(ibordeizq1) = []; vecizq1(ibordeizq1)=[];
bordown1(ibordown1) = []; vecdown1(ibordown1)=[];
bordup1(ibordup1) = []; vecup1(ibordup1)=[];



%quizas mirar que haya un minimo de puntos, o un minimo de puntos entre uno
%y otro? sino mirar que no haya pixeles que crezcan

%primero calcular la media 
 meanxder1(q)=mean(border1);
meanxizq1(q)=mean(bordeizq1);
meanyup1(q)=mean(bordup1);
meanydown1(q)=mean(bordown1);

% COLONY 2
two=zeros(size(res_EGT));
two(down) = 1;
mask2=bwperim(two); maskcoor2=find(mask2);

col2=fix(maskcoor2/2048)+1;  
 row2=mod(maskcoor2,2048); 
 n2=find(row2==0); 
 if isempty(n2)
 else  
     row2(n2)=2048;
 end
rowmin2=min(row2); rowmax2=max(row2);
colmin2=min(col2); colmax2=max(col2);

%con esto, golpeando desde la izquierda, encontramos el pixel de la columna
%mas pequeña para cada row, y el mas grande tambien
out=round((rowmax2-rowmin2)*0.1); 
vec2=rowmin2+out:1:rowmax2-out; border2=zeros(1,length(vec2)); bordeizq2=zeros(1,length(vec2));
for jj=1:1:length(vec2)
     cizq2=find(two(vec2(jj),:),1);
     bordeizq2(jj)=cizq2;
     cder2=find(two(vec2(jj),:),1,'last');
     border2(jj)=cder2;
end

%con esto, golpeando desde ARRIBA, encontramos el pixel de la row
%mas pequeña para cada columna, y el mas grande tambien
out=round((colmax2-colmin2)*0.05);
vec22=colmin2+out:1:colmax2-out; bordup2=zeros(1,length(vec22)); bordown2=zeros(1,length(vec22));
for jjj=1:1:length(vec22)
     cup2=find(two(:,vec22(jjj)),1);
     bordup2(jjj)=cup2;
     cdown2=find(two(:,vec22(jjj)),1,'last');
     bordown2(jjj)=cdown2;
end



 %se deberia coger un valor el mas similar al anterior y dejarlo como ref//
 %o mirar la media de este (probamos el primero)
 
  %coger un valor de en medio, y quedarse con ese como referencia, calcular
  %distancia relativa y eliminar los que no
 if q==1
     meanizq2=mean(bordeizq2);
       meander2=mean(border2); 
        meanup2=mean(bordup2);
        meandown2=mean(bordown2);
        limh=87; limv=70;
%      if mod(jj,2)==0
%       meanizq2=bordeizq2(jj/2);
%        meander2=border2(jj/2);
%   
%      else
%           meanizq2=bordeizq2((jj+1)/2);
%          meander2=border2((jj+1)/2);
%      end
%      if mod(jjj,2)==0
%             meanup2=bordup2(jjj/2);
%        meandown2=bordown2(jjj/2);
%      else
%          meanup2=bordup2((jjj+1)/2);
%          meandown2=bordown2((jjj+1)/2);
%      end
 else
     meanizq2=meanxizq2(q-1); %referencia la media del paso anterior
     meander2=meanxder2(q-1);
     meanup2=meanyup2(q-1);
     meandown2=meanydown2(q-1);
     limh=40; limv=80;
 end
 %acceptaremos aprox 100 pixeles de diferencia, nos tenemos que quedar con
 %los indices
 vecder2=vec2; vecizq2=vec2; vecup2=vec22; vecdown2=vec22;
 iborder2=find(abs(meander2-border2)>limh);
 ibordeizq2=find(abs(meanizq2-bordeizq2)>limh);
 ibordup2=find(abs(meanup2-bordup2)>limv);
ibordown2=find(abs(meandown2-bordown2)>limv);


border2(iborder2) = []; vecder2(iborder2)=[];
bordeizq2(ibordeizq2) = []; vecizq2(ibordeizq2)=[];
bordown2(ibordown2) = []; vecdown2(ibordown2)=[];
bordup2(ibordup2) = []; vecup2(ibordup2)=[];


%quizas mirar que haya un minimo de puntos, o un minimo de puntos entre uno
%y otro? sino mirar que no haya pixeles que crezcan

%primero calcular la media 
 meanxder2(q)=mean(border2);
meanxizq2(q)=mean(bordeizq2);
meanyup2(q)=mean(bordup2);
meanydown2(q)=mean(bordown2);




  plot(border1,vecder1,'o')
% % figure(1)
   hold on
    plot(bordeizq1,vecizq1,'o','DisplayName','izq1')
   hold on
    plot(border2,vecder2,'o')
    hold on
    plot(bordeizq2,vecizq2,'o')
   hold on
   plot(vecup1,bordup1,'o','DisplayName','up1')
   hold on 
   plot(vecup2,bordup2,'o','DisplayName','up2')
   hold on 
   plot(vecdown1,bordown1,'o','DisplayName','down1')
   hold on 
    plot(vecdown2,bordown2,'o','DisplayName','down2')
   % legend
   hold on 
   f=gcf;
   pause(10)
   exportgraphics(f,'f0007.gif',"Append",true)
  % exportgraphics(f,'AnnotatedPlot.tif')
  hold off

% [rectx,recty,area,perimeter] = minboundrect(col1,row1);
% plot(rectx,recty,'b','LineWidth',2)

% S=FindLargestRectangles(one);
% %[C, H, W]=FindLargestRectangles(two);
% [~,pos]=max(S(:));
% [r c]=ind2sub(size(S), pos);
% rectangle('Position',[c,r,S(r,c),S(r,c)], 'EdgeColor','r', 'LineWidth',3);
% 
% %this way non-tilted
% LRout=LargestRectangle(one,0,0,0,0,0);
% xco=LRout(2:end,1); xco=[xco; LRout(2,1)]; yco=LRout(2:end,2); yco=[yco; LRout(2,2)];
% %plot(xco,yco,'b','LineWidth',2.5)
% 
%  meanxder1(q)=LRout(3,1);
% meanxizq1(q)=LRout(2,1);


end %main loop

matpos=[meanxizq1; meanyup1; meanxder1; meanydown1; meanxizq2; meanyup2; meanxder2; meanydown2];
save('f0007.mat',"matpos");
%% DADES, faltan en y// colonia 2
%lets pass it to micrometres
% meanxizq=meanxizq*0.65;
% meanxizq1=meanxizq1*0.65; es un 10x
% close all
% meanxderr1=flip(meanxder1);%*0.01625;
% 
% figure(3)
% plot(1:1:K,meanxizq1,'-ro','DisplayName','method1')
% legend
% title('pos izq')
% 
% 
% figure(4)
% plot(1:1:K,meanxderr1,'-ro','DisplayName','method1')
% legend
% title('pos der')
% 
% %supose time units of one between images
% velder1=(meanxderr1(2:end)-meanxderr1(1:end-1))/At; meanvder1=mean(velder1);
% velizq1=(meanxizq1(2:end)-meanxizq1(1:end-1))/At; meanvizq1=mean(velizq1);
% 
% 
% figure(5)
% plot(1:1:K-1,velder1,'-bo','DisplayName','der')
% hold on 
% plot(1:1:K-1,meanvder1*ones(1,K-1),'c')
% hold on
% plot(1:1:K-1,velizq1,'-ro','DisplayName','izq')
% hold on
% plot(1:1:K-1,meanvizq1*ones(1,K-1),'m')
% hold on
% plot(1:1:K-1,zeros(1,K-1),'g')
% legend
% title('v')

%matpos=[meanxizq1; meanyup1; meanxder1; meanydown1; meanxizq2; meanyup2; meanxder2; meanydown2];
%save('f0002.mat',"matpos");

 
 %% not successful trials
 
 %{
%no esta ben fet pero maso
 plot(rowb,colb,'o')

 %CH_objects = bwconvhull(def2,'objects');
%imshow(CH_objects);
xk=[rowb, colb];
[k,av] = convhull(xk);
plot(xk(:,1),xk(:,2))
hold on
plot(xk(k,1),xk(k,2))
 
xx=[rowb' ; colb'];
[z, a, b, alpha] = fitellipse(xx, 'linear', 'constraint', 'trace');
 plotellipse(z, a, b, alpha)



%THIS FIT IS NOT WELL DONE
ellipse_tb = fit_ellipse(rowb,colb);
%%
                 % Parameterize the equation.
long_axisb= ellipse_tb.long_axis;
short_axisb= ellipse_tb.short_axis;
t = linspace(0, 360,720);
phaseShiftb = ellipse_tb.phi;
xAmplitudeb = cos(phaseShiftb)*long_axisb; %1269.27894590088;
yAmplitudeb = sin(phaseShiftb)*short_axisb; %723.919768548135;
xb = xAmplitudeb * sind(t + phaseShiftb);
yb = yAmplitudeb * cosd(t);


%draw ellipse
xcb=ellipse_tb.X0_in; ycb=ellipse_tb.Y0_in;
a=ellipse_tb.a; %xRad
b=ellipse_tb.b; %yRad
m = 1000;
x = zeros(m,1);
y = zeros(m,1);
theta = linspace(0,2*pi,m);
for k = 1:m
        x(k) = a * cos(theta(k));
        y(k) = b * sin(theta(k));
end
alpha = phaseShiftb;
R  = [cos(alpha) -sin(alpha); sin(alpha)  cos(alpha)];
rCoords = R*[x' ; y'];   
xr = rCoords(1,:)';      
yr = rCoords(2,:)';      
%plot(x+xcb,y+ycb,'r');
grid on;
hold on;
plot(xr+xcb,yr+ycb,'b');



%%


%{
% Now plot the rotated ellipse.
plot(xb, yb, 'b-', 'LineWidth', 2);
axis equal
grid on;
xlabel('X', 'FontSize', fontSize);
ylabel('Y', 'FontSize', fontSize);
%title('Rotated Ellipses', 'FontSize', fontSize);
hold off
%}
 %}
 %}
         %% Tendre que separar los cuadrados para fittearlos
%{
         [rectx,recty,area,perimeter] = minboundrect(row2,col2,'a');
% rectfit=regionprops(l, 'BoundingBox'); fit1=rectfit.BoundingBox(1,1); fit2=rectfit.BoundingBox(2,1);
ellipse_t = fit_ellipse( row2,col2);

         % Parameterize the equation.
long_axis= ellipse_t.long_axis;
short_axis= ellipse_t.short_axis;
t = linspace(0, 360,500);
phaseShift = ellipse_t.phi;
xAmplitude = cos(phaseShift)*long_axis; %1269.27894590088;
yAmplitude = sin(phaseShift)*short_axis; %723.919768548135;
x = xAmplitude * sind(t + phaseShift);
y = yAmplitude * cosd(t);

% Now plot the rotated ellipse.
plot(x, y, 'b-', 'LineWidth', 2);
axis equal
grid on;
xlabel('X', 'FontSize', fontSize);
ylabel('Y', 'FontSize', fontSize);
%title('Rotated Ellipses', 'FontSize', fontSize);
hold off
%}

%{
if q==1
imwrite(imBF1,newfilename,'png');%print('-dtiff',newfilename) %imwrite(fig,file,"tif");
imwrite(l,newfilename,'png','WriteMode', 'append');
else
%imwrite(a,newfilename,'tif', ']WriteMode', 'append');%print('-dtiff','-append',newfilename);%imwrite(fig,"tif",file, 'WriteMode', 'append');
end

end
%}
%{
t=Tiff('Stacked.tiff','w');
        tagstruct.ImageLength = y1; % image height
        tagstruct.ImageWidth = x1; % image width
        tagstruct.Photometric = Tiff.Photometric.RGB; % https://de.mathworks.com/help/matlab/ref/tiff.html
        tagstruct.BitsPerSample = 8;
        tagstruct.SamplesPerPixel = 3;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; % groups rgb values into a single pixel instead of saving each channel separately for a tiff image
        tagstruct.Software = 'MATLAB';
        setTag(t,tagstruct)
        write(t,squeeze(im2uint8(Imagelayer1)));
        
        writeDirectory(t);
        setTag(t,tagstruct)
        write(t,squeeze(im2uint8(Imagelayer2))) %%%appends the next layer to the same file t
        
        % do this for as many as you need, or put it in a loop if you can
        close(t) %%% this is necessary otherwise you won't be able to open it in imageJ etc to double check, unless you close matlab
%}






