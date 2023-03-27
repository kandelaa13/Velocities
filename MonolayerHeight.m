vh=[];
for i=1:7
     for j=0:1
path = ['C:\Users\kballerini\Desktop\Lab\11112022\20221111_noY27\RoiSet_f000', num2str(i), '_', num2str(j), '.zip'];
%path='C:\Users\kballerini\Desktop\Lab\11112022\20221111_noY27\RoiSet_f0002_0.zip';
sROI = ReadImageJROI(path);


%ROI_i = sROI{1,i}; %Still need to loop for all i's, if needed.
ROI_1 = sROI{1,1}.mnCoordinates;
figure(j+1)
plot(ROI_1(:,1),ROI_1(:,2)); axis equal; hold on

center = mean(ROI_1);

% Spot points over/below center of ROI

is_above = ROI_1(:,2)>center(2); %logical
is_below = ~is_above;

% Spot leftmost and rightmost points and the 20um around them

max_right = max(ROI_1(:,1));
is_right = ROI_1(:,1)>(max_right-20);
min_left = min(ROI_1(:,1));
is_left = ROI_1(:,1)<(min_left+20);

% Now it's just time to make the appropiate mean
% ~ makes the complementary 0 --> 1

above_line = logical(is_above.*(~is_right).*(~is_left));
below_line = logical(is_below.*(~is_right).*(~is_left));
ROI_1_y = ROI_1(:,2);

top = mean(ROI_1_y(above_line)); 
bottom = mean(ROI_1_y(below_line));
h=abs(top-bottom) %this result is in um
vh=[vh h];

x = min_left:max_right;
plot(x,ones(1,numel(x)).*top,'r');
plot(x,ones(1,numel(x)).*bottom,'r');
axis([min_left-30 max_right+30 0 200])
hold off
     end 
end
save('height.mat',"vh");