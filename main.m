clc
tic
I1= imread('C:\Users\Administrator\Desktop\demo2\data\gang1.bmp');  
I1=rgb2gray(I1);  %把RGB图像变成灰度图像
I2= imread('C:\Users\Administrator\Desktop\demo2\data\gang2.bmp');   
I2=rgb2gray(I2);

%寻找特征点
points1 = detectKAZEFeatures(I1);  %读取特征点,'MinQuality',0.005
points2 = detectKAZEFeatures(I2);
% points1 = detectMSERFeatures(I1); 
% points2 = detectMSERFeatures(I2); 
% points1 = detectBRISKFeatures(I1); 
% points2 = detectBRISKFeatures(I2);  
 
%Extract the features.计算描述向量  
[f1, vpts1] = extractHOGFeatures(I1, points1);  
[f2, vpts2] = extractHOGFeatures(I2, points2);  
 
%进行匹配  
indexPairs = matchFeatures(f1, f2,'Method','NearestNeighborSymmetric','MatchThreshold',100) ;
%用'Method','NearestNeighborSymmetric'，阈值调节用'MatchThreshold',范围0-100，表示选择最强的匹配的百分比，越大匹配点越多
%用'Method','NearestNeighborRatio'，阈值调节'MaxRatio',0.7，范围0-1，默认0.6，较大该值获得较多匹配
matched_pts1 = vpts1(indexPairs(:, 1));
matched_pts2 = vpts2(indexPairs(:, 2));
 
resultpairs1 = matched_pts1.Location; %存储匹配点坐标
resultpairs2 = matched_pts2.Location;
 
%显示匹配
figure;
showMatchedFeatures(I1,I2,matched_pts1,matched_pts2,'montage');  
legend('matched points 1','matched points 2'); 

% % ransac去除误匹配点
%  Ipts1=(resultpairs1)';
%  Ipts2=(resultpairs2)';
%  [H, inliers] = ransacfithomography(Ipts1,Ipts2, .05);
%  inlierBoxPoints = [Ipts1(2,inliers)' Ipts1(1,inliers)'];
%  inlierScenePoints = [Ipts2(2,inliers)' Ipts2(1,inliers)'];
% figure;
% showMatchedFeatures(I1, I2, inlierBoxPoints, inlierScenePoints, 'montage','Parent',axes);
% title('RANSAC去除误匹配点');
% %MSAC去除误匹配点
[tform, inBoxPoints, inScenePoints] = estimateGeometricTransform(matched_pts1, matched_pts2, 'projective');
figure;
showMatchedFeatures(I1, I2, inBoxPoints, inScenePoints,'montage','PlotOptions',{'ro','y+','g'});
toc