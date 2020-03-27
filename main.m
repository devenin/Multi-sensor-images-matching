clc
tic
I1= imread('C:\Users\Administrator\Desktop\demo2\data\gang1.bmp');  
I1=rgb2gray(I1);  %��RGBͼ���ɻҶ�ͼ��
I2= imread('C:\Users\Administrator\Desktop\demo2\data\gang2.bmp');   
I2=rgb2gray(I2);

%Ѱ��������
points1 = detectKAZEFeatures(I1);  %��ȡ������,'MinQuality',0.005
points2 = detectKAZEFeatures(I2);
% points1 = detectMSERFeatures(I1); 
% points2 = detectMSERFeatures(I2); 
% points1 = detectBRISKFeatures(I1); 
% points2 = detectBRISKFeatures(I2);  
 
%Extract the features.������������  
[f1, vpts1] = extractHOGFeatures(I1, points1);  
[f2, vpts2] = extractHOGFeatures(I2, points2);  
 
%����ƥ��  
indexPairs = matchFeatures(f1, f2,'Method','NearestNeighborSymmetric','MatchThreshold',100) ;
%��'Method','NearestNeighborSymmetric'����ֵ������'MatchThreshold',��Χ0-100����ʾѡ����ǿ��ƥ��İٷֱȣ�Խ��ƥ���Խ��
%��'Method','NearestNeighborRatio'����ֵ����'MaxRatio',0.7����Χ0-1��Ĭ��0.6���ϴ��ֵ��ý϶�ƥ��
matched_pts1 = vpts1(indexPairs(:, 1));
matched_pts2 = vpts2(indexPairs(:, 2));
 
resultpairs1 = matched_pts1.Location; %�洢ƥ�������
resultpairs2 = matched_pts2.Location;
 
%��ʾƥ��
figure;
showMatchedFeatures(I1,I2,matched_pts1,matched_pts2,'montage');  
legend('matched points 1','matched points 2'); 

% % ransacȥ����ƥ���
%  Ipts1=(resultpairs1)';
%  Ipts2=(resultpairs2)';
%  [H, inliers] = ransacfithomography(Ipts1,Ipts2, .05);
%  inlierBoxPoints = [Ipts1(2,inliers)' Ipts1(1,inliers)'];
%  inlierScenePoints = [Ipts2(2,inliers)' Ipts2(1,inliers)'];
% figure;
% showMatchedFeatures(I1, I2, inlierBoxPoints, inlierScenePoints, 'montage','Parent',axes);
% title('RANSACȥ����ƥ���');
% %MSACȥ����ƥ���
[tform, inBoxPoints, inScenePoints] = estimateGeometricTransform(matched_pts1, matched_pts2, 'projective');
figure;
showMatchedFeatures(I1, I2, inBoxPoints, inScenePoints,'montage','PlotOptions',{'ro','y+','g'});
toc