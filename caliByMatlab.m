imgcount = 5;%需要处理的图片数量
squareSize = 0.031;
imageFileNames = {1,imgcount};
for i = 1:imgcount
  imageFileNames{i} = sprintf('images/%d.jpg', i);
end
[imagePoints,boardSize,imagesUsed] = detectCheckerboardPoints(imageFileNames);
m=permute(imagePoints,[2 1 3]);
worldPoints= generateCheckerboardPoints(boardSize,squareSize);
I = imread(imageFileNames{1});
imageSize = [size(I,1),size(I,2)];
cameraParams = estimateCameraParameters(imagePoints,worldPoints, ...
                                       'ImageSize',imageSize);
cameraParams.RadialDistortion
(cameraParams.IntrinsicMatrix)'

