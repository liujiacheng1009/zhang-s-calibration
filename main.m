imgcount = 5;%需要处理的图片数量
squareSize = 0.031;
imageFileNames = {1,imgcount};
for i = 1:imgcount
  imageFileNames{i} = sprintf('images/%d.jpg', i);
end
[imagePoints,boardSize,imagesUsed] = detectCheckerboardPoints(imageFileNames);
m=permute(imagePoints,[2 1 3]);
M = generateCheckerboardPoints(boardSize,squareSize)';
%M是棋盘格角点世界坐标系坐标
%m是5张图片中的棋盘格角点图像坐标系坐标
[k1,k2,A]=Zhang(M,m);
k1,k2,A