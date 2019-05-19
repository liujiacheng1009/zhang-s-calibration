%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%估计单应性矩阵H,参考原文附录A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = homography2d(x1,x2)
    M=x1;                             
    m=x2;
    %数据规范化
    [x1, T1] = normalise2dpts(x1);
    [x2, T2] = normalise2dpts(x2);
    Npts = length(x1);
    A = zeros(3*Npts,9);   
    O = [0 0 0];
    for n = 1:Npts
        X = x1(:,n)';
        x = x2(1,n);y = x2(2,n); w = x2(3,n);
        A(3*n-2,:) = [  O  -w*X  y*X];
        A(3*n-1,:) = [ w*X   O  -x*X];
        A(3*n  ,:) = [-y*X  x*X   O ];
    end
    [~,~,V] = svd(A);
    H1 = reshape(V(:,9),3,3)';
    H2= T2\H1*T1;
    H=H2/H2(3,3);
    options = optimset('Algorithm', 'levenberg-marquardt');
    [x,~,~,~,~]  = lsqnonlin( @simon_H, reshape(H,1,9) , [],[],options,m, M);
    H=reshape(x,3,3);
    H=H/H(3,3);
end
    
    

