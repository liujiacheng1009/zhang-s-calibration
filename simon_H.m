function f = simon_H(H, m, M)
    H=reshape(H,3,3);
    h=H;
    m=[m([1:2],:); ones(1,size(m,2))];
    M=[M([1:2],:); ones(1,size(M,2))];
    
    X=h*M;
    X=[X(1,:)./X(3,:) ; X(2,:)./X(3,:); X(3,:)./X(3,:)];
    
    res=m-X;
    req=[res(1,:), res(2,:)]; 
    f = req;
end

