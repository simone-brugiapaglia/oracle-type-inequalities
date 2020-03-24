function pts=Draw_iid_Points_i(p,n)
%% Draw n iid points according to the probability p (which can be d-dimensional)
%
% Author : Pierre Weiss

dimension=length(size(p));
if (dimension==2)
    if (size(p,1)==1 || size(p,2)==1)
        dimension=1;
    end
end
distrib=cumsum(p(:));

pp=rand(n,1);
pp=sort(pp(:));


pts_i=zeros(n,1);
i=1;
iprec=1;
ind=1;
while (i<=n)
    while (distrib(ind)<pp(i))
        ind=ind+1;
    end
    while (i<=n && pp(i)<=distrib(ind))
        i=i+1;
    end
    pts_i(iprec:i-1)=ind;

    iprec=i;
    ind=ind+1;
end

pts=zeros(dimension,n);
if dimension==1
    pts=pts_i;
elseif dimension==2
    [p1,p2]=ind2sub(size(p),pts_i);
    pts(1,:)=p1;pts(2,:)=p2;
elseif dimension==3
    [p1,p2,p3]=ind2sub(size(p),pts_i);
    pts(1,:)=p1;pts(2,:)=p2;pts(3,:)=p3;
end
