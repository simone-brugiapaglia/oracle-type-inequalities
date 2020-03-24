%Computes the 1D Legendre matrix of size length(pts) x k, where t is the set of
%data points

function A=LegMat(pts,k)
pts = pts(:);
A=zeros(length(pts),k);
A(:,1)=1;
A(:,2)=pts*sqrt(3);
for i=2:k-1
    A(:,i+1)=(pts.* (2*i - 1).*A(:,i)./sqrt(i-1/2)- (i - 1).*A(:,i-1)./sqrt(i-3/2)).*sqrt(i+1/2)/i;
end
end

