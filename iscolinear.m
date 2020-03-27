function r = iscolinear(p1,p2,p3,flag)

    if nargin ==3
        flag='inhomog';
    end
    if ~all(size(p1)==size(p2))|~all(size(p1)==size(p3))|...
       ~(length(p1)==2|length(p1)==3)
        error('points must have the same dimension of 2 or 3');
    end
    if length(p1)==2
         p1(3)=1;p2(3)=1;p3(3)=1;
    end
    if flag(1)=='h'
        r=abs(dot(cross(p1,p2),p3))<eps;
    else
        r=norm(cross(p2-p1,p3-p1))<eps;
    end
end