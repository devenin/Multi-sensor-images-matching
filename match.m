%%%file2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ncc算法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res=match(a1,cnt1,r1,c1,a2,cnt2,r2,c2)
% res=match(a1,a2)
% 将从a1寻找a2中的最佳匹配点，得到从a2中抽取的res，也就是单向搜索
% [result1,cnt1,r11,c11]=harris(a1);
% [result2,cnt2,r22,c22]=harris(a2);%可以保证想匹配哪些点就匹配哪些
% figure;
%  imshow(result1);title('result1角点位置');
%  figure;title('result2角点位置');
%  imshow(result2);   
 
% win=[1/9 1/9 1/9;1/9 1/9 1/9;1/9 1/9 1/9];
% u1=filter2(win,a1);
% u2=filter2(win,a2);  %求均值
%                           
%       
% 
a1=double(a1);
a2=double(a2);
% A=filter2(win,(a1-u1).^2);%求方差
% B=filter2(win,(a2-u2).^2);
[m1,n1]=size(a1);
[m2,n2]=size(a2);
res1=zeros(m1,n1);
res2=zeros(m2,n2);           %寻找的匹配的点

for s=1:cnt1           
        max=0; p=0;q=0;i=r1(s,1);j=c1(s,1);      %p.q存放坐标
        for v=1:cnt2
            m=r2(v,1);n=c2(v,1);
            u1(i,j)=(a1(i-1,j-1)+a1(i-1,j)+a1(i-1,j+1)+a1(i,j-1)+a1(i,j)+a1(i,j+1)+a1(i+1,j-1)+a1(i+1,j)+a1(i+1,j+1))/9;%%求均值
            u2(m,n)=(a2(m-1,n-1)+a2(m-1,n)+a2(m-1,n+1)+a2(m,n-1)+a2(m,n)+a2(m,n+1)+a2(m+1,n-1)+a2(m+1,n)+a2(m+1,n+1))/9;
            
      x1=(a1(i-1,j-1)-u1(i,j));y1=(a2(m-1,n-1)-u2(m,n));
      x2=(a1(i-1,j)-u1(i,j));y2=(a2(m-1,n)-u2(m,n));
      x3=(a1(i-1,j+1)-u1(i,j));y3=(a2(m-1,n+1)-u2(m,n));
      x4=(a1(i,j-1)-u1(i,j));y4=(a2(m,n-1)-u2(m,n));
      x5=(a1(i,j)-u1(i,j));y5=(a2(m,n)-u2(m,n));
      x6=(a1(i,j+1)-u1(i,j));y6=(a2(m,n+1)-u2(m,n));
      x7=(a1(i+1,j-1)-u1(i,j));y7=(a2(m+1,n-1)-u2(m,n));
      x8=(a1(i+1,j)-u1(i,j));y8=(a2(m+1,n)-u2(m,n));
      x9=(a1(i+1,j+1)-u1(i,j));y9=(a2(m+1,n+1)-u2(m,n));
      num1=x1^2+x2^2+x3^2+x4^2+x5^2+x6^2+x7^2+x8^2+x9^2;
      num2=y1^2+y2^2+y3^2+y4^2+y5^2+y6^2+y7^2+y8^2+y9^2;
      
      
      k1=x1*y1;  %用result是找不到什么信息，3*3对于如此稀疏点没有用
      k2=x2*y2;
      k3=x3*y3;  %用循环要好一些
      k4=x4*y4;
      k5=x5*y5;
      k6=x6*y6;
      k7=x7*y7;
      k8=x8*y8;
      k9=x9*y9;
      num=k1+k2+k3+k4+k5+k6+k7+k8+k9;
       den=sqrt(num1*num2);
       ncc=num/den;
        if ncc>max
            max=ncc;
            p=m;
            q=n;
        end
       
        end
 
  res2(p,q)=1;
 
end
                %做不出来可以先搜索一个点
res=res2;
return