function correctPoints = match_ransac(P1,P2, it, N, t)
% MATCH_RANSAC:RANSAC outliner detector for matched point-pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%       P1        - coordinates of points in the 1st image, [pointnum x 2]
%       P2        - coordinates of points in the 2nd image, [pointnum x 2]
%       it        - iterations 迭代次数(the number of subset to try)
%       N         - number of points in a consensus set to determine a correct model has been found
%       t         - error tolerance（错误容忍度） to determine whether a point is
%       compatible for the model,门限
% 
% Output:
%       correctPoints - indexes of point-pairs corrected matched determined by RANSAC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    %switch nargin
       %case 2
          % it =size(P2,1);
            %N = fix(0.4*size(P2,1));
           %t = 0.04;
        %case 3
           % N = fix(0.4*size(P2,1));
           % t = 0.04;
       % case 4
           % t = 0.04;
   % end
 
    correctPoints = [];
    pointNum = size(P2,1);
    history = [];
 
    for i = 1:it
        %initialize consensus set with 4 points randomly selected 
        S = randperm(pointNum);
        S = S(1:4);
        %check whether every 3 points are noncollinear
        if check_collineation(P1(S, :))
            it = it+1;
            continue;
        end
 
        %check whether the selected group haven't been used before
        historyNum = size(history,1);
        exist = 0;
        for j = 1:historyNum
            equal = 1;
            for k = 1:4
                if history(j,k) ~= S(k)
                    equal = 0;
                    break;
                end
            end
            if equal
                exist = 1;
                break;
            end
        end
        if exist
            %reject used group
            it = it+1;
            continue;
        else
            history(end+1, :) = S;
        end
    end
        %compute projective transformation matrix
        H=projectivematrix(P2(S,:),P1(S,:));
        %get full consensus set
        S =  get_consensus_set(P1, P2, S, H, t); 
        %check if the model(namely the matrix) is correct or better than that we've found
        if size(S,2) >N
            N = size(S,2);
            correctPoints = S;
        end
   
 
    %update projective transformation matrix变换矩阵
    H=projectivematrix(P2(correctPoints,:),P1(correctPoints,:));
    %get all points compatible with the best model
    correctPoints =  get_consensus_set(P1, P2, correctPoints, H, t);
end
 

 

