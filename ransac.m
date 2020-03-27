% RANSAC - Robustly fits a model to data with the RANSAC algorithm
%
% Usage:
%
% [M, inliers] = ransac(x, fittingfn, distfn, degenfn, s, t, maxDataTrials, maxTrials)
%
% Arguments:
%     x         - Data sets to which we are seeking to fit a model M
%                 It is assumed that x is of size [d x Npts]
%                 where d is the dimensionality of the data and Npts is
%                 the number of data points.
%
%     fittingfn - Handle to a function that fits a model to s
%                 data from x.  It is assumed that the function is of the
%                 form:
%                    M = fittingfn(x)
%                 Note it is possible that the fitting function can return
%                 multiple models (for example up to 3 fundamental matrices
%                 can be fitted to 7 matched points).  In this case it is
%                 assumed that the fitting function returns a cell array of
%                 models.
%                 If this function cannot fit a model it should return M as
%                 an empty matrix.
%
%     distfn    - Handle to a function that evaluates估计 the
%                 distances距离 from the model to data x.
%                 It is assumed that the function is of the form:
%                    [inliers, M] = distfn(M, x, t)
%                 This function must evaluate the distances between points
%                 and the model returning the indices of elements in x that
%                 are inliers内点, that is, the points that are within distance
%                 't' of the model.  Additionally, if M is a cell array of
%                 possible models 'distfn' will return the model that has the
%                 most inliers.  If there is only one model this function
%                 must still copy the model to the output.  After this call M
%                 will be a non-cell object representing only one model.
%
%     degenfn   - Handle to a function that determines whether a
%                 set of datapoints will produce a degenerate model.
%                 This is used to discard丢弃 random samples that do not
%                 result in useful models.
%                 It is assumed that degenfn is a boolean function of
%                 the form:
%                    r = degenfn(x)
%                 It may be that you cannot devise想到 a test for degeneracy in
%                 which case you should write a dummy假的 function that always
%                 returns a value of 1 (true) and rely on 'fittingfn' to return
%                 an empty model should the data set be degenerate.
%
%     s         - The minimum number of samples采样最小值 from x required by
%                 fittingfn to fit a model.
%
%     t         - The distance threshold距离阈值 between a data point and the model
%                 used to decide whether the point is an inlier or not.
%
%     maxDataTrials - Maximum number of attempts尝试最大值 to select a non-degenerate
%                     data set. This parameter is optional and defaults to 100.
%
%     maxTrials - Maximum number of iterations迭代次数最大值. This parameter is optional and
%                 defaults to 1000.
%
%
% Returns:
%     M         - The model having the greatest number of inliers.内点数量
%     inliers   - An array of indices of the elements of x that were
%                 the inliers for the best model.
%
% For an example of the use of this function see RANSACFITHOMOGRAPHY or
% RANSACFITPLANE

% References:
%    M.A. Fishler and  R.C. Boles. "Random sample concensus: A paradigm
%    for model fitting with applications to image analysis and automated
%    cartography". Comm. Assoc. Comp, Mach., Vol 24, No 6, pp 381-395, 1981
%
%    Richard Hartley and Andrew Zisserman. "Multiple View Geometry in
%    Computer Vision". pp 101-113. Cambridge University Press, 2001

% Copyright (c) 2003-2006 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
% May      2003 - Original version
% February 2004 - Tidied up.
% August   2005 - Specification of distfn changed to allow model fitter to
%                 return multiple models from which the best must be selected
% Sept     2006 - Random selection of data points changed to ensure duplicate
%                 points are not selected.
% February 2007 - Jordi Ferrer: Arranged warning printout.
%                               Allow maximum trials as optional parameters.
%                               Patch the problem when non-generated data
%                               set is not given in the first iteration.

function [M, inliers] = ransac(x, fittingfn, distfn, degenfn, s, t, ...
                               maxDataTrials, maxTrials)
  % Test number of parameters
  error ( nargchk ( 6, 8, nargin ) ); %nargchk确认输入参数的数量
  error ( nargoutchk ( 2, 2, nargout ) ); %nargoutchk确认输出参数的数量

  [rows, npts] = size(x);                 

  p = 0.99;    % Desired要求 probability机率 of choosing at least one sample
               % free from outliers
               % 要求至少有一个好样本的概率阈值为99%

  if nargin < 8; maxTrials = 1000;    end; % Maximum number of trials before we give up.
                                           % 在我们放弃前实验次数最大值
  if nargin < 7; maxDataTrials = 100; end; % Max number of attempts to select a non-degenerate
                                           % data set.
                                           % 尝试去选择一个非退化数据集的最大次数 
  bestM = NaN;      % Sentinel value allowing detection发现 of solution failure.
  trialcount = 0;
  bestscore =  0;
  N = 1;            % Dummy initialisation for number of trials.

  while N > trialcount && trialcount <= maxTrials;
    % Select at random s datapoints to form a trial model,
    % M.随机选择s个数据点来形成一个实验模型，M
    % In selecting these points we have to check that they are not in
    % a degenerate configuration. 我们必须在选择的点中检查他们并不是在一个退化的结构中。
    degenerate = 1;
    count = 1;
    while degenerate && count <= maxDataTrials;
      % Generate s random indicies in the range 1..npts
      % 在范围1...npts中形成s个随机指标。
      % (If you do not have the statistics toolbox, or are using Octave八个一组的,
      % use the function RANDOMSAMPLE from my webpage)
      ind = randsample(npts, s); %随机采样

      % Test that these points are not a degenerate configuration.
      degenerate = feval(degenfn, x(:,ind));
	    
      if ~degenerate;
        % Fit model to this random selection of data points.
        % Note that M may represent a set of models that fit the data in
        % this case M will be a cell array of models
        M = feval(fittingfn, x(:,ind));
		
        % Depending on your problem it might be that the only way you
        % can determine whether a data set is degenerate or not is to
        % try to fit a model and see if it succeeds.  If it fails we
        % reset degenerate to true.
        if isempty(M), degenerate = 1; end;
      end

      % Safeguard保护 against being stuck in this loop forever
      count = count + 1;
    end

    if degenerate;
      warning ( 'MATLAB:ransac:Output', ...
                'Unable to select a nondegenerate data set!' );
      trialcount = trialcount + 1;
      break;
    end

    % Once we are out here we should have some kind of model...        
    % Evaluate distances估计距离 between points and model returning the indices指标
    % of elements in x that are inliers. 
    % 估计点和模型之间的距离，返回x中的元素指标，它们是内点。
    % Additionally, if M is a cell array of possible models 'distfn' will return the model that has
    % the most inliers.  After this call M will be a non-cell object
    % representing only one model.
    [inliers, M] = feval(distfn, M, x, t);

    % Find the number of inliers to this model.
    ninliers = length(inliers);

    if ninliers > bestscore;   % Largest set of inliers so far...
      bestscore = ninliers;    % Record data for this model
      bestinliers = inliers;
      bestM = M;

      % Update estimate of N, the number of trials to ensure we pick,
      % with probability p, a data set with no outliers.
      fracinliers =  ninliers/npts;
      pNoOutliers = 1 -  fracinliers^s;
      pNoOutliers = max(eps, pNoOutliers);  % Avoid division by -Inf
      pNoOutliers = min(1-eps, pNoOutliers);% Avoid division by 0.
      N = log(1-p)/log(pNoOutliers);
    end

    trialcount = trialcount + 1;  
  end

  % Safeguard against being stuck in this loop forever
  if trialcount > maxTrials;
    warning ( 'MATLAB:ransac:Output', ...
              sprintf ( 'Ransac reached the maximum number of %d trials!', maxTrials ) );
  end   

  if ~isnan(bestM)   % We got a solution
    M = bestM;
    inliers = bestinliers;
  else
    M = [];
    inliers = [];
    warning ( 'MATLAB:ransac:Output', '\nRansac was unable to find a useful solution!' );
  end
end