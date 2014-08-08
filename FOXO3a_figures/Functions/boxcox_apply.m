% BOXCOX: Box-Cox maximum-likelihood transformation to find the optimal power 
%         transformation of an empirical distribution to symmetry.  The transformed 
%         data are scaled to have the same mean and standard deviation as the original 
%         data.  Ignores missing data.
%
%     Syntax:  [xp,lambda,c,Lmax] = boxcox(x,{plotflag})
%
%         x =        vector of observations for a single variable.
%         plotflag = optional boolean flag producing, if true, histograms of
%                      original and transformed data and plot of L vs lambda 
%                      [default = 0].
%         -------------------------------------------------------------------------
%         xp =       corresponding vector of transformed variable.
%         lambda =   Box-Cox parameter.
%         c =        value added to data before transforming to ensure all positive 
%                      values.
%         Lmax =     max log-likelihood function value.
%

% RE Strauss, 2/10/97
%    9/3/99 - changed plot colors for Matlab v5.
%   12/8/02 - change fmin() to fminbnd(); output the constant c;
%             added error message for missing data.
%   1/2/03 -  added optional histograms.
%   1/3/03 -  ignore missing data; produce separate histograms rather than subplots. 

function xp = boxcox_apply(x,lambda,c)
  if (nargin < 3) help boxcox; return; end;
  
  norig = length(x);
  i = find(isfinite(x));                  % Remove missing data
  x = x(i);
  n = length(x);

  xmin = min(x);
  if (xmin <= 0)                          % Distribution must be positive
    x = x+c;
  end;
  
  if (abs(lambda) > eps)                
    xp = ((x.^lambda)-1)/lambda;
  else
    xp = log(x);
  end;

  return;
