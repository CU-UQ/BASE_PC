% Solves BPDN
% code modified from spgl1 https://www.math.ucdavis.edu/~mpf/spgl1/

% -----
% x = bpdn(A, b, sigma)
% -----
% Input
% -----
% A = lhs matrix
% b = rhs vector
% sigma = regularization parameter that may be cross-validated
% ------
% Output
% ------
% x = coefficient vector

% % This code modifes spgl1.m and spg_bpdn.m as part of SPGL1 http://www.cs.ubc.ca/~mpf/spgl1/index.html having the following copyright
%   SPGL1 copyright
%   -------------------------------------------------------------------
%   Copyright (C) 2007 Ewout van den Berg and Michael P. Friedlander,
%   Department of Computer Science, University of British Columbia, Canada.
%   All rights reserved. E-mail: <{ewout78,mpf}@cs.ubc.ca>.
%
%   SPGL1 is free software; you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as
%   published by the Free Software Foundation; either version 2.1 of the
%   License, or (at your option) any later version.
%
%   SPGL1 is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General
%   Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with SPGL1; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
%   USA
%   ----------------------------------------------------------------------
function x = bpdn(A,b,sigma)
    m = length(b);
    x = A\b;
    
    %%%%%%%%%%%%%%%%%%%%
    % function options %
    %%%%%%%%%%%%%%%%%%%%
    tau =  0; % Initialization of regularization parameter tau.
    maxIts        = 1000*m; % 100 times standard size. Some problems require large numbers of iterations. However some problems seem to require very large iteration counts
    nPrevVals     = 6; %  twice standard size. Doesn't cost much. Can improve accuracy/speed.
%    bpTol         = 1e-06; % We choose not to do this check (overlap with SUBOPTIMAL_BP SOL check)
    lsTol         = 0.5e-09; % 0.0005 standard size.
    optTol        = 1e-07; % reduced by 0.001
    decTol        = 1e-07; % reduced by 0.001
    stepMin       = 1e-16;
    stepMax       = 1e+05;
    maxLineErrors = 20;     % Maximum number of line-search failures. Twice typical value

    
    % Appropriate relativization of sigma for b.
    b_norm = norm(b);
    sigma = sigma*b_norm; % Normalizes sigma based on b_norm to be compatible across multiple rhs.
    %----------------------------------------------------------------------
    % Initialize local variables.
    %----------------------------------------------------------------------
    iter          = 0;  % Total SPGL1 iterations.
    lastFv        = -inf(nPrevVals,1);  % Last m function values.   
    stat          = false;
    testUpdateTau = 0;              % Previous step did not update tau
    
    % Exit conditions (constants).
    EXIT_ROOT_FOUND    = 1;
    %EXIT_BPSOL_FOUND   = 2;
    EXIT_LEAST_SQUARES = 3;
    EXIT_ITERATIONS    = 5;
    EXIT_LINE_ERROR    = 6;
    EXIT_SUBOPTIMAL_BP = 7;

    % Project the starting point and evaluate function and gradient.
    x = project(x,tau); % Soft Thresholding with parameters determined by x,tau
    r = b - A*x;  % Residual vector
    g = - A'*r;  %Correlations with residuals
    f = r'*r / 2;  % Residual Functional for many tests

    % Required for nonmonotone strategy.
    lastFv(1) = f;
    fBest     = f;
    fOld      = f;

    % Compute projected gradient direction and initial steplength.
    dx     = project(x - g, tau) - x;
    dxNorm = norm(dx,inf);
    if dxNorm < (1 / stepMax)
       gStep = stepMax;
    else
       gStep = min(stepMax, max(stepMin, 1/dxNorm)); % stepMin > stepMax is possible.
    end

    while 1 % main loop

        % Compute quantities needed for log and exit conditions.
        gNorm   = norm(g,inf);
        r_norm   = norm(r, 2);
        gap     = r'*(r-b) + tau*gNorm;
        rGap    = abs(gap) / max(1,f);
        aError1 = r_norm - sigma;
        aError2 = f - sigma^2 / 2;
        rError1 = abs(aError1) / max(1,r_norm);
        rError2 = abs(aError2) / max(1,f);
        
        % If rNorm is small, least squares solution. Can occur
        % if oversampled or exceptional approximation.
        if gNorm <= lsTol * r_norm % Occurs if no column has significant correlation with residual
            stat = EXIT_LEAST_SQUARES;
        end
        
        if rGap <= max(optTol, rError2) || rError1 <= optTol
            % The problem is nearly optimal for the current tau.
            % Check optimality of the current root.
            %if r_norm  <=  bpTol, stat=EXIT_BPSOL_FOUND; end  % Resid minim'zd -> BP sol.
            if rError1 <=  optTol, stat=EXIT_ROOT_FOUND; end  % Found approx root.
            if r_norm <=  sigma, stat=EXIT_SUBOPTIMAL_BP; end  % Found suboptimal BP sol.
        end
        
        testRelChange1 = (abs(f - fOld) <= decTol * f);
        testRelChange2 = (abs(f - fOld) <= 1e-1 * f * (abs(r_norm - sigma)));
        testUpdateTau  = ((testRelChange1 && r_norm >  2 * sigma) || ...
            (testRelChange2 && r_norm <= 2 * sigma)) && ~stat && ~testUpdateTau;
        
        if testUpdateTau
            % Update tau.
            tauOld   = tau;
            tau      = max(0,tau + (r_norm * aError1) / gNorm);           
            if tau < tauOld % if decrease, project
                x = project(x,tau);
            end
        end
        
        if ~stat  &&  iter >= maxIts % iteration limit
            stat = EXIT_ITERATIONS;
        end
        
        if stat
            %disp(stat);
            break; 
        end % If done, break.
        
        iter = iter + 1;
        xOld = x;  fOld = f;  gOld = g; % in case search fails
        %---------------------------------------------------------------
        % Projected gradient step and linesearch.
        %---------------------------------------------------------------
        [f,x,r,lnErr] = spgLineCurvy(x,gStep*g,max(lastFv),A,b,tau);
        if lnErr % Projected backtrack failed. Retry with feasible dir'n linesearch.
            x    = xOld;
            f    = fOld;
            dx   = project(x - gStep*g, tau) - x;
            gtd  = g'*dx;
            [f,x,r,lnErr] = spgLine(f,x,dx,gtd,max(lastFv),A,b);
        end
        if lnErr % Failed again. Revert to previous iterates and damp max BB step.            
            x = xOld;
            f = fOld;
            if maxLineErrors <= 0
                stat = EXIT_LINE_ERROR;
            else
                stepMax = stepMax / 10;
                maxLineErrors = maxLineErrors - 1;
            end
        end
        
        %---------------------------------------------------------------
        % Update gradient and compute new Barzilai-Borwein scaling.
        %---------------------------------------------------------------
        if (~lnErr)
            g    = - A'*r;
            s    = x - xOld;
            y    = g - gOld;
            sts  = s'*s;
            sty  = s'*y;
            if   sty <= 0,  gStep = stepMax;
            else            gStep = min( stepMax, max(stepMin, sts/sty) );
            end
        else
            gStep = min( stepMax, gStep );
        end
        
        %------------------------------------------------------------------
        % Update function history.
        %------------------------------------------------------------------
        if f > sigma^2 / 2 % Don't update if superoptimal.
            lastFv(mod(iter,nPrevVals)+1) = f;
            if fBest > f
                fBest = f;
            end
        end        
    end
end

function [fNew,xNew,rNew,err] = spgLineCurvy(x,g,fMax,A,b,tau)
    % Projected backtracking linesearch.
    % On entry,
    % g  is the (possibly scaled) steepest descent direction.
    EXIT_CONVERGED  = 0;
    EXIT_ITERATIONS = 1;
    EXIT_NODESCENT  = 2;
    gamma  = 1e-4;
    maxIts = 10;
    step   =  1;
    sNorm  =  0;
    scale  =  1;      % Safeguard scaling. Reduced if needed.
    nSafe  =  0;      % No. of safeguarding steps.
    iter   =  0;
    n      =  length(x);

    while 1

        % Evaluate trial point and function value.
        xNew     = project(x - step*scale*g, tau);
        rNew     = b - A*xNew;
        fNew     = rNew'*rNew / 2;
        s        = xNew - x;
        gts      = scale * real(g' * s); % real suffices
        if gts >= 0 % can happen in poorly conditioned problem
           err = EXIT_NODESCENT;
           break
        end

        if fNew < fMax + gamma*step*gts % check convergence
           err = EXIT_CONVERGED;
           break
        elseif iter >= maxIts % check iteration limit
           err = EXIT_ITERATIONS;
           break
        end

        % New linesearch iteration.
        iter = iter + 1;
        step = step / 2;

        % If same point after projection, drastically damp next search.
        sNormOld  = sNorm;
        sNorm     = norm(s) / sqrt(n);
        if abs(sNorm - sNormOld) <= 1e-6 * sNorm
           gNorm = norm(g) / sqrt(n);
           scale = sNorm / gNorm / (2^nSafe);
           nSafe = nSafe + 1;
        end
    end
end

function [fNew,xNew,rNew,err] = spgLine(f,x,d,gtd,fMax,A,b)
    % Nonmonotone linesearch.
    EXIT_CONVERGED  = 0;
    EXIT_ITERATIONS = 1;
    maxIts = 10;
    step   = 1;
    iter   = 0;
    gamma  = 1e-4;
    gtd    = -abs(gtd); % 03 Aug 07: If gtd is complex,
                        % then should be looking at -abs(gtd).
    while 1

        % Evaluate trial point and function value.
        xNew = x + step*d;
        rNew = b - A*xNew;
        fNew = rNew'*rNew / 2;

        % Check exit conditions.
        if fNew < fMax + gamma*step*gtd  % Sufficient descent condition.
           err = EXIT_CONVERGED;
           break
        elseif  iter >= maxIts           % Too many linesearch iterations.
           err = EXIT_ITERATIONS;
           break
        end

        % New linesearch iteration.
        iter = iter + 1;

        % Safeguarded quadratic interpolation.
        if step <= 0.1
           step  = step / 2;
        else
           tmp = (-gtd*step^2) / (2*(fNew-f-step*gtd));
           if tmp < 0.1 || tmp > 0.9*step || isnan(tmp)
              tmp = step / 2;
           end
           step = tmp;
        end
    end
end

function x = project(b,tau)   

    % Get sign of b and set to absolute values
    s = sign(b);
    b = abs(b);

   % Initialization
   n     = length(b);
   b_norm = norm(b,1);
   
    % Check for quick exit.
   if (tau >= b_norm), x = b; return; end   
   x     = zeros(n,1); % Initialize x here.      
   if (tau <  eps), return; end

   % Preprocessing (b is assumed to be >= 0)
   [b,idx] = sort(b,'descend'); % Descending.

   csb       = -tau;
   alphaPrev = 0;
   for j= 1:n
      csb       = csb + b(j);
      alpha     = csb / j;
   
      % We are done as soon as the constraint can be satisfied
      % without exceeding the current minimum value of b
      if alpha >= b(j)
         break;
      end   
      alphaPrev = alpha;
   end

   % Set the solution by applying soft-thresholding with
   % the previous value of alpha
   x(idx) = max(0,b - alphaPrev);
   
   % Restore signs in x
   x = x.*s;
end
