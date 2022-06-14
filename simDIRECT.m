function output = simDIRECT(fun,con,xl,xu,feps,fu,maxfn)
                         
% SIMDIRECT is an extension of the DIRECT Global Optimization Algorithm to
% handle multiple objectives, nonlinear constraints, and missing data.
%
% It solves the problem
%
%      Min f = {f_1(x),...,f_nobj(x)}  % find Pareto optimal points
% 
%  subject to:
% 
%      f(x) <= fu                      % objectives <= worst allowed values
%      g(x) <= 0                       % side constraints satisfied     
%      xl <= x <= xu                   % variables meet lower/upper bounds                
% 
% INPUTS:
% 
%      fun     function that can be invoked as 
%                  f = fun(x)
%              where
%                  x is the [1 x nvar] point to be evaluated 
%                  f is the [1 x nobj] vector of objective values at x
%              We assume all objectives are to be minimized, and that the
%              goal is to find points that approximate the Pareto set.
%              If the calculation of any objective fails, and so no value
%              is available, set the missing value to NaN (not a number).
%              Note that it is ok to have nobj=1, so simDIRECT can also
%              solve traditional single-objective optimization problems.
%
%      con     function that can be invoked as
%                  g = con(x)
%              where
%                  x is the [1 x nvar] point to be evaluated 
%                  g is the [1 x ncon] vector of constraint values at x
%              We assume all constraints are of the form g(x) <= 0.
%              If the calculation of any constraint fails, and so no value
%              is available, set the missing value to NaN (not a number).
%              If there are no side constraints, set con=[].
%                     
%      xl      [1 x nvar]  lower bounds on inputs
%
%      xu      [1 x nvar]  upper bounds on inputs
%
%      feps    [1 x nobj]  feps(i) = desired absolute accuracy objective i.
%                          Differences on objective i of feps(i) or less
%                          are considered insignificant.
%
%      fu      [1 x nobj]  worst allowed objective value ("nadir point"),  
%                          used in multiobjective optimization.  
%                          We require f <= fu in order for a point to be
%                          included in the Pareto set.  We also use fu
%                          to compute the "dominated hypervolume" after
%                          each iteration, which allows us to track
%                          progress converging to the Pareto front.  For
%                          single-objective optimization, set fu=Inf.
%      
%      maxfn   [ scalar ]  maximum allowed function evaluations
%
% OUTPUT:
%      
%      output   structure whose fields have the results from optimization:
%               output.x        [maxfn x nvar] all evaluated points
%               output.f        [maxfn x nobj] objectives at output.x
%               output.g        [maxfn x ncon] constraints at output.x
%               output.isPareto [maxfn x 1] isPareto(i)=true if point is
%                               feasible and non-dominated, i.e, on the
%                               Pareto front. It can happen that no 
%                               feasible point will be found.  If this 
%                               happens, then the Pareto Front will be 
%                               empty with isPareto(i)=false for all i. 
%
%               The Pareto front then has the points
%                   x_Pareto = output.x(output.isPareto,:)
%               with corresponding objectives and constraints given by
%                   f_Pareto = output.f(output.isPareto,:)
%                   g_Pareto = output.g(output.isPareto,:)
%               In single-objective optimization, x_Pareto will usually be
%               a single point unless there are points tied for the best
%               objective value.  If no feasible point was found, we will
%               have any(isPareto)=false.
%
% VERSION:  
%
%        This is version 1.1 (June 13, 2022).  Please direct questions or
%        comments on this software to Don Jones at donjon@umich.edu.
%
% IMPORTANT NOTE:
% 
%         To compute the dominated hypervolume, simDIRECT expects the user
%         to have downloaded the open-source code of Fonseca et al that
%         computes this quantity.  Instructions for getting this code and
%         compiling it into a Matlab MEX file are given below.  If you 
%         don't do this, simDIRECT will still run, but will just report the
%         dominated hypervolume as zero -- which means you won't be able to
%         track progress converging to the Pareto front.
% 
%         How to get the Hypervolume code for Matlab:
% 
%         * Go to:  http://lopez-ibanez.eu/hypervolume#download  
% 
%         * Download the Version 1.3 source code and unzip the tar file.  
%           You should get a directory called hv-1.3-src.
%           Open MATLAB and go to the hv-1.3-src directory.  
% 
%         * At the MATLAB command prompt, run:  
%           mex -DVARIANT=4 Hypervolume_MEX.c hv.c avl.c
% 
%         * This should give you a file called
%           Hypervolume_MEX.mexw64 (Windows) or
%           Hypervolume_MEX.mexmaci64 (Mac)
% 
%         * Put this file on the path when using simDIRECT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Evaluate the midpoint and get problem dimensions (nobj, ncon, nvar)
    
    x = (xl + xu)/2;
    fmid = fun(x);
    nobj = length(fmid);
    if isempty(con)
        ncon = 0;
    else
        gmid = con(x);
        ncon = length(gmid);
    end
    nvar = size(xl,2);

    % If nobj>1, make sure the size of both feps and fu are 1 x nobj, since
    % users sometimes think feps and fu are scalars.
    
    if nobj > 1 && (size(feps,2) ~= nobj || size(fu,2) ~= nobj)
        fprintf('ERROR:  The size of feps and fu should be 1 x nobj.\n');
        return
    end
    
    % Issue warning if nobj > 1 and Hypervolume_MEX is NOT installed
    
    if nobj > 1 && ~(exist('Hypervolume_MEX','file')==3)
        fprintf('WARNING:  You have multiple objectives, but\n');
        fprintf('Hypervolume_MEX is not found.  You can optimize,\n');
        fprintf('but you will not be able to track your progress\n');
        fprintf('in converging to the Pareto front.\n');
    end
    
    % simDIRECT assumes xl < xu.  You can't set xl(i)=xu(i) to hold
    % variable i constant.  Issue an error if this assumption is violated.
    
    if any(xl>=xu)
        fprintf('ERROR:  All variable lower bounds must be STRICTLY\n');
        fprintf('less than the upper bound.  You have violated this.\n');
        return
    end
    
    % Initialize / allocate space
    
    rc = zeros(maxfn,nvar);     % rectangle center points
    rl = zeros(maxfn,nvar);     % rectangle lower bounds
    ru = zeros(maxfn,nvar);     % rectangle upper bounds
    SmallestAllowedSide = 1e-10;% Smallest allowed normalized side length
    fathomed = false(maxfn,1);  % true if all rect sides < the above
    d = zeros(maxfn,1);         % rectangle center-vertex distances
    f = zeros(maxfn,nobj);      % objective values at rectangle centers
    g  = zeros(maxfn,ncon);     % constraint values at rectangle centers
    nsplit = zeros(maxfn,nvar); % nsplit(i,k) = times rect i split on dim k
    fsum = zeros(1,nobj);       % sum of observed rates of change for f
    gsum = zeros(1,ncon);       % sum of observed rates of change for g
    nfsum = zeros(1,nobj);      % number observed rates of change for f
    ngsum = zeros(1,ncon);      % number observed rates of change for g  
    tiebreak = zeros(1,nvar);   % tiebreak(k)=times ANY rect split on dim k
    iter = 0;                   % iteration counter
    nrec_last = 0;              % number rectangles at end of previous iter

    % Write output header
    
    if nobj > 1
        fprintf('         new    total   number     hyper-                       \n');
        fprintf(' iter   evals   evals   Pareto     volume      min(sumviolation)\n');
        fprintf('----------------------------------------------------------------\n');
    else
        fprintf('         new    total   number                                  \n');
        fprintf(' iter   evals   evals   Pareto     fmin        min(sumviolation)\n');
        fprintf('----------------------------------------------------------------\n');
    end
    
    % Form root rectangle
    
    nrec = 1; 
    rl(1,:) = zeros(1,nvar);
    ru(1,:) = ones(1,nvar);
    rc(1,:) = 0.5*ones(1,nvar);
    nsplit(1,:) = zeros(1,nvar);
    d(1) = CenterVertexDistance(nsplit(1,:));     
    fathomed(1) = false;
    f(1,:) = fmid; % (from our earlier evaluation of the midpoint)
    if ncon > 0
        g(1,:) = gmid;
    end
    
    % Iterate until reach max allowed function evaluations.  We break out
    % of the infinite loop below after reporting the iter with nrec=maxfn.
    
    while true
        
        iter = iter + 1;

        % Get status of search: isPareto, MinSumVio, dominated hypervolume
        % (hv) and fmin (best feasible objective value when nobj==1)

        [isPareto, hv, MinSumVio, fmin] = status(nrec,nobj,ncon,f,g,fu);    

        % Report results for this iteration 

        if nobj > 1
            fprintf(' %4d   %5d  %6d    %5d    %11.4e    %15.9e\n',...
                    iter,nrec-nrec_last,nrec,sum(isPareto),hv,MinSumVio);
        else
            fprintf(' %4d   %5d  %6d    %5d    %11.4e    %15.9e\n',...
                    iter,nrec-nrec_last,nrec,sum(isPareto),fmin,MinSumVio);
        end        

        % Record the number of rectangles as of the above progress report

        nrec_last = nrec;
                
        % Break out of loop if nrec = maxfn
        
        if nrec == maxfn
            break;
        end

        % Break out of loop if all rectangles are fathomed
        
        if all(fathomed(1:nrec))
            fprintf('All rectangles are fathomed!\n')
            break;
        end
        
        % Update average rates of change of objective and constraints
        
        fslope = fsum ./ max(nfsum, 1e-10);
        fslope = max(fslope, 1e-10);
        if ncon > 0
            gslope = gsum ./ max(ngsum, 1e-10);
            gslope = max(gslope, 1e-10);
        else
            gslope = [];
        end

        % Process failed runs (if any) and select rectangles for search
  
        selected = process_failed_runs_and_select(nrec,nobj,ncon,f,g,d,...
                               rc,fu,feps,fslope,gslope,isPareto,fathomed);
            
        % Trisect selected rectangles and sample center of children
        
        ns = length(selected);
        for i=1:ns
            
            % p = index of parent rectangle

            p = selected(i);
            
            % Find k = dimension for trisection.   This is the longest side
            % (or least split side) with smallest value of tiebreak 

            k = 1;
            for j=2:nvar
                if nsplit(p,j) < nsplit(p,k)
                    k = j;
                elseif nsplit(p,j)==nsplit(p,k) && tiebreak(j)<tiebreak(k)
                    k = j;
                end
            end
            
            % Increment tiebreak(k)

            tiebreak(k) = tiebreak(k) + 1;

            % Make short variables for parent bounds on dim k
            
            L = rl(p,k);
            U = ru(p,k);
            
            % Get distance from center of parent rectangle to center of
            % left/right children rectangles
            
            delta = 3^(-nsplit(p,k)-1);

            % Left child will be rectangle nrec+1 and right child will
            % be rectangle nrec+2

            n1 = nrec+1;
            n2 = nrec+2;

            % Form left child

            rl(n1,:) = rl(p,:);
            ru(n1,:) = ru(p,:);
            rc(n1,:) = rc(p,:);
            ru(n1,k) = (2*L+U)/3;
            rc(n1,:) = (rl(n1,:) + ru(n1,:))/2;
            fathomed(n1) = all( ru(n1,:)-rl(n1,:) < SmallestAllowedSide );
            nsplit(n1,:) = nsplit(p,:);
            nsplit(n1,k) = nsplit(p,k) + 1;
            d(n1) =  CenterVertexDistance(nsplit(n1,:));
                        
            % Evaluate left child
            
            f(n1,:) = fun( xl + rc(n1,:).*(xu-xl) );
            if ncon > 0
                g(n1,:) = con( xl + rc(n1,:).*(xu-xl) );
            end
            
            % Update sums of observed, absolute rates of change

            for t=1:nobj
                if ~isnan( f(n1,t) ) && ~isnan( f(p,t) )
                    fsum(t) = fsum(t) + abs(f(n1,t)-f(p,t))/delta;
                    nfsum(t) = nfsum(t) + 1;
                end
            end
            for t=1:ncon
                if ~isnan( g(n1,t) ) && ~isnan( g(p,t) )
                    gsum(t) = gsum(t) + abs(g(n1,t)-g(p,t))/delta;
                    ngsum(t) = ngsum(t) + 1;
                end
            end

            % Exit loop over selected rectangles if we have reached the 
            % max function evaluations
            
            if n1 == maxfn  
                nrec = n1;
                break
            end

            % Form right child as rectangle n2
            
            rl(n2,:) = rl(p,:);
            ru(n2,:) = ru(p,:);
            rc(n2,:) = rc(p,:);
            rl(n2,k) = (L+2*U)/3;
            rc(n2,:) = (rl(n2,:) + ru(n2,:))/2;
            fathomed(n2) = all( ru(n2,:)-rl(n2,:) < SmallestAllowedSide );
            nsplit(n2,:) = nsplit(p,:);
            nsplit(n2,k) = nsplit(p,k) + 1;
            d(n2) = CenterVertexDistance(nsplit(n2,:));
                        
            % Evaluate right child
            
            f(n2,:) = fun( xl + rc(n2,:).*(xu-xl) );
            if ncon > 0
                g(n2,:) = con( xl + rc(n2,:).*(xu-xl) );
            end
            
            % Update sums of observed, absolute rates of change

            for t=1:nobj
                if ~isnan( f(n2,t) ) && ~isnan( f(p,t) )
                    fsum(t) = fsum(t) + abs(f(n2,t)-f(p,t))/delta;
                    nfsum(t) = nfsum(t) + 1;
                end
            end
            for t=1:ncon
                if ~isnan( g(n2,t) ) && ~isnan( g(p,t) )
                    gsum(t) = gsum(t) + abs(g(n2,t)-g(p,t))/delta;
                    ngsum(t) = ngsum(t) + 1;
                end
            end

            % Exit loop over selected rectangles if we have reached the 
            % max function evaluations
                        
            if n2 == maxfn
                nrec = n2;
                break
            end
            
            % Modify parent rectangle to become center child
            
            rl(p,k) = (2*L+U)/3;
            ru(p,k) = (L+2*U)/3;
            fathomed(p) = all( ru(p,:)-rl(p,:) < SmallestAllowedSide );
            nsplit(p,k) = nsplit(p,k) + 1;
            d(p) = CenterVertexDistance(nsplit(p,:));
            
            % Update rectangle counter
            
            nrec = n2;
            
        end % end of loop over selected rectangles
     
    end % end of infinite loop (we break out when nrec=maxfn)
        
    % We are done!

    output.x = xl + rc(1:nrec,:) .* (xu-xl); % put points on original scale
    output.f = f(1:nrec,:);
    output.g = g(1:nrec,:);
    output.isPareto = isPareto(1:nrec,1);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = CenterVertexDistance(nsplit)

    % CENTERVERTEXDISTANCE computes the center-vertex distance based on 
    % how many times the rectangle has been split on each dimension.
    % The formula below is equivalent to 0.5 * norm( 3.^(-nsplit) ).
    % Using this formula ensures that two rectangles that have been split
    % the same total number of times will have EXACTLY the same d value.
    %
    % INPUT:
    %     nsplit   [1 x nvar] nsplit(k) is # times rec trisected on dim k
    % 
    % OUTPUT:
    %     d        [scalar  ] distance from rec center to any vertex 
    
    nvar = length(nsplit);
    T = sum(nsplit);
    j = mod(T,nvar);
    k = (T-j)/nvar;
    d = 0.5 * 3^(-k) * sqrt( j/9 + nvar - j );
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [isPareto, hv, MinSumVio, fmin] = status(nrec,nobj,ncon,f,g,fu)

    % STATUS finds the Pareto set - all points that have no missing values,
    % no constraint violations, and are not dominated on the objectives by
    % another feasible point.  If the Pareto set is non-empty, and if
    % nobj>1, it also computes the dominated hypervolume.  If there are no
    % feasible points, it computes the minimum, over all points, of the
    % sum-of-constraint-violations. Finally, when nobj==1, it computes the
    % best feasible objective function value (fmin).
    %
    % INPUTS
    %
    %     nrec      number of rectangles
    %     nobj      number of objectives 
    %     ncon      number of constraints
    %     f         [nrec x nobj] objective values, assume lower-is-better
    %     g         [nrec x ncon] side constraint values
    %     fu        [1    x nobj] objective upr bnd (require f(i)<=fu(i))
    %
    % OUTPUTS
    %
    %     IsPareto  [nrec x    1] true if rect i is feasible & nondominated
    %     hv        dominated hypervolume (or 0 if no feas pt or nobj==1)
    %     MinSumVio least sum of constraint vio (0 if have feasible pt)
    %     fmin      if nobj==1, fmin = best feas objective; else inf
    
    % We initialize all points that have no missing values and no
    % constraint violations as nondominated.  
    
    isPareto = true(nrec,1);
    for i=1:nrec
        if any( isnan(f(i,:)) )
            isPareto(i) = false; % missing objectives
        end
        if any( f(i,:) > fu )
            isPareto(i) = false; % objectives violate upper limit
        end
        if ncon > 0 && any( isnan(g(i,:)) )
            isPareto(i) = false; % missing constraints
        end
        if ncon > 0 && any( g(i,:) > 0 )
            isPareto(i) = false; % violate constraints
        end
    end
    
    % Of those points that pass the above filters, we set isPareto = false 
    % if the point is dominated on the objectives by some other
    % point that passes the filters
    
    for i=1:nrec
        if isPareto(i)
            for j=1:nrec
                if isPareto(j) && j ~= i
                    if all( f(j,:) <= f(i,:) ) && any( f(j,:) < f(i,:) )
                        isPareto(i) = false; % j dominates i
                    end
                end
            end
        end
    end
    
    % Get MinSumVio and fmin
    
    MinSumVio = Inf;
    fmin = Inf;
    for i=1:nrec
        if all( ~isnan(f(i,:)) )
            if ncon == 0
                % no missing values, so compute SumVio
                SumVio = sum( max(f(i,:)-fu,0) );
                MinSumVio = min(MinSumVio,SumVio);
                if SumVio == 0 && f(i,1) < fmin
                    fmin = f(i,1);
                end
            elseif all( ~isnan(g(i,:)) )
                % No missing values, so compute SumVio
                SumVio = sum( max(f(i,:)-fu,0) ) + sum( max(g(i,:),0) );
                MinSumVio = min(MinSumVio,SumVio);
                if SumVio == 0 && f(i,1) < fmin
                    fmin = f(i,1);
                end
            end
        end
    end

    % Get the dominated hypervolume if nobj>1 & the Pareto set is nonempty
    
    if nobj > 1 && any(isPareto) && exist('Hypervolume_MEX','file')==3
        hv = Hypervolume_MEX( f(isPareto,:), fu);
    else
        hv = 0;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selected = process_failed_runs_and_select(nrec,nobj,ncon,f,g,d,...
                                rc,fu,feps,fslope,gslope,isPareto,fathomed)

% PROCESS_FAILED_RUNS_AND_SELECT identifies points with missing data (NaN)
% and replaces the missing values with the values from the nearest point  
% that has no missing data.   It also creates another constraint whose 
% value is the distance to the nearest point with no missing data.  For
% points with no missing values, the constraint is zero.  For points with
% missing data, the constraint will be positive.   The average rate of
% change of this new constraint is set to 1.   Having filled in the missing
% values and added this new constraint, the routine then calls the SELECT
% function to identify rectangles for subdivision and sampling.
% Of course, if there are no points with missing data, the routine simply
% calls SELECT without adding the new constraint.
%
% INPUTS
%
%     nrec     number of rectangles
%     nobj     number of objectives
%     ncon     number of constraints
%     f        [nrec x nobj] objective values, assume lower-is-better
%     g        [nrec x ncon] side constraint values
%     d        [nrec x    1] center-vertex distances for rectangles
%     rc       [nrec x nvar] coordinates of rectangle center points
%     fu       [1    x nobj] worst allowed obj values (require f(i)<=fu(i))
%     feps     [1    x nobj] desired accuracy for objectives 
%     fslope   [1    x nobj] average absolute rate of change of objectives
%     gslope   [1    x ncon] average absolute rate of change of constraints
%     isPareto [nrec x    1] isPareto(i)=true if rec i cntr is nondominated
%     fathomed [nrec x    1] fathomed(i)=true if rec i too small to select
%
% OUTPUT
%
%     selected [ns   x    1] indices of selected rectangles (ns=#selected)

    % Set the center-vertex distance for any fathomed rectangle to zero.
    % As a result, these fathomed rectangles won't be able to dominate any
    % non-fathomed (positive size) rectangle.   They can't be selected
    % either because, with feps > 0, they can't beat any "epsilon-improved"
    % Pareto points -- even if they are one of them.

    d(fathomed) = 0;

    % Create "failed" where failed(i)=true if i-th point has NaN values

    if ncon > 0
        failed = any( isnan([f(1:nrec,:),g(1:nrec,:)]), 2 );
    else
        failed = any( isnan(f(1:nrec,:)), 2 );
    end
    
    if all(failed)
        
        % All runs failed; select all rectangles
        
        selected = (1:nrec)';
        return
        
    elseif all(~failed)
        
        % No failures, select as usual

         selected =  select(nrec,nobj,ncon,f,fu,feps,g,fslope,gslope,...
                            d,isPareto,fathomed);
                        
    else
        
        %   Preprocess data to account for failed runs.
        %   For each failed run, replace any NaN with values for f and g
        %   from the nearest non-failed point.  The nearest non-failed
        %   point is jmin and the distance is dminNofail(i).   
        %   For nonfailed points i, dminNoFail(i)=0.
        %   Having done this, add the constraint dminNofail <= 0 
        %   with gslope of 1.
        
        dminNofail = inf(nrec,1);
        for i=1:nrec
            if ~failed(i)
                dminNofail(i) = 0;
            else
                % Find closest non-failed point jmin and distance dmin
                dmin = inf;
                for j=1:nrec
                    if ~failed(j)
                        d = norm( rc(i,:) - rc(j,:) );
                        if d < dmin
                            dmin = d;
                            jmin = j;
                        end
                    end
                end
                % Store distance to nearest non-failed point
                dminNofail(i) = dmin;
                % Replace any NaN in f(i,:) and g(i,:) with values from
                % f and g for this closest point jmin
                for t=1:nobj
                    if isnan(f(i,t))
                        f(i,t) = f(jmin,t);
                    end
                end
                for t=1:ncon
                    if isnan(g(i,t))
                        g(i,t) = g(jmin,t);
                    end
                end
            end
        end        
        % add dist2nonfail <= 0 as another constraint with slope 1
        g      = [g(1:nrec,:) , dminNofail];
        gslope = [gslope, 1];
        ncon   = ncon + 1;
        
        % Call selection with the modified data
        
        selected =  select(nrec,nobj,ncon,f,fu,feps,g,fslope,gslope,...
                           d,isPareto,fathomed);
                        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selected = select(nrec,nobj,ncon,f,fu,feps,g,fslope,gslope,...
                           d,isPareto,fathomed)
    
% SELECT selects rectangles for subdivision based on multiple objectives,  
% side constraints, rectangle size, and average absolute 
% rates of change of the objectives and constraints.
%
% INPUTS
%
%     nrec     number of rectangles
%     nobj     number of objectives
%     ncon     number of constraints
%     f        [nrec x nobj] objective values, assume lower-is-better
%     fu       [1    x nobj] worst allowed obj values (require f(i)<=fu(i))
%     feps     [1    x nobj] desired accuracy for objectives 
%     g        [nrec x ncon] side constraint values
%     fslope   [1    x nobj] average absolute rate of change of objectives
%     gslope   [1    x ncon] average absolute rate of change of constraints
%     d        [nrec x    1] center-vertex distances for rectangles
%     isPareto [nrec x    1] isPareto(i)=true if rec i cntr is nondominated
%     fathomed [nrec x    1] fathomed(i)=true if rec i too small to select
%
% OUTPUT
%
%     selected [ns   x    1] indices of selected rectangles (ns=#selected)
    
    % If only one rectangle, select it and return
     
    if nrec == 1
        selected(1) = 1;
        return
    end

    % Make sure d is of length nrec
    
    d = d(1:nrec);
    
    % Initialize all non-fathomed rectangles as Potentially Pareto Optimal
    
    PPO = ~fathomed(1:nrec);
     
    % Find lower bound on alpha (amin) needed for each rect to be feasible.
    % Also find lower bound (alow) so feasible *and* nondominated  
    % by any "epsilon-improved" Pareto point.  This requires the lower
    % bound for some objective t to be feps(t) epsilon better than the
    % Pareto point's value of objective t.
    
    amin = zeros(nrec,1); % alpha needed to pass constraints
    alow = zeros(nrec,1); % alpha as above AND be strictly nondominated
        
    for r=1:nrec
                       
        % Find smallest alpa needed for r to be meet side constraints

        for k=1:ncon
            vio = max( 0, g(r,k) ) ;
            amin(r) = max( amin(r), vio/(gslope(k)*d(r)) );
        end

        % Make sure alpha high enough so f(r,:) <= fu

        for i=1:nobj
            vio  = max(0,f(r,i)-fu(i));
            amin(r) = max( amin(r), vio/(fslope(i)*d(r)) );  
        end

        % Finally make sure alpha big enough so strictly nondominated
        % by epsilon-improved Pareto points

        alow(r) = amin(r);
        for i=1:nrec
            if isPareto(i)
                if all( f(i,:)-feps <= f(r,:) )
                    % Epsilon-improved Pareto point i dominates r on obj's.
                    % Find lowest alpha so that r is feps(m) better than i
                    % on SOME objective m (so Strictly Non-Dominated)
                    % Below, "asnd" stands for "alpha so non dominated"
                    asnd = min((f(r,:)-f(i,:)+feps)./(fslope*d(r)));
                    alow(r) = max( alow(r), asnd );
                end
            end
        end
        
    end
    
    % Cycle through all rectangles r checking if there is an alpha
    % for which they are feasible and nondominated
   
    for r=1:nrec
        
        if ~PPO(r)
            % We already ruled out r as Potentially Pareto Optimal
            continue
        end
        
        % Initialize alpha range where nondominated to [alow(r),inf],
        % that is, the range where feasible and nondominated by
        % the "epsilon-improved" Pareto points

        validAlpha = set_of_intervals(alow(r),inf);

        % Cycle through all s ~= r, 
        % subtracting the interval of alpha values where s dominates r
       
        for s=1:nrec
            if s ~= r

                if d(s) > d(r) 
                    % s is bigger than r, so s will dominate r if alpha 
                    % is big enough.  So we find the alpha just large
                    % enough so that (1) s is feasible and (2)
                    % s dominates r on the objectives.   Call this 
                    % value a.  We substract (a,inf].
                    a1 = amin(s);
                    a2 = max((f(s,:)-f(r,:))./(fslope*(d(s)-d(r))));
                    a  = max( a1, a2 );
                    validAlpha = Subtract(validAlpha,a,inf);
                    if validAlpha.n == 0
                        % validAlpha is now empty, there is no alpha where
                        % r is feasible and nondominated.
                        PPO(r) = false;
                        % break out of loop checking if some s dominates r
                        break;
                    end
                elseif d(s) == d(r) 
                    % s has the same size as r.
                    % Differences between s & r stay constant as alpha
                    % changes.  If s dominates r on objectives, then r
                    % will be dominated by s as soon as alpha is large
                    % enough that s is feasible.
                    if all(f(s,:)<=f(r,:)) && any(f(s,:)<f(r,:))
                        % s dominates r on objectives; as soon as s is
                        % feasible, s dominates, so s dominates if 
                        % alpha > amin(s).  We subtract (amin(s), inf]:
                        validAlpha = Subtract(validAlpha,amin(s),inf);
                        if validAlpha.n == 0
                            % validAlpha is now empty, there is no alpha 
                            % where r is feasible and nondominated.
                            PPO(r) = false;
                            % break out of loop over s
                            break;
                        end
                    end
                elseif d(s) < d(r)
                    % s is smaller than r.
                    % So it only imposes a potential restriction on 
                    % alpha if it already dominates on the objectives.
                    % If it dominates on the objectives, it will rule
                    % out some alphas if there exist alphas such that:
                    % (1) alpha is big enough so s is feasible and
                    % (2) it is small enough that s still dominates r
                    % on the objectives
                    if all(f(s,:)<=f(r,:)) && any(f(s,:)<f(r,:))
                        % s dominates r on objectives as long as 
                        % amin(s) < alpha < b where
                        b=min((f(r,:)-f(s,:))./(fslope*(d(r)-d(s))));
                        if b > amin(s)
                            validAlpha=Subtract(validAlpha,amin(s),b);
                            if validAlpha.n == 0
                                % validAlpha is now empty, there is no 
                                % alpha where r is feasible and 
                                % nondominated.
                                PPO(r) = false;
                                % break out of loop over s
                                break;
                            end
                        end
                    end                        
                end 

            end
        end
        
    end
    
    % Set the selected rectangles to list of r with PPO(r)=true
    
    rnum = (1:nrec)';
    selected = rnum(PPO);
    
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = set_of_intervals(a0,b0)

    % SET_OF_INTERVALS initializes a structure "obj" to contain a single
    % interval [a0,b0].   Later, the structure "obj" will be operated on by
    % invoking "obj = Subtract(obj,a,b)" which will replace each interval
    % in obj with what is left after subtracting [a,b].  Subtraction can 
    % result in fewer intervals, the same number of intervals, or more 
    % intervals.  It can also happen that no intervals remain, in which 
    % case we will have obj.n==0.
    
    nmax     = 10;
    obj.L    = zeros(nmax,1);   % If more intervals than nmax end up being
    obj.U    = zeros(nmax,1);   % needed, Matlab will automatically resize
    obj.n    = 1;
    obj.L(1) = a0;
    obj.U(1) = b0;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = Subtract(obj,a,b)

    % SUBTRACT subtracts [a,b] from all intervals in the set contained in
    % the structure "obj", returning the modified structure.

    i = 1;
    while i <= obj.n

        % Do the subtraction of [a,b] from [L(i),U(i)]

        if a <= obj.L(i) && b >= obj.L(i) &&  b < obj.U(i)
            % Relative order:  a L b U
            % Replace [L,U] with [b,U]
            obj.L(i) = b;
            % Move on to the next interval in the set
            i = i + 1;
            continue
        elseif a <= obj.L(i) && b >= obj.U(i)
            % Relative order:  a L U b 
            % Delete [L,U]
            if i == obj.n
                % interval i is last, delete by decrementing n
                obj.n = obj.n - 1;
                % No need to increment i.
                % We now have i > obj.n so loop will terminate.
            else
                % Interval i is not last; delete it by replacing it
                % with the n-th interval and decrementing n.
                obj.L(i) = obj.L(obj.n);
                obj.U(i) = obj.U(obj.n);
                obj.n    = obj.n - 1;
                % Do not increment i as we need to now process the 
                % interval which has been copied to slot i
            end
            continue
        elseif a > obj.L(i) && b < obj.U(i)
            % Relative order:  L a b U  
            % Add new interval [b,U] and replace [L,U]  with [L,a]
            % Note that, if obj.n is now bigger than the original allocated
            % size of obj.L and obj.U, Matlab automatically resizes those
            % arrays.
            obj.n = obj.n +  1;
            obj.L(obj.n) = b;
            obj.U(obj.n) = obj.U(i);
            obj.U(i) = a;
            % Move on to the next interval in the set
            i = i + 1;
            continue
        elseif a > obj.L(i) && b >= obj.U(i) && a < obj.U(i)
            % Relative order:  L  a U b  
            % Replace [L,U] with [L,a]
            obj.U(i) = a;
            % Move on to the next interval in the set
            i = i + 1;
        else
            % [a,b] has no intersecton with [ L(i), U(i) ] so
            % nothing need be done.  Move on to the next interval.
            i = i + 1;
        end        
        
    end % of while i <= obj.n

end
