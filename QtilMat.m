classdef QtilMat
  %QtilMat class
  %  A QtilMat object is used to represent a matrix Qtil,
  %       where Q is accessed via function evaluations for both
  %       matrix-vector and matrix-transpose-vector multiplications
  %     Qtil = [Q2 0 ] + [W] Q1 [W' I]
  %            [0  0 ]   [I]
  %
  %  The QtilMat class has inputs:
  %     Q1  - matrix or object
  %     Q2  - matrix or object
  %     W   - matrix or object
  %
  %
  %   and is based on a structure with the following fields:
  %     Q1
  %     Q2
  %     W
  %     transpose - indicates if the matrix has been transposed
  %
  %  Calling Syntax:
  %
  %    P = QtilMat    (returns object with empty fields)
  %    P = QtilMat(QtilMat) (returns QtilMat)
  %    P = QtilMat(Q1, Q2, W)
  %
  % J. Chung, 1/2/2020
  
  properties
    Q1
    Q2
    W
    transpose
  end
  
  methods
    function P = QtilMat(varargin) % Constructor
      switch nargin
        case 0
          P.transpose = false;
          P.Q1 = [];
          P.Q2 = [];
          P.W = [];
        case 1
          if isa(varargin{1}, 'QtilMat')
            P = varargin{1};
          else
            error('Incorrect input arguments')
          end
        otherwise
          P.transpose = false;
          if nargin >= 3
            P.Q1 = varargin{1};
            P.Q2 = varargin{2};
            P.W = varargin{3};
          else
            error('Incorrect number of input arguments')
          end
      end
    end
    
    function P = ctranspose(P) % Overload transpose
      P.transpose = not(P.transpose); % switches booleen transpose flag
    end
    
    function y = mtimes(arg1, arg2) % Overload matrix vector multiplication
      
      if isa(arg1,'QtilMat') && isa(arg2,'double')
        %   Implement A*s and A'*s for QtilMat object A.
        if isscalar(arg2)
          error('Matrix times a scalar not implemented yet')
        elseif isvector(arg2)          
          if arg1.transpose 
            % Matrix transpose times vector
            [n,l] = size(arg1.W);
            y1 = zeros(n+l,1);
            y1(1:n) = arg1.Q2'*arg2(1:n);
            
            tmp = arg1.Q1'*(arg1.W'*arg2(1:n) + arg2(n+1:end));
            y2 = [arg1.W*tmp;tmp];
            y = y1 + y2;            
          else
            % Matrix times vector
            [n,l] = size(arg1.W);
            y1 = zeros(n+l,1);
            y1(1:n) = arg1.Q2*arg2(1:n);
            
            tmp = arg1.Q1*(arg1.W'*arg2(1:n) + arg2(n+1:end));
            y2 = [arg1.W*tmp;tmp];
            y = y1 + y2;
          end  
        end
      else
        error('Multiplication is not implemented.')
      end
      end
    
    function varargout = size(A, dim) % Overload size
      [n,l] = size(A.W);
      d = [n+l, n+l];
      if nargin == 2
        d = d(dim);
      end
      
      if nargout == 1 || nargout == 0
        varargout{1} = d;
      else
        for i = 1:length(d)
          varargout{i} = d(i);
        end
      end
      
    end

    function n = sizen(A)
      n = size(A.W,1);      
    end;

    function p = sizep(A)
      p = size(A.W,2); 
    end;

  end % methods
end % classdef

