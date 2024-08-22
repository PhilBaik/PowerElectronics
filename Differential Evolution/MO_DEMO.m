function f = MO_DEMO(varargin)
    if length(varargin)==1
        x = varargin{1};
        x1 = x(1,1);
        x2 = x(1,2);
    else
        x1 = varargin{1};
        x2 = varargin{2};
    end
    
    f = (x1 - 2).^2 + (x2 - 3).^2 + sin(5 * x1) .* sin(5 * x2);
end