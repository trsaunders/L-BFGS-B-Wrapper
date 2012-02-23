function [x,fx,res] = lbfgsb(f,x0,lb, ub, varargin)
	% L-BGFS-B

	% check that the bounds vectors are the same length as input
	bc = @(x) isnumeric(x) && length(x) == length(x0);
	p = inputParser;
	o = inputParser;
	p.addRequired('f', @(x) isa(x, 'function_handle'));
	p.addRequired('x0', @isnumeric);
	p.addRequired('lb', bc);
	p.addRequired('ub', bc);
	p.parse(f,x0,lb,ub);
	
	% optional params
	o.addParamValue('maxiter', 100, @(x) isnumeric(x) && x > 0);
	o.addParamValue('factr', 1e7, @isnumeric);
	o.addParamValue('m', 5, @isnumeric);
	o.addParamValue('pgtol', 1e-5, @isnumeric);
	%o.addParamValue('verbose', 1, @isnumeric);
	o.parse(varargin{:});

	% Call mex function
	[x] = lbfgsb_(f, x0(:), lb(:), ub(:), o.Results);
end
