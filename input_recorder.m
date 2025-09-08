classdef input_recorder < handle
properties
%list that store function inputs
input_list;
end
methods
function obj=input_recorder()
obj.input_list = [];
end
%resets input_list
function clear_input_list(obj)
obj.input_list = [];
end
%returns input_list
function input_list = get_input_list(obj)
input_list = obj.input_list;
end
%creates a version of function fun(x) that will store
%x into the input_list whenever fun is called
function f_record = generate_recorder_fun(obj,fun)
f_record = @(x) obj.record_wrapper(fun,x);
end
function varargout = record_wrapper(obj,fun,x)
nOutputs = nargout;
varargout = cell(1,nOutputs);
if nargout == 1
a = fun(x);
varargout{1} = a;
end
if nargout == 2
[a, b] = fun(x);
varargout{1} = a;
varargout{2} = b;
end
if nargout == 3
[a, b, c] = fun(x);
varargout{1} = a;
varargout{2} = b;
varargout{3} = c;
end
obj.input_list(end+1) = x;
end
end
end