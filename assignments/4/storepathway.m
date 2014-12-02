function stop = storepathway(x, optimValues, state, varargin)
global pathway
pathway.x(end+1) = x(1);
pathway.y(end+1) = x(2);
pathway.z(end+1) = optimValues.resnorm;
stop = 0;