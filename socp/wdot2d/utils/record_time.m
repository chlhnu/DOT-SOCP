function obj = record_time(times, names)
%% Record the running time of algorithm

obj = array2table(times);
obj.Properties.VariableNames = names;

end