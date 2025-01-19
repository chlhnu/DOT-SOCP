function obj = record_time(times, names)
%% Record the running time of algorithm
% Input
%   times: array of execution times [time1, time2, ...]
%   names: corresponding labels for the times {'name1', 'name2', ...}
% Output
%   obj: table containing execution times and their labels

obj = array2table(times);
obj.Properties.VariableNames = names;

end