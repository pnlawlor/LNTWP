%% Warped signals

signal_true = [0 0 0 .1 .2 .3 .4 .5 .4 .3 .2 .1 0 0 0 0];

signal_warped = [0 0 0 0 0 0 .1 .2 .3 .4 .5 .4 .3 .2 .1 0]; 


%% Match them up 

sim_func = @(s1,s2) exp(-(s1 - s2).^2);


d1 = 1;
d2 = 1;

for i = 1:length(signal_true)
    options = [ ...
        d1 + 1, d2; ...     % over
        d1, d2 + 1; ...     % up
        d1 + 1, d2 + 1 ...  % over and up
        ];
        
    scores = sim_func( ...
        signal_true(options(:,1)), ...gi
        signal_warped(options(:,2)));
    
    [~, choice] = max(scores);
    choice = choice(1); % if there are multiple good options, choose "over"
    
end