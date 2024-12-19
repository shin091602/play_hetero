% find the point where dot{r2} == 0
function [value,isterminal,direction] = odestop_hetero_5(x,target_x)

value = x(1) - target_x; % dot{r}, the value that we want to be zero
isterminal = 1; % if isterminal == 1, the calcuration will stop when the event occurs.
direction = 0; % 2021-04-27  commment by Kyo Tetsu:  if it's 0 here, Apoaposis will also be outputed.
end