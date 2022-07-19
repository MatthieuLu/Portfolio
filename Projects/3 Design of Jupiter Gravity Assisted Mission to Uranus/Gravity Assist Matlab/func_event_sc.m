% Stops integration when s/c reaches a particular Planet's orbit.

function [value,isterminal,direction] = func_event_sc(r_Planet,~,p)

r = sqrt(p(1,1)^2+p(2,1)^2);

value = r-r_Planet;
isterminal = 1; % stop the integration
direction = 0; % 1 pos direc, 0 both direcs, -1 neg direc
