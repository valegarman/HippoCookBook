function [n_jumpy, led_pos] = led_speed_filter(led_pos, max_pix_per_sample, led)

% Filters out short runs of data caused by tracker picking up an incorrect distant point.
% Resets led_pos_x to NaN

% Find first OK point & assume that this is OK
ok_pos = find( ~isnan(led_pos(:,led,1)) );

if length(ok_pos) < 2
    warning(' < 2 tracked points for led %d\n', led);
end
mpps_sqd = max_pix_per_sample^2;

n_jumpy = 0;
prev_pos = ok_pos(1);
for i = 2:length(ok_pos)
    pos = ok_pos(i);
    % Get speed of shift from prev_pos in pixels per sample (squared)
    pix_per_sample_sqd = ((led_pos(pos,led,1)-led_pos(prev_pos,led,1))^2+...
        (led_pos(pos,led,2)-led_pos(prev_pos,led,2))^2)/(pos-prev_pos)^2;
    if pix_per_sample_sqd > mpps_sqd
        led_pos(pos, led, :) = NaN;
        n_jumpy = n_jumpy+1;
    else
        prev_pos = pos;
    end
end

% ----------------------------------------------------------------------------------------------------------------------
