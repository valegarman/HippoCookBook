function swap_list = led_swap_filter(led_pos, led_pix)
% Checks for instances of two leds swapping or big one replacing little one
% when the big one gets obscured.
% Input xy posiiton of each led and
% and number of pixels in each. Big light is light number 1
% format: led_pos(1:n_pos, 1:num_cols, x-y), npix(1:n_pos, 1:num_cols)

thresh = 5;

mean_npix = nanmean(led_pix);
std_npix = nanstd(led_pix);

% Check if big light closer to small light at t-1 than to big light at t-1
% and small light either not found or closer to big light at t-1 than small light at t-1
pos = 2:size(led_pix,1);

% Use one of the two following blocks of code - the first calculates a city
% block metric and the second euclidian distance. For most applications the
% latter is correct.

% Calculate city block metric
% dist12 = abs(led_pos(pos, 1, 1)-led_pos(pos-1, 2, 1)) + abs(led_pos(pos, 1, 2)-led_pos(pos-1, 2, 2));
% dist11 = abs(led_pos(pos, 1, 1)-led_pos(pos-1, 1, 1)) + abs(led_pos(pos, 1, 2)-led_pos(pos-1, 1, 2));
% dist21 = abs(led_pos(pos, 2, 1)-led_pos(pos-1, 1, 1)) + abs(led_pos(pos, 2, 2)-led_pos(pos-1, 1, 2));
% dist22 = abs(led_pos(pos, 2, 1)-led_pos(pos-1, 2, 1)) + abs(led_pos(pos, 2, 2)-led_pos(pos-1, 2, 2));

%Calculate eucldian
dist12 = sqrt(nansum(((squeeze(led_pos(pos,1,:))-squeeze(led_pos(pos-1,2,:))).^2),2));
dist11 = sqrt(nansum(((squeeze(led_pos(pos,1,:))-squeeze(led_pos(pos-1,1,:))).^2),2));
dist21 = sqrt(nansum(((squeeze(led_pos(pos,2,:))-squeeze(led_pos(pos-1,1,:))).^2),2));
dist22 = sqrt(nansum(((squeeze(led_pos(pos,2,:))-squeeze(led_pos(pos-1,2,:))).^2),2));

switched = (dist12 < dist11-thresh) & ( isnan(led_pos(pos, 2, 1)) |(dist21 < dist22-thresh) );

% Check if size of big light has shrunk to be closer to that of small light (as Z score)
z11 = (mean_npix(1) - led_pix(pos, 1))/std_npix(1);
z12 = (led_pix(pos, 1) - mean_npix(2))/std_npix(2);
shrunk = z11 > z12;
swap_list = find( switched & shrunk ) + 1;

% ----------------------------------------------------------------------------------------------------------------------
