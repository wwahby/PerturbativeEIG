function vec_flipped = flip_vector(vec_toflip, vec_ref)
% Figure out which way to orient vec_toflip based on how closely it matches
% vec_ref in either orientation

vec_diff = sum( (vec_toflip - vec_ref).^2);
vec_flipped = fliplr(vec_toflip);
vec_diff_2 = sum( (vec_flipped - vec_ref).^2);

if(vec_diff_2 > vec_diff)
    vec_flipped = vec_toflip;
end