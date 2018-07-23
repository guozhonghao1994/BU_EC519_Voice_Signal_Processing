% Main function for proj2
% Recommended frame #:
% Whole Figure: 1:90
% Vowel Window: 31:32
% Unvoiced Window: 66:67
% VFV Window: 50:54
p_opt;
figure()
start_frame = 1; % start_frame #
end_frame = 90; % end_frame #
mode = 'w'; % 'w' for whole, 'v' for vowel, 'u' for unvoiced, 'c' for v+f+v combination
[SNR_e,SNR_u] = LP(start_frame,end_frame,mode) % SNR is shown on command window