function [ data_whole,data_out,loc ] = cutfile( filename )
%CUTFILE Roughly cut the audio to phonemes based on energy
%   [x,y,z] = CUTFILE(FILENAME) filename is the name of audio file, in string format
%   y is a 1xnum(peaks)+1 cell, containing phonemes divided by energy, x is
%   the original data as a whole. z is the location of through of wave.
%   Author:Zhonghao Guo
%   Date:2018.3
[data,fs] = wavread(filename);
data = resample(data,10000,fs);
data_o = data;
data_whole = data;
data = envelope(data);
data = smooth(data,2500,'lowess');
[~,loc] = findpeaks(-data);
data_out = cell(1,length(loc)+1);
data_out{1,1} = (data_o(1:loc(1)));
for i=2:length(loc)
    data_out{1,i} = data_o((1+loc(i-1)):loc(i));
end
data_out{1,i+1} = data_o((1+loc(i)):length(data));
end

