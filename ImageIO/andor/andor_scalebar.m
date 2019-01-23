function [barlength_cm, fov_um] = andor_scalebar(image_in_cm, magnification, bar_in_um, varargin)
%function [barlength_cm, fov_um] = andor_scalebar(image_in_cm, magnification, bar_in_um, varargin)
%calculates the length of a scale bar given the size of the image in the presentation (poster etc) in cm, 
%the magnification used when recording the image, and the desired length of the scale bar in micrometers. 
%
%input: imagesize in cm, magnificatioin, desired barlength in µm
%output: barlength in cm, total FOV in µm
%
%Parameters are preset for the Andor Ixon 888 (ie. 1024 pixels, 16um pixel size),
%but this can be changed with 
%optional input:
%               nr_pixels
%               pixelsize (in micrometers on the chip)
%
%for scanned images where the pixel size in the sample is given you can use:
%magnfication = 1, nr_pixels (optional input) = Your nr of pixels, pixelsize (optional input) = pixel size given by the scanning software

params.nr_pixels = 1024; %Andor Ixon 888
params.pixelsize = 13; %um, Andor Ixon 888
params = omex_read_params(params, varargin);


fov_um= params.nr_pixels * params.pixelsize / magnification; %total FOV in um

barlength_cm = image_in_cm / fov_um * bar_in_um;


