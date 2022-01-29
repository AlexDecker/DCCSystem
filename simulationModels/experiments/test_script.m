clc;
clear all;

name = 'wpt';
nt = 3;
nr = 3;
top_rect = []; %default

wpt = WPTSystem(name, nt, nr, top_rect);