clear
clc
close all
%%
load('one.txt')
load('two.txt')
load('three.txt')


figure 
hold on
plot([1 2 3],[mean(one) mean(two) mean(three)],'*r')
% plot(ones(25,1),one,"r")
% plot(ones(25,1)*2,two,'.r')
% plot(ones(25,1)*3,three,'.r')
set(gca, 'YScale', 'log')
