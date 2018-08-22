

clear all
close all
nt = 801;
time = (0:nt-1)*0.004;
fid = fopen('in/sv1','r');
data = fread(fid,'single');
sv = reshape(data,161,2,801);

fid = fopen('out/tinti','r');
data = fread(fid,'single');
svf = reshape(data,161,2,801);

ix = 61;

figure(1)
plot(time,squeeze(sv(ix,1,:)));
hold on
plot(time,squeeze(svf(ix,1,:)));


