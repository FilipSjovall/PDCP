% Postprocessing
% import results from Fortran program
clc; clear; clf; close all;
load data.mat

ed=extract(edof,a);
nnod = length(coord);
for i = 1:nnod
   dof(i,:) = [2*i-1 2*i]; 
end
edof = double(edof);
[ex,ey]=coordxtr(edof,coord,dof,6);

nelem = length(edof);

figure(15)
    clf
    for el = 1:nelem
%       indx = mesh.enod(el,2:end);
      ix   = edof(el,2:2:6);
      iy   = edof(el,3:2:7);
      mex(el,:) = ex(el,1:3) + a(ix)';
      mey(el,:) = ey(el,1:3) + a(iy)';
    end
    
    patch(mex',mey',1000*ones(nelem,3)')