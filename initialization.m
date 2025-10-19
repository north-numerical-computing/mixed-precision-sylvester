format compact
format long
warning off

% External projects (git submodules).
addpath('deps/cpfloat/bin');
addpath('deps/anymatrix');
addpath('utils');

% Install and set up SLICOT group if not already available.
if all(~strcmp(anymatrix('g'), 'slicot'))
  anymatrix('g', 'slicot', 'mfasi/slicot');
  run([fileparts(which('anymatrix')) '/slicot/private/setup'])
  anymatrix('s');
end

if all(~strcmp(anymatrix('g'), 'sylvester_equations'))
  anymatrix('g', 'sylvester_equations', 'mfasi/sylvester_equations');
end

% Local folders.
addpath('methods')

datfolder = './datfiles/';
if ~exist(datfolder, 'dir')
  mkdir(datfolder)
end
