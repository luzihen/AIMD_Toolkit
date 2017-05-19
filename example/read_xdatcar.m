%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIT License
%
% Copyright (c) 2017 Ziheng LU
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,lattice] = read_xdatcar(filename_in,fframes_in)
% this function reads in the XDATCAR of a vasp AIMD simulations

% Input: fileanme_in - the name of your vasp XDATCAR in a string format.
% Usually it is 'XDATCAR'. This version does not support varible cell MD.
% Pls carefully check the XDATCAR format before using it.

% Output: data - a 3-dimentional matrix(N,M,t). N is the number of
% particles. M = 1 is the atom type. M = 2:4 are the x,y,z fractional
% coordinates. t is the simulation steps. lattice - a 3x3 matrix
% of the lattice vector.


filename=filename_in;
fframes=fframes_in;


fid = fopen(filename);

fprintf('start reading XDATCAR file\n');
% line 1 system_name
sys_name = fgetl(fid);
% line 2 lattice_mutiplier
lattice_mutiplier = fgetl(fid);
% lattice vector 1
tline = fgetl(fid);
lattice_a=strread(tline);
% lattice vector 2
tline = fgetl(fid);
lattice_b=strread(tline);
% lattice vector 3
tline = fgetl(fid);
lattice_c=strread(tline);
% number of atoms
tline = fgetl(fid);
tline = fgetl(fid);
num_atoms=strread(tline);
if num_atoms(end)==0
    num_atoms(end)=[];
end
tot_num=sum(num_atoms);

%creat an array for coordinate
data=[];



frames=1;
while frames<=fframes
    data_temp=[];
    %read corredinate
    tline = fgetl(fid);%line 'Direct'
    for num_type=1:size(num_atoms,2)%read type by type
        for atom=1:num_atoms(num_type)
            tline = fgetl(fid);
            coordinate=strread(tline);
            %if num_atoms(end)==0
            if size(coordinate,2)==4
                coordinate(4)=[];
            end
            coordinate=[num_type coordinate];
            data_temp=[data_temp;coordinate];
        end
    end
    data(:,:,frames)=data_temp;
%     for ii=1:6
%         tline = fgetl(fid);
%     end
    frames=frames+1;
end
fclose(fid);

lattice=[lattice_a;lattice_b;lattice_c];
end