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
%
% To use the code, please copy the POSCAR and XDATCAR file
% to matlab working directory. Set the user defined parameters and run the code.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

%user defined parameters
trajname='XDATCAR'; %trajectory name in lammpstrj format
N_step=2000; %number of simulation steps
N=30;                     %define N grid
Nparticle=1:8;          %define particles to be included in the density calculation
smoothing_tag=2;          %1-call smoothing, 2-no smoothing
r=0.8; %radius






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%starting reading data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[data,lattice]=read_poscar(trajname,N_step);
fprintf(' done reading data\n start calculating density...\n');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%starting reading data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spacing=1.0/(N-1);
Nsteps=N_step;%size(data);
density=zeros(N,N,N);

%Doing the counting
for s_index=1:Nsteps
    for p_index=Nparticle
        for x_ii=1:N
            for y_ii=1:N
                for z_ii=1:N
                    x_coor=rem(data(p_index,2,s_index)+999,1);
                    y_coor=rem(data(p_index,3,s_index)+999,1);
                    z_coor=rem(data(p_index,4,s_index)+999,1);
                    vector=[((x_ii)/N-x_coor) ((y_ii)/N-y_coor) ((z_ii)/N-z_coor)];
                    
                    if norm(lattice*vector')<=r
                    density(x_ii,y_ii,z_ii)=density(x_ii,y_ii,z_ii)+1;
                    end
                end
            end
        end
    end
    s_index
end
density=density/N_step;
fprintf(' done calculating raw density (normalized by simulation steps)\n');

temp=density(:,:,1)+density(:,:,N);
density(:,:,1)=temp;
density(:,:,N)=temp;
temp=density(:,1,:)+density(:,N,:);
density(:,1,:)=temp;
density(:,N,:)=temp;
temp=density(1,:,:)+density(N,:,:);
density(1,:,:)=temp;
density(N,:,:)=temp;


if smoothing_tag==1
fprintf(' smoothing called...\n');
density_new=zeros(N+1,N+1,N+1);
	for x_i=2:N
		for y_i=2:N
			for z_i=2:N
				density_new(x_i,y_i,z_i)=density(x_i,y_i,z_i)...
				+density(x_i+1,y_i,z_i)+density(x_i-1,y_i,z_i)...
				+density(x_i,y_i+1,z_i)+density(x_i,y_i-1,z_i)...
				+density(x_i,y_i,z_i+1)+density(x_i,y_i,z_i-1);
			end
		end
	end
	density_new=density_new./7;
	for x_i=1:N+1
		for y_i=[1 N+1]
			for z_i=[1 N+1]
				density_new(x_i,y_i,z_i)=density(x_i,y_i,z_i);
			end
		end
	end
	for x_i=[1 N+1]
		for y_i=1:N+1
			for z_i=[1 N+1]
				density_new(x_i,y_i,z_i)=density(x_i,y_i,z_i);
			end
		end
	end
	for x_i=[1 N+1]
		for y_i=[1 N+1]
			for z_i=1:N+1
				density_new(x_i,y_i,z_i)=density(x_i,y_i,z_i);
			end
		end
    end
    density=density_new;
	fprintf(' smoothing finished\n writing data to CHG file...\n');
	else
	fprintf(' smoothing not called\n writing data to CHG file...\n');
end
	



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% writing to CHG file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%append the POSCAR to CHGCAR first
system('del CHG');
system('type POSCAR>>CHG');

%open file and write the grid number
fileID = fopen('CHG','a');
fprintf(fileID,'\n\r');
%fprintf(fileID,'   ',size(density));
fprintf(fileID,'%d   ',size(density));
fprintf(fileID,'\n');

% four data in a line
count=0;
for k=1:N
    for j=1:N
        for i=1:N
            count=count+1;
            if rem(count,4)==0
            fprintf(fileID,' %E\n',density(i,j,k));
            else
            fprintf(fileID,' %E',density(i,j,k));    
            end
        end
    end
end
fclose(fileID);
