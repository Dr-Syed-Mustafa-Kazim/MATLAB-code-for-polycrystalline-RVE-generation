% # MATLAB-code-for-polycrystalline-RVE-generation # 

%% Author: Dr. Syed Mustafa Kazim (Ph.D., IIT Kanpur) %%

% A MATLAB code is developed as a part of my Ph.D. work to implement a 
% novel algorithm for generating a 3D Voronoi tessellated polycrystalline 
% RVE microstructure for Ti-alloys in which the primary-$\alpha$ grains are
% made to nucleate at the triple points and grain boundaries of the 
% transformed-$\beta$ colonies. The generated microstructure retains the 
% experimentally observed ratio of sizes and volume fraction of the 
% primary-$\alpha$ grains and transformed-$\beta$ colonies. The algorithm 
% assigns Euler angles and also distributes the two phases randomly to the 
% grains in the RVE. It works by importing the mesh details from the ABAQUS 
% input file and then forming the two kinds of grain clusters out of the
% RVE elements. This code is generic and can be extended to generate a 
% similar 3D Voronoi tessellated polycrystalline RVE microstructure for 
% any other material system. 
% 
% The novel MATLAB code is uploaded on GitHub on the following path and 
% can be availed without any charges:  
% 
% \textbf{https://github.com/Dr-Syed-Mustafa-Kazim/MATLAB-code-for-polycrystalline-RVE-generation
% 

%%%%%%%%%%%% SMK CODE FOR 3D VORONAI GENERATION %%%%%%%%%%%%
% \begin{verbatim}

% ##### CODE FOR 3D VORONAI GENERATION #####
clc
clear
close all

%% SMK CODE 
%% LOGIC 1
no_of_grains =30 ;  %20   %6; %good 10 , 12
volume=1*1*1; % micrometer
rad = ((0.75*volume)/no_of_grains/pi)^(1/3);  % SMK: L^3=1
dia = 1.7*rad;
max_iter=10000000000000000000;
x = zeros(no_of_grains,1);
y=zeros(no_of_grains,1);
z=zeros(no_of_grains,1);
seed=1;
rng(seed)
x(1)=rand(1);
y(1)=rand(1);
z(1)=rand(1);
all_prev = zeros(no_of_grains,3);
iter=zeros(no_of_grains,1);
gap=zeros(no_of_grains,1);

for nel=2:no_of_grains
    all_prev(nel-1,:) = [x(nel-1), y(nel-1), z(nel-1)];
    d=dia/2;
    i=0;
    while d<dia && i<max_iter
        i=i+1;
        x(nel)=rand(1);
        y(nel)=rand(1);
        z(nel)=rand(1);
        current = [x(nel),y(nel),z(nel)];
        k = dsearchn(all_prev,current);
        
        d = norm(current-all_prev(k,:));
    end
    iter(nel)=i;
    gap(nel)=d;
end
%plot(iter)
plot(gap)
hold on
plot([1:no_of_grains],dia*ones(no_of_grains))
hold off
all_prev(nel,:)=[x(1),y(1),z(1)];
scatter3(x,y,z);
cent_vor_x = x;
cent_vor_y = y;
cent_vor_z = z;
clear x y z


count=0;
%% Loading the input files and reading Nodes and element set
% inpFileName1 = 'struc_mesh_voronoi_2_cpe4r.inp';
inpFileName1 = '2197e_CMNAME_RF_U.inp';
[node, element1, elementType1] = readinp(inpFileName1);
elem_needed = [1:length(element1)];
%% Finding the centroids of C3D8R elements
for elem=1:length(element1)
    if ismember(element1(elem,1),elem_needed)
        count=count+1;
        
        node_set_1(count,:)=element1(elem,2:end);
        for i=1:size(element1,2)-1
            x(i)=node(node_set_1(count,i),2);
            y(i)=node(node_set_1(count,i),3);
            z(i)=node(node_set_1(count,i),4);
        end
        %polyin=polyshape(x,y,z);
% %         [cent_elem_x(count),cent_elem_y(count),cent_elem_z(count)]=centroid(polyin);        
%polyin=polyshape(x,y,z);    

% SMK: centriod calculation
cent_elem_x(count)=mean(x);
cent_elem_y(count)=mean(y);
cent_elem_z(count)=mean(z);
        
    end
end

clear x y z k
P = [cent_vor_x(:,1),cent_vor_y(:,1),cent_vor_z(:,1)]; % Centroid set of voronoi
cent_elem_x=cent_elem_x';
cent_elem_y=cent_elem_y';
cent_elem_z=cent_elem_z';
Q = [cent_elem_x(:,1),cent_elem_y(:,1),cent_elem_z(:,1)]; % Centroid set of Elements
offset=[0,0,0];
P(:,:,:)=P(:,:,:)+offset; %% Offsetting the Voronoi set to match the abaqus coordinates

k = dsearchn(P,Q); %% k contains the closest voronoi centroid of elements

%% Creating random euler angles = (no_of_grains)
for i=1:no_of_grains
    euler1(i)=randi([0 359]);
    euler2(i)=randi([0 89]);
    euler3(i)=randi([0 359]);
    %     phase(i)=randi([1 2]);
end  

%%%%%%%%
no_of_grains_alpha=round(0.15*no_of_grains)  % SMK: PERCENTAGE OF ALPHA GRAINS =15% as per serrated paper


for i=1:no_of_grains_alpha    
    phase_alpha(i)=1;
end  

for i=1:(no_of_grains-no_of_grains_alpha)    
    phase_beta(i)=2;  %randi([1 2]);
end  

phase1=[phase_alpha(randperm(length(phase_alpha))) phase_beta(randperm(length(phase_beta)))]

phase2=phase1(randperm(length(phase1)))

phase=phase2(randperm(length(phase2)))

%%%%%%%%%%

%% Writing euler angles to the respective voronoi cells
fileID = fopen('texture_2197.dat','w');
     

    fprintf(fileID,'       xeul(1,%d)=%.2fd0\n',elem,euler1(k(elem)));
    fprintf(fileID,'       xeul(2,%d)=%.2fd0\n',elem,euler2(k(elem)));
    fprintf(fileID,'       xeul(3,%d)=%.2fd0\n',elem,euler3(k(elem)));    
    
%   fprintf(fileID,'%d,%d,0,0,%d \n',elem,euler(k(elem)),r(elem));
    
end


fclose(fileID);


%% SMK: Grouping elements of same phase(MAT1 or MAT2) for INPUT FILE


file_phase = fopen('elem_set_alpha_beta_CMNAME.txt','w');


iread_alpha=0;
iread_beta=0;
alpha_count=1;
beta_count=1;

for elem=1:length(k)
    fprintf('For element %d, PHASE is= %d, VORONAI is = %d \n',elem,phase(k(elem)),k(elem))
    
    if (phase(k(elem))==1)
        
    
        alpha_elem(alpha_count)=elem   ;     
        alpha_count=alpha_count+1   ; 
        
    elseif(phase(k(elem))==2)
        
        beta_elem(beta_count)=elem ;       
        beta_count=beta_count+1  ;  

    end
    
end



for i=1:length(alpha_elem)       
       
        if(iread_alpha==0)
            fprintf(file_phase,'\n*Elset, elset=ALPHA\n');
            iread_alpha=1;
        end        
     
                if(mod(i,10)==0)            
                    fprintf(file_phase,'%d\n',alpha_elem(i));            
                else            
                    fprintf(file_phase,'%d,',alpha_elem(i));            
                end
end


for i=1:length(beta_elem)       
       
        if(iread_beta==0)
            fprintf(file_phase,'\n*Elset, elset=BETA\n');
            iread_beta=1;
        end        
     
                if(mod(i,10)==0)            
                    fprintf(file_phase,'%d\n',beta_elem(i));            
                else            
                    fprintf(file_phase,'%d,',beta_elem(i));            
                end
end

fclose(file_phase);
phase
length(phase);
percentage_alpha=(length(alpha_elem)/((length(alpha_elem)+length(beta_elem))))*100


% \end{verbatim}
%%%

    


% \end{appendices}
