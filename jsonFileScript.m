clear; close all; clc;
% load spline data 
topSplineData=load('splines/topSpline_centerUFO.mat');
botSplineData=load('splines/botSpline_centerUFO.mat');

% upload a .json file and chekc its contents
fname = 'input.json';
val = jsondecode(fileread(fname));

% create grid vectors from file data
xVec = 0:(val.L/(val.M-1)):val.L;
yVec = xVec;

% loop over coords and determine m_i_j vals
nx=length(xVec);
ny=length(yVec);

mij=zeros(nx,ny);

alpha_1=zeros(nx,ny);    % horizontal distance to boundary for boundary points ( mij = -1 )
alpha_2=zeros(nx,ny);    % vertical distance to boundary for boundary points ( mij = -1 )

for i=1:nx
    for j=1:ny        
        tmpPoint=[xVec(i) yVec(j)];
        tempBool=get_bool_m_i_j_value_at_point( tmpPoint, topSplineData, botSplineData );
        
        if ( tempBool==0 )    % point is outside            
            mij(i,j)=0;
        else                % point is inside, need to check neighbors
            if     ( i==1 || i==nx )    % boundary point, default to outside
                mij(i,j)=0;
            elseif ( j==1 || j==ny )    % boundary point, default to outside
                mij(i,j)=0;
            else                
                leftPoint=[xVec(i-1) yVec(j)];
                lftPtBool=get_bool_m_i_j_value_at_point( leftPoint, topSplineData, botSplineData );
                
                rightPoint=[xVec(i+1) yVec(j)];
                rhtPtBool=get_bool_m_i_j_value_at_point( rightPoint, topSplineData, botSplineData );
                
                topPoint=[xVec(i) yVec(j+1)];
                tpPtBool=get_bool_m_i_j_value_at_point( topPoint, topSplineData, botSplineData );
                
                botPoint=[xVec(i) yVec(j-1)];
                btPtBool=get_bool_m_i_j_value_at_point( botPoint, topSplineData, botSplineData );
                
                if( lftPtBool==1 && rhtPtBool==1 && tpPtBool==1 && btPtBool==1 )
                    mij(i,j)=1;
                else
                    mij(i,j)=-1;
                    % ----- calculate alpha1 and alpha2 for this point ----- % 
                    
                    % alpha_1 calculation
                    if ( lftPtBool==0 )     % left point is outside
                        % interpolate to find the distance to the boundary on the left                        
                        LS_left  = get_ls_val_at_point( leftPoint, topSplineData, botSplineData );
                        LS_right = get_ls_val_at_point( tmpPoint, topSplineData, botSplineData );   % this is the current grid point                        
                        LS_bound = 0; % LS = 0 on boundary
                        
                        x_bound = leftPoint(1) + (LS_bound - LS_left)/(LS_right - LS_left)*(tmpPoint(1)-leftPoint(1));
                        
                        alpha_1(i,j) = tmpPoint(1) - x_bound;
                        
                    elseif ( rhtPtBool==0 ) % right point is outside
                        % interpolate to find the distance to the boundary on the right                                              
                        LS_left  = get_ls_val_at_point( tmpPoint, topSplineData, botSplineData );   % this is the current grid point                        
                        LS_right = get_ls_val_at_point( rightPoint, topSplineData, botSplineData );   
                        LS_bound = 0; % LS = 0 on boundary
                        
                        x_bound = leftPoint(1) + (LS_bound - LS_left)/(LS_right - LS_left)*(rightPoint(1)-tmpPoint(1));
                        
                        alpha_1(i,j) = x_bound - tmpPoint(1);
                        
                    else                    % left and right points are either on the boundary or inside domain
                        alpha_1(i,j)=0;
                    end                    
                    
                    % alpha_2 calculation
                    if ( tpPtBool==0 )     % top point is outside         
                        % interpolate to find the distance to the boundary on the top                                              
                        LS_top = get_ls_val_at_point( topPoint, topSplineData, botSplineData );
                        LS_bot = get_ls_val_at_point( tmpPoint, topSplineData, botSplineData );   % this is the current grid point                        

                        LS_bound = 0; % LS = 0 on boundary
                        
                        y_bound = topPoint(2) + (LS_bound - LS_top)/(LS_bot - LS_top)*(tmpPoint(2)-topPoint(2));
                        
                        alpha_2(i,j) = y_bound - tmpPoint(2);
                        
                    elseif ( btPtBool==0 ) % right point is outside
                        % interpolate to find the distance to the boundary on the bottom                                              
                        LS_top = get_ls_val_at_point( tmpPoint, topSplineData, botSplineData ); % this is the current grid point                        
                        LS_bot = get_ls_val_at_point( botPoint, topSplineData, botSplineData );   

                        LS_bound = 0; % LS = 0 on boundary
                        
                        y_bound = tmpPoint(2) + (LS_bound - LS_top)/(LS_bot - LS_top)*(botPoint(2)-tmpPoint(2));
                        
                        alpha_2(i,j) = tmpPoint(2) - y_bound;
                        
                    else                    % top and bottom points are either on the boundary or inside domain
                        alpha_2(i,j)=0;
                    end
                end                
            end            
        end        
    end
end

imagesc(mij);
colormap(gray);
colorbar
%% save results
disp(size(alpha_1))
alpha_1 = reshape(alpha_1, val.M^2, 1);
alpha_2 = reshape(alpha_2, val.M^2, 1);
mij = reshape(mij, val.M^2, 1);
save('mij.out','mij',"-ascii");
save('alpha_1.out','alpha_1',"-ascii");
save('alpha_2.out','alpha_2',"-ascii");

