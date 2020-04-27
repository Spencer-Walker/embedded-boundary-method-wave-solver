function [ val ] = get_bool_m_i_j_value_at_point_circle( point )
% returns the level-set value for a circle centered in the grid
    LSval = ( point(1)-0.5 )^2 + ( point(2)-0.5 )^2 - 0.1^2;
    
    if     ( LSval < 0.0 )
        % inside
        val = 1;
    else
        % assuming outside (could be on the boundary)
        val = 0;
    end

end

