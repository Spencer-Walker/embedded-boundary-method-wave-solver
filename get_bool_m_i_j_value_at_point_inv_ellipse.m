function [ val ] = get_bool_m_i_j_value_at_point_inv_ellipse( point )
% returns the level-set value for a circle centered in the grid
    LSval = (point(1)-1.1)^2 + ((point(2)-1.1)/0.75)^2 - 1;
    
    if ( LSval < 0.0 )
        % inside
        val = 0;
    else
        % assuming outside (could be on the boundary)
        val = 1;
    end

end

