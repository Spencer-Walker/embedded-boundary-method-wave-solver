function [val] = get_bool_m_i_j_value_at_point( point, topSplineData, botSplineData )
% returns the 'pseudo' m_i_j value as it does not check neighboring point
%
% i.e. val = 0 if point outside 
%          = 1 if point inside

topPP = mkpp(topSplineData.topSpline.breaks, topSplineData.topSpline.coefs);
LSTop = ppval(topPP,point(1))-point(2);

botPP = mkpp(botSplineData.botSpline.breaks, botSplineData.botSpline.coefs);
LSBot = ppval(botPP,point(1))-point(2);

if     ( LSTop > 0.0 && LSBot < 0.0 )
    % inside 
    val = 1;
else
    % assuming outside (could be on the boundary)
    val = 0;
end

end

