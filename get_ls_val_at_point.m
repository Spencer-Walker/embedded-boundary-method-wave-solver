function [val] = get_ls_val_at_point( point, topSplineData, botSplineData )

topPP = mkpp(topSplineData.topSpline.breaks, topSplineData.topSpline.coefs);
LSTop = ppval(topPP,point(1))-point(2);

botPP = mkpp(botSplineData.botSpline.breaks, botSplineData.botSpline.coefs);
LSBot = ppval(botPP,point(1))-point(2);

% point is near the boundary of the level set with smaller magnitude
if (abs(LSTop) < abs(LSBot))
    % near top boundary
    val = LSTop;
else
    % must be near bottom boundary
    val = LSBot;
end

end
