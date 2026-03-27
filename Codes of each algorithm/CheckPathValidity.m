function isValid = CheckPathValidity(pos, model)
    n = model.n;
    waypoints = reshape(pos, [3, n])';

    isValid = true;
    margin = model.Safeh;

    for i = 1:n
        x = waypoints(i, 1);
        y = waypoints(i, 2);
        z = waypoints(i, 3);

        % 海拔高度 + 安全高度限制
        terrain_z = interp2(model.x_data, model.y_data, model.z_data, x, y, 'linear', Inf);
        if z < terrain_z + margin
            isValid = false;
            return;
        end

        % 检查是否碰撞障碍物（以圆柱形近似）
        for b = 1:model.Num_Barrier
            bx = model.Barrier(b,1);
            by = model.Barrier(b,2);
            br = model.Barrier(b,3);

            if norm([x - bx, y - by]) < br && z <= model.zmax
                isValid = false;
                return;
            end
        end
    end
end
