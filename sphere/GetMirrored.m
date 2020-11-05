
function [theta_mirrored, phi_mirrored] = GetMirrored(theta, phi)
    if(theta >= 0 && theta < pi/2)
        theta_mirrored = -theta + pi;
    elseif(theta >= pi/2 && theta < pi)
        d = pi-theta;
        theta_mirrored = theta - pi + d * 2;
    end
    if(phi <= pi && phi > 0)
        phi_mirrored = phi + pi;
    elseif(phi >= pi && phi < 2 * pi)
        phi_mirrored = phi - pi;
    end
end
