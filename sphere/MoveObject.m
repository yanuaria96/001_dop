
function obj = MoveObject(obj)
    v = obj.velocity;
    obj.phi = obj.phi + v.phi;
    obj.theta = obj.theta + v.theta;
end
