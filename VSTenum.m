classdef VSTenum < double
    % VSTree enumerations:
    %   Vnode           - Octree node denoted as -1
    %   Tnode           - Transition node denoted as -2
    %   Snode           - Quadtree node denoted as -3
    %   SnodeL          - Leaf node of Quadtree denoted as -4
    
    enumeration
        Vnode (-1)
        Tnode (-2)
        Snode (-3)
        SnodeL(-4)
        Empty (-5)
    end
end

