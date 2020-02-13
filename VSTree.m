classdef VSTree < handle
    %VSTree implementation for 3D Point subdivision. 
    %   자세한 설명 위치
    % Property values : 
    %   ptsnormal       - Normal vector of pts. If not given, normal is
    %                   estimated for each pts from PCA method
    %   nodeCapacity    - Maximum # of points a node may contain 
    %                   Defaults to ceil(numPts/10)
    %   maxDepth        - Maximum # of times a node may be subdivided. 
    %                   Defaults to (INF, INF) - V, S node each
    %   minSize         - Minimum size of a cell, this prohibit infinite
    %                   division due to duplicate points
    %   errorFun        - Stopping criteria for V-node construction
    %                   Defaults to L^2 metric
    %   kappaFun        - Transition crteria for T-node construction
    %                   Defaults to height field indicator from the 
    %                   original paper
    %   delta_a         - Angular threshold for height field indicator
    %   delta_d         - Plane distance threshold for height field
    %                   indicator
    %   errorTh         - Threshold of error metric (Default 1e-3)
    %
    % VSTree methods(T.B.D.):    
    %   query           - Ask which nodes a new set of points belongs to
    %   plot, plot3     - Plots node bounding boxes to the current axes
    % VSTree properties:
    %   Points          - The coordinate of points in the subdivision
    %   PointNodes      - Indices of the leaf node that each point belongs
    %   NodeCount       - Total # of nodes created
    %   NodeBoundaries  - NodeCount-by-6 [MIN MAX] coordinates of node
    %                   edges 
    %   NodeDepth       - The # of subdivisions to reach. 2nd row
    %                   represents type ((-1) for V, (-2) for T, (-3) for
    %                   S, (-4) for leaf node of S)
    %   NodeParents     - Indices of the node that each node belongs to
    %   NodeFrame       - Parameter of local frame (plane) for each node,
    %                   [ui', vi', ni' | Center(1x3)] 
    %   Properties      - Name/Val pairs used for creation
    % See also, 
    %   "2006_EUROGRAPHICS_T. Boubekuer_Volume-Surface Trees"

    properties
        Points;
        PointNodes;
        NodeCount;
        NodeBoundaries;
        NodeDepths;
        NodeParents = zeros(0,1);
        NodeFrames;
        Properties;
    end    
    
    methods
        function obj = VSTree(pts,varargin)
            % VSTree constructor
            validateattributes(pts, {'numeric'},...
                {'real','finite','nonnan','ncols',3},...
                mfilename,'PTS');
            % Initialize a single node surrounding all given points
            numPts = size(pts,1);
            obj.NodeBoundaries = [min(pts,[],1), max(pts,[],1)];
            obj.Points = pts;
            obj.PointNodes = ones(numPts,1);
            obj.NodeDepths = [0; double(VSTenum.Vnode)];
            obj.NodeParents(1) = 0;
            obj.NodeCount = 1;
            
            % Allow custom setting of Properties
            IP = inputParser;
            IP.addParameter('nodeCapacity', ceil(numPts)/10);
            IP.addParameter('maxDepth', [inf, inf]);
            IP.addParameter('minSize', 1000*eps);
            IP.addParameter('errorFun', @Geometric);
            IP.addParameter('kappaFun', @HfieldIndi);
            IP.addParameter('delta_a', 0);
            IP.addParameter('delta_d', 1/6);
            IP.addParameter('errorTh', 1e-3);
            % If point normals are not provided, these are automatically
            % estimated by using built-in function pcnormals
            % Note. pcnormals should be altered to a function which solves
            % sign ambiguity by using smoothness (See H. Hoppe`s paper)
            assert(numPts>9, 'Error : Too few points, more than 10 points needed');
            IP.addParameter('ptsnormal', pcnormals(pointCloud(obj.Points), 10));
            IP.parse(varargin{:});
            obj.Properties = IP.Results;
                                    
            % Estimate local frame of the first node
            obj.NodeFrames = computeLocalFrame(obj, 1);
            
            % Start dividing
            obj.preallocateSpace;
            obj.divide(1);
            obj.deallocateSpace;
        end
        
        function preallocateSpace(obj)
            numPts = size(obj.Points,1);
            numNodes = numPts;
            if isfinite(obj.Properties.nodeCapacity)
                numNodes = ceil(2*numPts/obj.Properties.nodeCapacity);
            end
            obj.NodeDepths(1, numNodes) = 0;
            obj.NodeParents(numNodes) = 0;
            obj.NodeBoundaries(numNodes,1) = 0;
            obj.NodeFrames(numNodes,1) = 0;
        end
        
        function deallocateSpace(obj)
            obj.NodeDepths(:, obj.NodeCount+1:end) = [];
            obj.NodeParents(obj.NodeCount+1:end) = [];
            obj.NodeBoundaries(obj.NodeCount+1:end,:) = [];
            obj.NodeFrames(obj.NodeCount+1:end,:) = [];
        end
        
        function divide(obj, startingNodes)
            %Loop over each node we will consider for division
            for i = 1:length(startingNodes)
                nodeNo = startingNodes(i);
                % If parent node is not T or S node
                if ((obj.NodeParents(nodeNo) == 0) || ...
                        (obj.NodeDepths(2, nodeNo) == VSTenum.Vnode))
                    % Do Octree partitioning
                    oldCount = obj.NodeCount;
                    obj.oct_divideNode(nodeNo);
                    obj.divide(oldCount+1:obj.NodeCount);                    
                elseif (obj.NodeDepths(2, nodeNo) ~= VSTenum.Empty) && ...
                        (obj.NodeDepths(2, nodeNo) ~= VSTenum.SnodeL)
                    % Do Quadtree partitioning
                    % Prevent dividing beyond the maximum depth
                    % Note. T-node has 0 value over node depths                    
                    
                    % Check whether it passes Error Criterion
                    if obj.Properties.errorFun(obj, nodeNo)
                        obj.NodeDepths(2, nodeNo) = VSTenum.SnodeL;
                        continue;
                    end

                    oldCount = obj.NodeCount;                    
                    obj.quad_divideNode(nodeNo);
                    obj.divide(oldCount+1:obj.NodeCount);                    
                else
                    % Do nothing
                end
            end
        end
        
        function oct_divideNode(obj, nodeNo)
            % Gather the new points
            nodePtMask = obj.PointNodes==nodeNo;
            objNodePoints = obj.Points(nodePtMask,:);
            
            % Get the old corner points and the new division point
            oldMin = obj.NodeBoundaries(nodeNo,1:3);
            oldMax = obj.NodeBoundaries(nodeNo,4:6);
            newDiv = mean([oldMin; oldMax], 1);
            % Build the new boundaries of 8 subdivisions
            % newBounds = [BLB, TLB, BRB, TRB, BLF, TLF, BRF, TRF]
            % Bottom/Top, Left/Right, Front/Back
            minMidMax = [oldMin, newDiv, oldMax];
            newBounds = minMidMax([...
                1 2 3 4 5 6;
                1 2 6 4 5 9;
                1 5 3 4 8 6;
                1 5 6 4 8 9;
                4 2 3 7 5 6;
                4 2 6 7 5 9;
                4 5 3 7 8 6;
                4 5 6 7 8 9]);
            
            % Determine to which of these 8 nodes each current point belongs
            nodeMap = cat(3,[0 0 0],[0 0 1],[0 1 0],[0 1 1],...
                [1 0 0],[1 0 1],[1 1 0],[1 1 1]);
            gtMask = bsxfun(@gt, objNodePoints, newDiv);
            [~,nodeAssignment] = max(all(bsxfun(@eq,gtMask,nodeMap),2),[],3);            
            % [~, nodeAssignment] = ismember(gtMask,nodeMap,'rows'); % A little slower than above.
            
            % Make the new nodes and reassign old points to them
            newNodeInds = obj.NodeCount+1:obj.NodeCount+8;
            obj.NodeBoundaries(newNodeInds,:) = newBounds;
            obj.PointNodes(nodePtMask) = newNodeInds(nodeAssignment);
            obj.NodeParents(newNodeInds) = nodeNo;
            obj.NodeCount = obj.NodeCount + 8;
            
            % Test each divided node whether it can be turned into T-node
            for cidx = 1:8
                cellNo = newNodeInds(cidx);
                cellPoints = obj.Points(obj.PointNodes == cellNo,:);
                obj.NodeDepths(1, cellNo) = obj.NodeDepths(1, nodeNo)+1;
                if isempty(cellPoints)                    
                    obj.NodeDepths(2, cellNo) = double(VSTenum.Empty);                
                    continue;
                end
                cellNormals = obj.Properties.ptsnormal(obj.PointNodes == cellNo, :);                
                objNodeFrame = computeLocalFrame(obj, cellNo);                
                obj.NodeFrames(cellNo,:) = objNodeFrame;
                obj.NodeDepths(2, cellNo) = double(VSTenum.Vnode);     
                % Kappa test, delta_a, delta_d are set to 0, 1/6. (can be
                % modified w.r.t. input data)                
                flag = obj.Properties.kappaFun(obj, cellPoints, cellNormals, obj.NodeFrames(cellNo,:));
                if flag || checkNodeValidity(obj, cellNo)                    
                    obj.NodeDepths(2, cellNo) = double(VSTenum.Tnode);
                    % Project 3D points onto the local frame
                    rotMat = reshape(objNodeFrame(1:9), 3, 3);
                    objNodePointsProj = (cellPoints - ...
                        repmat(objNodeFrame(10:12), size(cellPoints,1), 1))*rotMat;
                    objNodePointsProj(:,3) = 0; %Orthographic projection
                    projMin = min(objNodePointsProj, [], 1);
                    projMax = max(objNodePointsProj, [], 1);
                    obj.NodeBoundaries(cellNo, :) = [projMin, projMax];                                    
                end
            end
        end
        
        function quad_divideNode(obj, nodeNo)
            % Gather the new points
            nodePtMask = obj.PointNodes==nodeNo;
            objNodePoints = obj.Points(nodePtMask,:);            
            objNodeFrame = obj.NodeFrames(nodeNo,:);
            
            % Project 3D points onto the local frame
            rotMat = reshape(objNodeFrame(1:9), 3, 3);
            objNodePointsProj = (objNodePoints - ...
                repmat(objNodeFrame(10:12), size(objNodePoints,1), 1))*rotMat;
            objNodePointsProj(:,3) = 0; %Orthographic projection
            
            % Compute plane boundary 
            projMin = obj.NodeBoundaries(nodeNo,1:3);
            projMax = obj.NodeBoundaries(nodeNo,4:6);
            projDiv = mean([projMin; projMax], 1);
            
            % Build the new boundaries of our 4 subdivisions
            % newBounds = [TL, BL, TR, BR]
            % Bottom/Top, Left/Right
            minMidMax = [projMin projDiv projMax];
            newBounds = minMidMax([...
                1 2 3 4 5 6;
                4 2 6 7 5 9;
                1 5 3 4 8 6;
                4 5 6 7 8 9]);
            
            % Determine to which of these 4 nodes each current point
            % belongs
            nodeMap = cat(3, [0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]);
            gtMask = bsxfun(@gt, objNodePointsProj, projDiv);
            [~, nodeAssignment] = max(all(bsxfun(@eq, gtMask, nodeMap),2),[],3);
            
            newNodeInds = obj.NodeCount+1:obj.NodeCount+4;
            obj.NodeBoundaries(newNodeInds,:) = newBounds;
            obj.PointNodes(nodePtMask) = newNodeInds(nodeAssignment);
            obj.NodeParents(newNodeInds) = nodeNo;
            obj.NodeCount = obj.NodeCount+4;
            obj.NodeDepths(1, newNodeInds) = obj.NodeDepths(1, nodeNo)+1;
            for cidx = 1:4
                cellNo = newNodeInds(cidx);
                obj.NodeDepths(2, cellNo) = double(VSTenum.Snode);
                if isempty(obj.Points(obj.PointNodes==cellNo,:))
                    obj.NodeDepths(2, cellNo) = double(VSTenum.Empty);
                elseif checkNodeValidity(obj, cellNo)
                    obj.NodeDepths(2, cellNo) = double(VSTenum.SnodeL);                
                end
            end
            obj.NodeFrames(newNodeInds,:) = repmat(objNodeFrame, length(newNodeInds), 1);
        end
        
        function [flag, error] = Geometric(obj, nodeNo)
            nodePtMask = obj.PointNodes==nodeNo;
            objNodePoints = obj.Points(nodePtMask,:);            
            objNodeFrame = obj.NodeFrames(nodeNo,:);
            
            a = objNodeFrame(7);
            b = objNodeFrame(8);
            c = objNodeFrame(9);
            d = -objNodeFrame(7:9)*objNodeFrame(10:12)';
            Q = [a^2 a*b a*c a*d;...
                a*b b^2 b*c b*d;...
                a*c b*c c^2 c*d;...
                a*d b*d c*d d^2];
            
            pij = [objNodePoints, ones(size(objNodePoints,1),1)];
            
            error = sum(sum((pij*Q).*pij,2),1);
            flag = error < obj.Properties.errorTh;
        end
        
        function flag = HfieldIndi(obj, cellPoints, cellNormals, cellFrame)
            center = cellFrame(10:12);
            ni = cellFrame(7:9);            
            normaltest = cellNormals*ni';
            
            dist = sqrt(sum(((cellPoints - repmat(center, size(cellPoints,1), 1))*ni').^2,2));            
            maxdist = max(sqrt(sum((cellPoints - repmat(center, size(cellPoints,1), 1)).^2,2)), [], 1);
            
            disttest = dist/(maxdist+eps);
            
            if all((disttest < obj.Properties.delta_d) & (normaltest > obj.Properties.delta_a))
                flag = true;
            else
                flag = false;
            end            
        end
        
        function flag = checkNodeValidity(obj, nodeNo)
            objBounds = obj.NodeBoundaries(nodeNo,:);            
            if obj.NodeDepths(2, nodeNo) == VSTenum.Vnode                
                % In the case of V-node
                minEdgeSize = min(diff(objBounds([1:3;4:6])));
                flag = (obj.NodeDepths(1, nodeNo)+1 >= obj.Properties.maxDepth(1)) || ...
                    minEdgeSize < obj.Properties.minSize || ...
                    nnz(obj.PointNodes==nodeNo) <= obj.Properties.nodeCapacity;
            else
                minEdgeSize = min(diff(objBounds([1:2;4:5])));
                flag = (obj.NodeDepths(1, nodeNo)+1 >= obj.Properties.maxDepth(2)) || ...
                    minEdgeSize < obj.Properties.minSize || ...
                    nnz(obj.PointNodes==nodeNo) <= obj.Properties.nodeCapacity;
            end
        end
        
        function Frame = computeLocalFrame(obj, nodeNo)
            % Gather the new points
            nodePtMask = obj.PointNodes==nodeNo;
            objNodesPoints = obj.Points(nodePtMask,:);
            objNodesNormals = obj.Properties.ptsnormal(nodePtMask,:);
            
            % Do PCA to the normal vector directly
            center = mean(objNodesPoints,1);
            
            [mu, pcs, ~] = obj.nestedPCA(objNodesNormals);
            ni = mu/norm(mu);
            [~, idx] = sort(abs(ni*pcs), 'ascend');
            proj = pcs(:,idx(1:2)) - repmat((ni*pcs(:,idx(1:2)))/norm(ni),3,1).*repmat(ni',1,2);
            locAxis = [proj, ni']; % x,y,z direction
            Frame = [locAxis(:)', center];
        end
        
        function [mu, pcs, lambda] = nestedPCA(obj, X)
            % Centering
            mu = mean(X,1);
            X = X - repmat(mu, size(X,1),1);            
            PHI = X'*X;
            [evec, eval] = eig(PHI);
            [~, sidx] = sort(diag(eval), 'descend');
            pcs = evec(:, sidx);
            lambda = diag(eval(:,sidx)); lambda = lambda(sidx);            
        end
        
        function h = plot(obj, mode, varargin)
            % OcTree.plot plots node bounding boxes of an OcTree object
            %
            % H = OT.plot(mode, 'name',value,...) allows you to specify any
            % properties of the bounding box lines that you would normally
            % supply to a plot(mode, ...,'name',value) command, and returns 
            % plot object handles (one per node) to H.
            % mode (tenary)
            %   1st digit - Visualize volume with none(0) / leaf(1) /
            %   all(2)
            %   2nd digit - Visualize surface with none(0) / leaf(1) /
            %   all(2)
            %   ex) Visualize volume only with leaf node - 01(1)
            %       Visualize surface only with leaf node - 10(3)
            
            VnodeSel = mod(mode,3);
            SnodeSel = ceil(mode/3);
            
            if VnodeSel == 2, selNode = (obj.NodeDepths(2,:) == VSTenum.Vnode) | (obj.NodeDepths(2,:) == VSTenum.Tnode);
            elseif VnodeSel == 1, selNode = obj.NodeDepths(2,:) == VSTenum.Tnode;
            else, selNode = false(1, size(obj.NodeDepths,2));
            end
            
            if SnodeSel == 2, selNode = selNode | (obj.NodeDepths(2,:) == VSTenum.Snode) | (obj.NodeDepths(2,:) == VSTenum.SnodeL);
            elseif SnodeSel == 1, selNode = selNode | (obj.NodeDepths(2,:) == VSTenum.SnodeL);
            else % Do nothing
            end
            
            hold on;
            h = zeros(sum(selNode),1);
            for i=1:obj.NodeCount
                if selNode(i)
                    nodeMinMax = obj.NodeBoundaries(i,:);
                    if (obj.NodeDepths(2,i) == double(VSTenum.Snode)) || ...
                            (obj.NodeDepths(2,i) == double(VSTenum.SnodeL))
                        rotMat = reshape(obj.NodeFrames(i, 1:9), 3, 3);
                        center = obj.NodeFrames(i, 10:12);
                        pts2d = nodeMinMax([1,2,3; 1,5,3; 4,2,6; 4,5,6]);
                        pts3d = rotMat*pts2d' + repmat(center', 1, 4);
                        pts = cat(1, pts3d([...
                            1,2,3; 4,5,6; 10,11,12; 7,8,9; 1,2,3]));
%                         pts = cat(1, nodeMinMax([...
%                             1,2,3; 1,2,6; 4,5,6; 4,5,3; 1,2,3]));
                    else
                        pts = cat(1, nodeMinMax([...
                            1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3;...
                            1 2 6; 4 2 6; 4 5 6; 1 5 6; 1 2 6; 1 2 3]),...
                            nan(1,3), nodeMinMax([4 2 3; 4 2 6]),...
                            nan(1,3), nodeMinMax([4 5 3; 4 5 6]),...
                            nan(1,3), nodeMinMax([1 5 3; 1 5 6]));
                    end
                    h(i) = plot3(pts(:,1),pts(:,2),pts(:,3),varargin{:});
                end
            end
        end
    end
end

