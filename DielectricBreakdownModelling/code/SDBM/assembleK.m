    function K = assembleK(V, E)
        %Assume Rotation matrices given lined up with `vertices
        %assembles K, matrix of undeformed edges that will be used to build
        %covar matrix
        Kijv = [];
        lengths = edge_lengths(V, E);
        eij = V(E(:, 2), :) - V(E(:, 1), :);
        eij = eij./(lengths);
        eji = V(E(:, 1), :) - V(E(:, 2), :);
        eji = eji./(lengths);
        
        eijv = [];
        ejiv = [];
        %K will be |V| x 3 |V| in size
        for i=1:length(E)
            %each edge has two half edges : ij and ji
            
            %half edge 1
            %vert 1 of half edge 1
            eijv = [eijv; E(i, 1), 2*E(i, 1) - 1, -eij(i, 1) ];
            eijv = [eijv; E(i, 1), 2*E(i, 1),-eij(i, 2) ];
            
            %vert 2 of half edge 1
            eijv = [eijv; E(i, 2), 2*E(i, 1) - 1, eij(i, 1)  ];
            eijv = [eijv; E(i, 2), 2*E(i, 1) , eij(i, 2)];
            %
            %half edge 2
            %vert 1 of half edge 2
            ejiv = [ejiv; E(i, 2), 2*E(i, 2) - 1,  -eji(i, 1)  ];
            ejiv = [ejiv; E(i, 2), 2*E(i, 2) ,  -eji(i, 2) ];
            %
            %           %vert 2 of half edge 2
            ejiv = [ejiv; E(i, 1), 2*E(i, 2) - 1,  eji(i, 1)  ];
            ejiv = [ejiv; E(i, 1), 2*E(i, 2),  eji(i, 2) ];
            
        end
        Kijv = [eijv; ejiv];
        
        K = sparse(Kijv(:, 1), Kijv(:, 2), Kijv(:, 3), length(V), 2*length(V));
    end