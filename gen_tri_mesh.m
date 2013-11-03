function [p1,p2,p3,p4,p5, pi1,pb1,pi2,pb2,pi3,pb3,pi4,pb4,pi5,pb5, dm1,dm2,dm3,dm4,dm5, np] = gen_tri_mesh(g)

    circ_flag = false;
    if strcmp(g, 'circleg')
        circ_flag = true;
    end
    %pdegplot('squareg'); %,'edgeLabels','on'
    [p1,e1,t1] = initmesh(g, 'Jiggle', 'on', 'Hmax',0.42);
    %pdemesh(p1, e1, t1), axis equal
    [p2,e2,t2] = refinemesh(g,p1,e1,t1);
    [p3,e3,t3] = refinemesh(g,p2,e2,t2);
    [p4,e4,t4] = refinemesh(g,p3,e3,t3);
    [p5,e5,t5] = refinemesh(g,p4,e4,t4);
    
    dm1 = DistanceMatrix(p1', p1');
    dm2 = DistanceMatrix(p2', p2');
    dm3 = DistanceMatrix(p3', p3');
    dm4 = DistanceMatrix(p4', p4');
    dm5 = zeros(length(p5), length(p5)); %
    %dm5 = get_dm(p5); %
    %dm5 = DistanceMatrix(p5', p5');
    %dist2C = bsxfun(@minus, point_set, point_set(:,center));
    %dmc = dot(dist2C,dist2C);
    
    [pi1, pb1, npi1, npb1] = get_i_e_pts(p1, circ_flag);
    [pi2, pb2, npi2, npb2] = get_i_e_pts(p2, circ_flag);
    [pi3, pb3, npi3, npb3] = get_i_e_pts(p3, circ_flag);
    [pi4, pb4, npi4, npb4] = get_i_e_pts(p4, circ_flag);
    [pi5, pb5, npi5, npb5] = get_i_e_pts(p5, circ_flag);
    
    np = [[npi1, npb1]; [npi2, npb2]; [npi3, npb3]; [npi4, npb4]; [npi5, npb5]];

end

function dm = get_dm(point_set)
    dm = zeros(length(point_set), length(point_set)); %DistanceMatrix(p5', p5');
    for center=1:length(point_set)
        dmc = bsxfun(@minus, point_set, point_set(:,center));
        dm(center,:) = dot(dmc,dmc);
    end
    %dmc = dot(dist2C,dist2C);
end

function [int_pts, bound_pts, npi, npb] = get_i_e_pts(points, circ_flag)
    if circ_flag
        int_pts = find(abs(dot(points,points)-1)>1e-14);
        bound_pts = setdiff([1:length(points')], int_pts);
    else
        bound_pts = union(find(points(1,:)==1), find(points(1,:)==-1));
        bound_pts = union(bound_pts, find(points(2,:)==1));
        bound_pts = union(bound_pts, find(points(2,:)==-1));
        int_pts = setdiff([1:length(points')], bound_pts);
    end
    npi = length(int_pts);
    npb = length(bound_pts);
end

function d1 = get_disk_with_hole()

gd = [4.000000000000000   3.000000000000000; ...
                   0   4.000000000000000; ...
                   0  -0.400000000000000; ...
   1.000000000000000   0.400000000000000; ...
   1.000000000000000   0.400000000000000; ...
                   0  -0.400000000000000; ...
                   0   0.400000000000000; ...
                   0   0.400000000000000; ...
                   0  -0.400000000000000; ...
                   0  -0.400000000000000];

sf = [E1-R1];
ns = [69    82; ...
    49    49];

d1=decsg(gd,sf,ns)

end