function make_msh(geo::String; clscale::Real=1, verbose::Bool=false, smooth::Int=5)
    fname = tempname()
    write("$(fname).geo", geo)
    run(`gmsh $(fname).geo -2 -o $(fname).msh -clscale $(clscale) -v $(verbose ? 1 : 0) -smooth $(smooth)`)
    return "$(fname).msh"
end

function make_msh_sqrlat_cylhole(lx::Real, ly::Real, Lx::Real, Ly::Real; clscale::Real=1, verbose::Bool=false, smooth::Int=5)
    @assert lx < Lx
    @assert ly < Ly
    template = """
    lx = $(lx);
    ly = $(ly);
    Lx = $(Lx);
    Ly = $(Ly);
    //+
    Point(1) = {0, 0, 0, $(min(Lx,Ly)/20)};
    Point(2) = {lx/2, 0, 0, $(min(Lx,Ly)/20)};
    Point(3) = {Lx/2, 0, 0, $(min(Lx,Ly)/20)};
    Point(4) = {Lx/2, Ly/2, 0, $(min(Lx,Ly)/20)};
    Point(5) = {0, Ly/2, 0, $(min(Lx,Ly)/20)};
    Point(6) = {0, ly/2, 0, $(min(Lx,Ly)/20)};
    //+
    Line(1) = {1, 2};
    Line(2) = {2, 3};
    Line(3) = {3, 4};
    Line(4) = {4, 5};
    Line(5) = {5, 6};
    Line(6) = {6, 1};
    Ellipse(7) = {2, 1, 1, 6};
    //+
    Line Loop(1) = {2, 3, 4, 5, -7};
    Plane Surface(1) = {1};
    //+
    Field[1] = Distance; // from cylinder surface
    Field[1].NNodesByEdge = 1000;
    Field[1].EdgesList = {7};
    Field[2] = Distance; // from end of simulation domain
    Field[2].NNodesByEdge = 1000;
    Field[2].EdgesList = {3, 4};
    //+
    Field[3] = MathEval;
    Field[3].F = "(F1+F2)/5"; // 5 points per gap
    Background Field = 3;
    //+
    Physical Surface("1") = {1};
    Physical Line("1") = {2};
    """
    return make_msh(template, clscale=clscale, verbose=verbose, smooth=smooth)
end

function make_msh_hexlat_cylhole(lx::Real, ly::Real, Lx::Real, Ly::Real; clscale::Real=1, verbose::Bool=false, smooth::Int=5)
    ly = lx*sqrt(3)#temp
    Ly = Lx*sqrt(3)#temp
    @assert lx < Lx
    @assert ly < Ly
    template = """
    lx = $(lx);
    ly = $(ly);
    Lx = $(Lx);
    Ly = $(Ly);
    //+
    Point(1) = {0, 0, 0, $(min(Lx,Ly)/20)};
    Point(2) = {Lx/2, 0, 0, $(min(Lx,Ly)/20)};
    Point(3) = {Lx/2, Ly/2, 0, $(min(Lx,Ly)/20)};
    Point(4) = {0, Ly/2, 0, $(min(Lx,Ly)/20)};
    Point(5) = {lx/2, 0, 0, $(min(Lx,Ly)/20)};
    Point(6) = {0, lx/2, 0, $(min(Lx,Ly)/20)};
    Point(7) = {Lx/2, Ly/2-lx/2, 0, $(min(Lx,Ly)/20)};
    Point(8) = {Lx/2-lx/2, Ly/2, 0, $(min(Lx,Ly)/20)};
    //+
    Line(1) = {1, 5};
    Line(2) = {5, 2};
    Line(3) = {2, 7};
    Line(4) = {7, 3};
    Line(5) = {3, 8};
    Line(6) = {8, 4};
    Line(7) = {4, 6};
    Line(8) = {6, 1};
    Circle(9) = {5, 1, 6};
    Circle(10) = {7, 3, 8};
    //+
    Line Loop(1) = {2, 3, 10, 6, 7, -9};
    Plane Surface(1) = {1};
    //+
    Field[1] = Distance; // from lower cylinder surface
    Field[1].NNodesByEdge = 100;
    Field[1].EdgesList = {9};
    Field[2] = Distance; // from upper cylinder surface
    Field[2].NNodesByEdge = 100;
    Field[2].EdgesList = {10};
    //+
    Field[3] = MathEval;
    Field[3].F = "(F1+F2)/5"; // 5 points per gap between cyl1 and cyl2
    Field[4] = MathEval;
    Field[4].F = Sprintf("(F1+%g-x)/2.5",Lx/2); // 2.5 points per gap between cyl1 and right end
    Field[5] = MathEval;
    Field[5].F = "(F2+x)/2.5"; // 2.5 points per gap between cyl2 and left end
    Field[6] = Min;
    Field[6].FieldsList = {3, 4, 5};
    Background Field = 6;
    //+
    Physical Surface("1") = {1};
    Physical Line("1") = {2};
    """
    return make_msh(template, clscale=clscale, verbose=verbose, smooth=smooth)
end

function make_msh_sqrlat_cyl(lx::Real, ly::Real, Lx::Real, Ly::Real, d_skin::Real, N_skin::Int; clscale::Real=1, verbose::Bool=false, smooth::Int=5)
    @assert lx < Lx
    @assert ly < Ly
    d_skin = min(d_skin, lx, ly)
    template = """
    lx = $(lx);
    ly = $(ly);
    Lx = $(Lx);
    Ly = $(Ly);
    //+
    Point(1) = {-Lx/2, -Ly/2, 0, $(min(Lx,Ly)/20)};
    Point(2) = {+Lx/2, -Ly/2, 0, $(min(Lx,Ly)/20)};
    Point(3) = {+Lx/2, +Ly/2, 0, $(min(Lx,Ly)/20)};
    Point(4) = {-Lx/2, +Ly/2, 0, $(min(Lx,Ly)/20)};
    Point(5) = {-Lx/2, 0, 0, $(min(Lx,Ly)/20)};
    Point(6) = {+Lx/2, 0, 0, $(min(Lx,Ly)/20)};
    Point(7) = {-lx/2, 0, 0, $(min(Lx,Ly)/20)};
    Point(8) = {0, -ly/2, 0, $(min(Lx,Ly)/20)};
    Point(9) = {+lx/2, 0, 0, $(min(Lx,Ly)/20)};
    Point(10) = {0, +ly/2, 0, $(min(Lx,Ly)/20)};
    Point(11) = {0, 0, 0, $(min(Lx,Ly)/20)};
    //+
    Line(1) = {1, 2};
    Line(2) = {2, 6};
    Line(3) = {6, 3};
    Line(4) = {3, 4};
    Line(5) = {4, 5};
    Line(6) = {5, 1};
    //+
    Ellipse(7) = {7, 11, 11, 8};
    Ellipse(8) = {8, 11, 11, 9};
    Ellipse(9) = {9, 11, 11, 10};
    Ellipse(10) = {10, 11, 11, 7};
    //+
    Line Loop(1) = {1, 2, 3, 4, 5, 6};
    Line Loop(2) = {7, 8, 9, 10};
    Plane Surface(1) = {1, 2};
    Plane Surface(2) = {2};
    //+
    Line(11) = {5, 7};
    Line(12) = {7, 11};
    Line(13) = {11, 9};
    Line(14) = {9, 6};
    Line{11} In Surface{1};
    Line{12} In Surface{2};
    Line{13} In Surface{2};
    Line{14} In Surface{1};
    //+
    Field[1] = Distance; // from cylinder surface
    Field[1].NNodesByEdge = 1000;
    Field[1].EdgesList = {7, 8, 9, 10};
    Field[2] = Distance; // from end of simulation domain
    Field[2].NNodesByEdge = 1000;
    Field[2].EdgesList = {1, 2, 3, 4, 5, 6};
    //+
    Field[3] = MathEval;
    Field[3].F = "(F1+F2)/5"; // 5 points per gap
    //Field[4] = Restrict;
    //Field[4].IField = 3;
    //Field[4].FacesList = {1};
    //Field[4].EdgesList = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14};
    //Field[4].VerticesList = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    //+
    Field[5] = MathEval;
    Field[5].F = "$(d_skin)/$(N_skin)*Exp(F1/$(d_skin))";
    //Field[6] = Restrict;
    //Field[6].IField = 5;
    //Field[6].FacesList = {2};
    //Field[6].EdgesList = {7, 8, 9, 10, 12, 13};
    //Field[6].VerticesList = {7, 8, 9, 10};
    //+

    Field[6] = MathEval;
    Field[6].F = "$(min(Lx,Ly)/20)";

    Field[7] = Min;
    Field[7].FieldsList = {3, 5, 6};
    Background Field = 7;
    //+
    Physical Surface("1") = {1};
    Physical Surface("2") = {2};
    Physical Line("1") = {11, 12, 13, 14};
    //+
    Mesh.CharacteristicLengthExtendFromBoundary = 0;
    """
    return make_msh(template, clscale=clscale, verbose=verbose, smooth=smooth)
end

function make_msh_hexlat_cyl(lx::Real, ly::Real, Lx::Real, Ly::Real, d_skin::Real, N_skin::Int; clscale::Real=1, verbose::Bool=false, smooth::Int=5)
    template = """
    lx = $(lx);
    ly = $(ly);
    Lx = $(Lx);
    Ly = $(Ly);
    //+
    Point(1) = {-Lx/2, -Ly/2, 0, $(min(Lx,Ly)/2)};
    Point(2) = {-Lx/2+lx/2, -Ly/2, 0, $(min(Lx,Ly)/2)};
    Point(3) = {+Lx/2-lx/2, -Ly/2, 0, $(min(Lx,Ly)/2)};
    Point(4) = {+Lx/2, -Ly/2, 0, $(min(Lx,Ly)/2)};
    Point(5) = {+Lx/2, -Ly/2+ly/2, 0, $(min(Lx,Ly)/2)};
    Point(6) = {+Lx/2, 0, 0, $(min(Lx,Ly)/2)};
    Point(7) = {+Lx/2, +Ly/2-ly/2, 0, $(min(Lx,Ly)/2)};
    Point(8) = {+Lx/2, +Ly/2, 0, $(min(Lx,Ly)/2)};
    Point(9) = {+Lx/2-lx/2, +Ly/2, 0, $(min(Lx,Ly)/2)};
    Point(10) = {-Lx/2+lx/2, +Ly/2, 0, $(min(Lx,Ly)/2)};
    Point(11) = {-Lx/2, +Ly/2, 0, $(min(Lx,Ly)/2)};
    Point(12) = {-Lx/2, +Ly/2-ly/2, 0, $(min(Lx,Ly)/2)};
    Point(13) = {-Lx/2, 0, 0, $(min(Lx,Ly)/2)};
    Point(14) = {-Lx/2, -Ly/2+ly/2, 0, $(min(Lx,Ly)/2)};
    Point(15) = {-lx/2, 0, 0, $(min(Lx,Ly)/2)};
    Point(16) = {0, -ly/2, 0, $(min(Lx,Ly)/2)};
    Point(17) = {+lx/2, 0, 0, $(min(Lx,Ly)/2)};
    Point(18) = {0, +ly/2, 0, $(min(Lx,Ly)/2)};
    Point(19) = {0, 0, 0, $(min(Lx,Ly)/2)};
    //+
    Line(1) = {1, 2};
    Line(2) = {2, 3};
    Line(3) = {3, 4};
    Line(4) = {4, 5};
    Line(5) = {5, 6};
    Line(6) = {6, 7};
    Line(7) = {7, 8};
    Line(8) = {8, 9};
    Line(9) = {9, 10};
    Line(10) = {10, 11};
    Line(11) = {11, 12};
    Line(12) = {12, 13};
    Line(13) = {13, 14};
    Line(14) = {14, 1};
    //+
    Ellipse(15) = {2, 1, 1, 14};
    Ellipse(16) = {5, 4, 4, 3};
    Ellipse(17) = {9, 8, 8, 7};
    Ellipse(18) = {12, 11, 11, 10};
    Ellipse(19) = {15, 19, 19, 16};
    Ellipse(20) = {16, 19, 19, 17};
    Ellipse(21) = {17, 19, 19, 18};
    Ellipse(22) = {18, 19, 19, 15};
    //+
    Line Loop(1) = {15, -13, -12, 18, -9, 17, -6, -5, 16, -2};
    Line Loop(2) = {19, 20, 21, 22};
    Line Loop(3) = {14, 1, 15};
    Line Loop(4) = {3, 4, 16};
    Line Loop(5) = {7, 8, 17};
    Line Loop(6) = {10, 11, 18};
    //+
    Plane Surface(1) = {1, 2};
    Plane Surface(2) = {2};
    Plane Surface(3) = {3};
    Plane Surface(4) = {4};
    Plane Surface(5) = {5};
    Plane Surface(6) = {6};
    //+
    Line(23) = {13, 15};
    Line(24) = {15, 19};
    Line(25) = {19, 17};
    Line(26) = {17, 6};
    Line{23} In Surface{1};
    Line{24} In Surface{2};
    Line{25} In Surface{2};
    Line{26} In Surface{1};
    //+

    // Create central cylinder of neighbouring cells
    Point(20) = {0-Lx, -ly/2, 0, 1.0};
    Point(21) = {+lx/2-Lx, 0, 0, 1.0};
    Point(22) = {0-Lx, +ly/2, 0, 1.0};
    Point(23) = {0-Lx, 0, 0, 1.0};
    Ellipse(27) = {21, 23, 23, 22};
    Ellipse(28) = {20, 23, 23, 21};
    Point(24) = {+lx/2, 0-Ly, 0, 1.0};
    Point(25) = {0, +ly/2-Ly, 0, 1.0};
    Point(26) = {-lx/2, 0-Ly, 0, 1.0};
    Point(27) = {0, 0-Ly, 0, 1.0};
    Ellipse(29) = {25, 27, 27, 26};
    Ellipse(30) = {24, 27, 27, 25};
    Point(28) = {0+Lx, +ly/2, 0, 1.0};
    Point(29) = {-lx/2+Lx, 0, 0, 1.0};
    Point(30) = {0+Lx, -ly/2, 0, 1.0};
    Point(31) = {0+Lx, 0, 0, 1.0};
    Ellipse(31) = {29, 31, 31, 30};
    Ellipse(32) = {28, 31, 31, 29};
    Point(32) = {-lx/2, 0+Ly, 0, 1.0};
    Point(33) = {0, -ly/2+Ly, 0, 1.0};
    Point(34) = {+lx/2, 0+Ly, 0, 1.0};
    Point(35) = {0, 0+Ly, 0, 1.0};
    Ellipse(33) = {33, 35, 35, 34};
    Ellipse(34) = {32, 35, 35, 33};

    // Distance to each arc of cylinder surface
    Field[1] = Distance;
    Field[1].NNodesByEdge = 1000;
    Field[1].EdgesList = {15};
    Field[2] = Distance;
    Field[2].NNodesByEdge = 1000;
    Field[2].EdgesList = {16};
    Field[3] = Distance;
    Field[3].NNodesByEdge = 1000;
    Field[3].EdgesList = {17};
    Field[4] = Distance;
    Field[4].NNodesByEdge = 1000;
    Field[4].EdgesList = {18};
    Field[5] = Distance;
    Field[5].NNodesByEdge = 1000;
    Field[5].EdgesList = {19};
    Field[6] = Distance;
    Field[6].NNodesByEdge = 1000;
    Field[6].EdgesList = {20};
    Field[7] = Distance;
    Field[7].NNodesByEdge = 1000;
    Field[7].EdgesList = {21};
    Field[8] = Distance;
    Field[8].NNodesByEdge = 1000;
    Field[8].EdgesList = {22};

    // Distance to *nearest* cylinder surface
    Field[9] = Min;
    Field[9].FieldsList = {1,2,3,4,5,6,7,8};

    // 5 points per gap
    Field[10] = MathEval;
    Field[10].F = "(F1+F2)/5";
    Field[11] = MathEval;
    Field[11].F = "(F2+F3)/5";
    Field[12] = MathEval;
    Field[12].F = "(F3+F4)/5";
    Field[13] = MathEval;
    Field[13].F = "(F4+F1)/5";
    Field[14] = MathEval;
    Field[14].F = "(F1+F5)/5";
    Field[15] = MathEval;
    Field[15].F = "(F2+F6)/5";
    Field[16] = MathEval;
    Field[16].F = "(F3+F7)/5";
    Field[17] = MathEval;
    Field[17].F = "(F4+F8)/5";

    // Distance to each arc of central cylinder in *neighbouring cell*
    Field[18] = Distance;
    Field[18].NNodesByEdge = 1000;
    Field[18].EdgesList = {29};
    Field[19] = Distance;
    Field[19].NNodesByEdge = 1000;
    Field[19].EdgesList = {30};
    Field[20] = Distance;
    Field[20].NNodesByEdge = 1000;
    Field[20].EdgesList = {31};
    Field[21] = Distance;
    Field[21].NNodesByEdge = 1000;
    Field[21].EdgesList = {32};
    Field[22] = Distance;
    Field[22].NNodesByEdge = 1000;
    Field[22].EdgesList = {33};
    Field[23] = Distance;
    Field[23].NNodesByEdge = 1000;
    Field[23].EdgesList = {34};
    Field[24] = Distance;
    Field[24].NNodesByEdge = 1000;
    Field[24].EdgesList = {27};
    Field[25] = Distance;
    Field[25].NNodesByEdge = 1000;
    Field[25].EdgesList = {28};

    // Gaps between central cylinder and central cylinder of next cell
    Field[26] = MathEval;
    Field[26].F = "(F18+F5)/5";
    Field[27] = MathEval;
    Field[27].F = "(F19+F6)/5";
    Field[28] = MathEval;
    Field[28].F = "(F20+F6)/5";
    Field[29] = MathEval;
    Field[29].F = "(F21+F7)/5";
    Field[30] = MathEval;
    Field[30].F = "(F22+F7)/5";
    Field[31] = MathEval;
    Field[31].F = "(F23+F8)/5";
    Field[32] = MathEval;
    Field[32].F = "(F24+F8)/5";
    Field[33] = MathEval;
    Field[33].F = "(F25+F5)/5";

    // Skin depth
    Field[34] = MathEval;
    Field[34].F = "$(d_skin)/$(N_skin)*Exp(F9/$(d_skin))";

    Field[35] = Min;
    Field[35].FieldsList = {10,11,12,13,14,15,16,17,26,27,28,29,30,31,32,33,34};
    Background Field = 35;

    //+
    Physical Surface("1") = {1};
    Physical Surface("2") = {2,3,4,5,6};
    Physical Line("1") = {23, 24, 25, 26};
    //+
    Mesh.CharacteristicLengthExtendFromBoundary = 0;
    """
    return make_msh(template, clscale=clscale, verbose=verbose, smooth=smooth)
end
