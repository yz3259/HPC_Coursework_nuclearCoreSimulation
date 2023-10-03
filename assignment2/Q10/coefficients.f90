!===================================================
! Functions for computing cross-section coefficients
! for 2D diffusion equations
!===================================================

function K1(x,y) result (v)
  real(kind=8), intent(in) :: x,y
  real(kind=8) :: v
  
  v = 2.939d-4 ! region 1
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>5d0/8).and.(y<=7d0/8)) &
     .or.((x>5d0/8).and.(x<=7d0/8).and.(y>2d0/8).and.(y<=5d0/8))  ) v = 4.245d-4 ! region 2
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>7d0/8)) &
     .or.((y>2d0/8).and.(y<=5d0/8).and.(x>7d0/8))  ) v = 4.359d-4 ! region 3
  if (   ((x>5d0/8).and.(x<=7d0/8).and.(y>7d0/8)) &
     .or.((y>5d0/8).and.(y<=7d0/8).and.(x>7d0/8))  ) v = 4.395d-4 ! region 4
  if ((x>5d0/8).and.(x<=7d0/8).and.(y>5d0/8).and.(y<=7d0/8)) v = 4.398d-4 ! region 5
  if ((x>7d0/8).and.(y>7d0/8)) v = 4.415d-4 ! region 6
end function K1
function K2(x,y) result (v)
  real(kind=8), intent(in) :: x,y
  real(kind=8) :: v
  
  v = 1.306d-4 ! region 1
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>5d0/8).and.(y<=7d0/8)) &
     .or.((x>5d0/8).and.(x<=7d0/8).and.(y>2d0/8).and.(y<=5d0/8))  ) v = 1.306d-4 ! region 2
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>7d0/8)) &
     .or.((y>2d0/8).and.(y<=5d0/8).and.(x>7d0/8))  ) v = 1.394d-4 ! region 3
  if (   ((x>5d0/8).and.(x<=7d0/8).and.(y>7d0/8)) &
     .or.((y>5d0/8).and.(y<=7d0/8).and.(x>7d0/8))  ) v = 1.355d-4 ! region 4
  if ((x>5d0/8).and.(x<=7d0/8).and.(y>5d0/8).and.(y<=7d0/8)) v = 1.355d-4 ! region 5
  if ((x>7d0/8).and.(y>7d0/8)) v = 1.345d-4 ! region 6
end function K2

function Sigma_A1(x,y) result (v)
  real(kind=8), intent(in) :: x,y
  real(kind=8) :: v
  
  v = 0.089 ! region 1
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>5d0/8).and.(y<=7d0/8)) &
     .or.((x>5d0/8).and.(x<=7d0/8).and.(y>2d0/8).and.(y<=5d0/8))  ) v = 0.105 ! region 2
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>7d0/8)) &
     .or.((y>2d0/8).and.(y<=5d0/8).and.(x>7d0/8))  ) v = 0.092 ! region 3
  if (   ((x>5d0/8).and.(x<=7d0/8).and.(y>7d0/8)) &
     .or.((y>5d0/8).and.(y<=7d0/8).and.(x>7d0/8))  ) v = 0.091 ! region 4
  if ((x>5d0/8).and.(x<=7d0/8).and.(y>5d0/8).and.(y<=7d0/8)) v = 0.097 ! region 5
  if ((x>7d0/8).and.(y>7d0/8)) v = 0.093 ! region 6
end function SIGMA_A1
function Sigma_A2(x,y) result (v)
  real(kind=8), intent(in) :: x,y
  real(kind=8) :: v
  
  v = 0.109 ! region 1
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>5d0/8).and.(y<=7d0/8)) &
     .or.((x>5d0/8).and.(x<=7d0/8).and.(y>2d0/8).and.(y<=5d0/8))  ) v = 0.025 ! region 2
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>7d0/8)) &
     .or.((y>2d0/8).and.(y<=5d0/8).and.(x>7d0/8))  ) v = 0.093 ! region 3
  if (   ((x>5d0/8).and.(x<=7d0/8).and.(y>7d0/8)) &
     .or.((y>5d0/8).and.(y<=7d0/8).and.(x>7d0/8))  ) v = 0.083 ! region 4
  if ((x>5d0/8).and.(x<=7d0/8).and.(y>5d0/8).and.(y<=7d0/8)) v = 0.098 ! region 5
  if ((x>7d0/8).and.(y>7d0/8)) v = 0.085 ! region 6
end function SIGMA_A2


function Sigma_S(x,y) result (v)
  real(kind=8), intent(in) :: x,y
  real(kind=8) :: v
  
  v = 0d0 ! region 1
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>5d0/8).and.(y<=7d0/8)) &
     .or.((x>5d0/8).and.(x<=7d0/8).and.(y>2d0/8).and.(y<=5d0/8))  ) v = 0d0 ! region 2
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>7d0/8)) &
     .or.((y>2d0/8).and.(y<=5d0/8).and.(x>7d0/8))  ) v = 0.0132 ! region 3
  if (   ((x>5d0/8).and.(x<=7d0/8).and.(y>7d0/8)) &
     .or.((y>5d0/8).and.(y<=7d0/8).and.(x>7d0/8))  ) v = 0.0114 ! region 4
  if ((x>5d0/8).and.(x<=7d0/8).and.(y>5d0/8).and.(y<=7d0/8)) v = 0.0132 ! region 5
  if ((x>7d0/8).and.(y>7d0/8)) v = 0.0114 ! region 6
end function SIGMA_S


function Sigma_F1(x,y) result (v)
  real(kind=8), intent(in) :: x,y
  real(kind=8) :: v
  
  v = 0d0 ! region 1
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>5d0/8).and.(y<=7d0/8)) &
     .or.((x>5d0/8).and.(x<=7d0/8).and.(y>2d0/8).and.(y<=5d0/8))  ) v = 0d0 ! region 2
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>7d0/8)) &
     .or.((y>2d0/8).and.(y<=5d0/8).and.(x>7d0/8))  ) v = 0.140 ! region 3
  if (   ((x>5d0/8).and.(x<=7d0/8).and.(y>7d0/8)) &
     .or.((y>5d0/8).and.(y<=7d0/8).and.(x>7d0/8))  ) v = 0.109 ! region 4
  if ((x>5d0/8).and.(x<=7d0/8).and.(y>5d0/8).and.(y<=7d0/8)) v = 0.124 ! region 5
  if ((x>7d0/8).and.(y>7d0/8)) v = 0.107 ! region 6
end function SIGMA_F1
function Sigma_F2(x,y) result (v)
  real(kind=8), intent(in) :: x,y
  real(kind=8) :: v
  
  v = 0.0079 ! region 1
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>5d0/8).and.(y<=7d0/8)) &
     .or.((x>5d0/8).and.(x<=7d0/8).and.(y>2d0/8).and.(y<=5d0/8))  ) v = 0.0222 ! region 2
  if (   ((x>2d0/8).and.(x<=5d0/8).and.(y>7d0/8)) &
     .or.((y>2d0/8).and.(y<=5d0/8).and.(x>7d0/8))  ) v = 0.0156 ! region 3
  if (   ((x>5d0/8).and.(x<=7d0/8).and.(y>7d0/8)) &
     .or.((y>5d0/8).and.(y<=7d0/8).and.(x>7d0/8))  ) v = 0.0159 ! region 4
  if ((x>5d0/8).and.(x<=7d0/8).and.(y>5d0/8).and.(y<=7d0/8)) v = 0.0151 ! region 5
  if ((x>7d0/8).and.(y>7d0/8)) v = 0.0157 ! region 6
end function SIGMA_F2
