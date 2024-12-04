use raylib::prelude::*;

pub enum PlaneIntersection {
    Line(Line3D),
    VerticalLine(VerticalLine3D),
    Infinite, // In this case the planes contain each other
}

pub enum LineIntersection {
    Point(Vector2),
    VerticalLine(VerticalLine2D),
    Infinite, // In this case the lines are parallel and contain each other
}

pub const fn cross_product(v1: Vector3, v2: Vector3) -> Vector3 {
    let x_component = v1.y * v2.z - v1.z * v2.y;
    let y_component = -1.0 * (v1.x * v2.z - v1.z * v2.x);
    let z_component = v1.x * v2.y - v1.y * v2.x;
    return Vector3::new(x_component, y_component, z_component);
}

pub struct Plane {
    // These are the variables necessary to define a plane in 3 space
    // ax + by + cz = d
    pub a: f32,
    pub b: f32,
    pub c: f32,
    pub d: f32,
    pub x0: f32,
    pub y0: f32,
    pub z0: f32,
}

impl Plane {
    // A plane can be uniquely defined with a normal vector and a point
    pub fn from_point_and_normal(p: Vector3, v: Vector3) -> Self {
        let d = p.x * v.x + p.y * v.y + p.z * v.z;
        Self {
            a: v.x,
            b: v.y,
            c: v.z,
            d,
            x0: p.x,
            y0: p.y,
            z0: p.z,
        }
    }
}

pub struct Line3D {
    // These are the variables necessary to define the parametric equations of a line in 3 space
    // x: x0 + at
    // y: y0 + yt
    // z: z0 + zt
    pub p0: Vector3,
    pub v: Vector3,
}

impl Line3D {
    pub fn from_point_and_parallel_vec(p0: Vector3, v: Vector3) -> Self {
        Self { p0, v }
    }
}

// These are the variables necessary to define vertical lines in the corresponding dimensions
pub struct VerticalLine3D {
    pub x0: f32,
    pub z0: f32,
}

pub struct VerticalLine2D {
    pub x0: f32,
}

impl VerticalLine3D {
    pub fn new(x0: f32, z0: f32) -> Self {
        VerticalLine3D { x0, z0 }
    }
}

impl VerticalLine2D {
    pub fn new(x0: f32) -> Self {
        VerticalLine2D { x0 }
    }
}

pub fn check_planes_intersect(plane1: &Plane, plane2: &Plane) -> bool {
    // Two planes intersect if they are not parallel, or if they are parallel then they exist at the same point
    if plane1.a.abs() != plane2.a.abs()
        || plane1.b.abs() != plane2.b.abs()
        || plane1.c.abs() != plane2.c.abs()
    {
        return true;
    }

    // Check if the planes contain the same point
    // ax0 = d0 = ax1 = d1
    return plane1.d / plane1.x0 == plane2.d / plane2.x0;
}

// This function assumes it has been checked that the planes intersect
pub fn calculate_line_intersection_between_planes(plane1: &Plane, plane2: &Plane) -> PlaneIntersection {
    if !check_planes_intersect(plane1, plane2) {
        panic!("The planes passed tot calculate_line_intersection_between_planes should intersect at some point, and this should be checked");
    }

    let parallel = cross_product(
        Vector3::new(plane1.a, plane1.b, plane1.c),
        Vector3::new(plane2.a, plane2.b, plane2.c),
    );

    // We just need a point, so we need to find a point that both plane1 and plane2 have

    // I'll call this TODO, but for now we can just say that they for sure they share a point a, so
    // let's just say it happens at z = 0
    // Now we have a0x + b0y = d0 = a1x + b1y = d1

    // We now have a system of equations that has a single solution instead of infinitely many
    // solutions (we just picked an arbitrary point)

    let line_intersection =
        find_intersection_of_lines_2d(plane1.a, plane1.b, plane1.d, plane2.a, plane2.b, plane2.d);
    
    match line_intersection {
        LineIntersection::Point(point) => {
            // This is the base case where we have an actual point
            // Remember z0 = 0
            let p0 = Vector3::new(point.x, point.y, 0.0);
            return PlaneIntersection::Line(Line3D::from_point_and_parallel_vec(p0, parallel));
        } 
        LineIntersection::VerticalLine(vertical_line) => {
            return PlaneIntersection::VerticalLine(VerticalLine3D::new(vertical_line.x0, 0.0));
        } 
        LineIntersection::Infinite => {
            return PlaneIntersection::Infinite;
        } 
    }
}

// The first three arguments are the x and y components of the first plane, and d is the d component of
// the plane
// The second three arguments correspond to the second plane
fn find_intersection_of_lines_2d(
    a0: f32,
    b0: f32,
    d0: f32,
    a1: f32,
    b1: f32,
    d1: f32,
) -> LineIntersection {
    // I just did y before x because when I did the math on paper I did it that
    if b0 - ((a0 / a1) / b1) == 0.0 {
        return LineIntersection::Infinite;
    }

    let y0 = (d0 - ((a0 / a1) * d1)) / (b0 - ((a0 / a1) / b1));

    // This would mean that the line points straight up, and because we can't represent this as a Vector2, we have to create an enum variant
    if a0 == 0.0 {
        return LineIntersection::VerticalLine(VerticalLine2D::new(a0));
    }
    let x0 = (d0 - b0 * y0) / a0;

    return LineIntersection::Point(Vector2::new(x0, y0));
}
