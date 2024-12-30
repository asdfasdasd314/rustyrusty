use std::cmp;
use crate::heap::custom_heap::*;
use crate::render::{MeshShape, Polygon};

use raylib::prelude::*;

pub const fn f32_min(a: f32, b: f32) -> f32 {
    if a < b {
        return a;
    }
    return b;
}

pub enum PlaneIntersection {
    Line(Line3D),
    Infinite, // In this case the planes contain each other
    None,
}

pub enum LineIntersection {
    Point(Vector2),
    Infinite, // In this case the lines are parallel and contain each other
    None,
}

// Theta sweeps across the x-z plane starting at the positive x-axis and goes up to the positive y-axis
// Phi comes down from the positive y-axis
pub struct SphericalAngle {
    pub theta: f32,
    pub phi: f32,
}

impl SphericalAngle {
    pub const fn new(theta: f32, phi: f32) -> SphericalAngle {
        SphericalAngle { theta, phi }
    }
}

#[derive(Debug, Clone)]
pub struct Plane {
    // These are the variables necessary to define a plane in 3 space
    // ax + by + cz = d
    // n = <a, b, c>
    pub n: Vector3,
    pub d: f32,

    // The point (x0, y0, z0) is an absolute point
    pub p0: Vector3,
}

impl Plane {
    // A plane can also be uniquely defined by three points
    pub fn from_three_points(points: (&Vector3, &Vector3, &Vector3)) -> Self {
        let u = *points.2 - *points.0;
        let v = *points.1 - *points.0;

        let normal = u.cross(v);

        Self {
            n: normal,
            d: normal.x * points.0.x + normal.y * points.0.y + normal.z * points.0.z,
            p0: *points.0,
        }
    }

    // A plane can be uniquely defined with a normal vector and a point
    pub fn from_point_and_normal(p0: Vector3, v: Vector3) -> Self {
        let d = p0.x * v.x + p0.y * v.y + p0.z * v.z;
        Self {
            n: v,
            d,
            p0,
        }
    }

    /**
    Projects a point in 3-space onto the plane
     */
    pub fn project_point(&self, point: &Vector3) -> Vector3 {
        let d_vec = *point - self.p0;

        let normal_magnitude = self.n.length();
        let d_magnitude = d_vec.length();

        if d_magnitude == 0.0 {
            return *point;
        }

        // theta is in radians
        let theta = (self.n.dot(d_vec) / (normal_magnitude * d_magnitude)).acos();

        let height = d_magnitude * (90.0 - theta.to_degrees()).to_radians().sin();
        let height_vec = self.n * (height / normal_magnitude);

        return *point - height_vec;
    }

    pub fn project_mesh(&self, mesh: &Box<dyn MeshShape>) -> Polygon {
        // Project each individual point onto the plane
        let mut projected_points: Vec<Vector3> = Vec::new();

        for polygon in mesh.get_polygons() {
            for point in &polygon.points {
                projected_points.push(self.project_point(point));
            }
        }

        // Convert the points to 2D
        let point_projector = TwoDimensionalPointProjector::new(self.clone());
        let mut two_dimensional_points = point_projector.project_into_2d(&projected_points);

        // Get the axis for which all points will be compared to
        let comparison_axis = ComparisonAxis::new(&two_dimensional_points);

        // Perform Graham's check to get the points of the polygon

        // Sort based on the root pointt
        heapsort(&comparison_axis, &mut two_dimensional_points);

        // I implemented graham's scan in such a way that sometimes the root isn't included
        let mut bounding_points = graham_scan(&two_dimensional_points);

        // It doesn't matter what order they are in the end at this point
        if bounding_points[0] != comparison_axis.p1 {
            bounding_points.push(comparison_axis.p1);
        }

        return Polygon::new(point_projector.project_into_3d(&bounding_points));
    }

    pub fn intersects(&self, other: &Plane) -> bool {
        // Two planes intersect if they are not parallel, or if they are parallel then they exist at the same point
        if self.n.x.abs() != other.n.x.abs()
            || self.n.y.abs() != other.n.y.abs()
            || self.n.z.abs() != other.n.z.abs()
        {
            return true;
        }

        // Check if the planes contain the same point
        // ax0 = d0 = ax1 = d1
        return self.d / self.p0.x == other.d / other.p0.x;
    }

    /**
    This function is basically useless, because it was part of the old collision detection but I don't really want to remove it
    This function assumes it has been checked that the planes intersect
     */
    pub fn calculate_line_intersection(&self, other: &Plane) -> PlaneIntersection {
        if !self.intersects(other) {
            panic!("The planes passed tot calculate_line_intersection_between_planes should intersect at some point, and this should be checked");
        }

        let parallel = self.n.cross(other.n);

        // We just need a point, so we need to find a point that both plane1 and plane2 have

        // I'll call this TODO, but for now we can just say that they for sure they share a point a, so
        // let's just say it happens at z = 0
        // Now we have a0x + b0y = d0 = a1x + b1y = d1

        // We now have a system of equations that has a single solution instead of infinitely many
        // solutions (we just picked an arbitrary point)

        let line1: Line2D;
        let line2: Line2D;

        if self.n.x == 0.0 {
            line1 = Line2D::from_two_points(Vector2::new(0.0, self.n.y), Vector2::new(1.0, self.n.y));
        } else {
            line1 = Line2D::from_two_points(Vector2::new(0.0, self.n.y), Vector2::new(self.n.x, 0.0));
        }
        if other.n.x == 0.0 {
            line2 = Line2D::from_two_points(Vector2::new(0.0, other.n.y), Vector2::new(1.0, other.n.y));
        } else {
            line2 = Line2D::from_two_points(Vector2::new(0.0, other.n.y), Vector2::new(other.n.x, 0.0));
        }

        let line_intersection = line1.find_intersection(&line2);

        match line_intersection {
            LineIntersection::Point(point) => {
                // This is the most basic case where we have an actual point
                // Remember z0 = 0
                let p0 = Vector3::new(point.x, point.y, 0.0);
                return PlaneIntersection::Line(Line3D::from_point_and_parallel_vec(p0, parallel));
            }
            LineIntersection::Infinite => {
                return PlaneIntersection::Infinite;
            }
            LineIntersection::None => {
                // We should have already panicked, there is no conceivable reason to end up here
                panic!("There is no possible way the logic should end up here because this test case was tested for twice now");
            }
        }
    }
}

#[derive(Debug)]
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

    pub fn from_line_segment(line_segment: (Vector3, Vector3)) -> Self {
        let vec = line_segment.1 - line_segment.0;
        return Line3D::from_point_and_parallel_vec(line_segment.1, vec);
    }

    pub fn project_polygon_onto_line(&self, polygon: &Polygon) -> (Vector3, Vector3) {
        let line_dir = self.v.normalized();

        let mut t_min = f32::INFINITY;
        let mut t_max = f32::NEG_INFINITY;

        for point in &polygon.points {
            let vec_to_point = *point - self.p0;

            let t = vec_to_point.dot(line_dir);

            if t < t_min {
                t_min = t;
            }
            if t > t_max {
                t_max = t;
            }
        }

        let start_point = self.p0 + line_dir * t_min;
        let end_point = self.p0 + line_dir * t_max;

        (start_point, end_point)
    }

    /**
    Helper function for finding a projection
     */
    fn project_point_onto_line(&self, point: Vector3) -> Vector3 {
        let u = point - self.p0;
        // u_initial + projv(u) = projected point
        return self.p0 + self.v * self.v.dot(u) / self.v.dot(self.v);
    }
}

// To avoid issues with divide by zero and whatnot we can just define every single like this where a line can uniquely be defined by any point on the line and a vector in the direction of the line
// This also works 100% perfect for vertical lines so there are no special cases, it is a bit more memory though (only maybe like 8 bytes though I think, 3 f32s vs 4 f32s + a little bit for pointers and whatnot)
#[derive(Debug)]
pub struct Line2D {
    p0: Vector2,
    v: Vector2,
}

impl Line2D {
    pub fn from_two_points(point1: Vector2, point2: Vector2) -> Self {
        Self {
            p0: point1,
            v: point2 - point1,
        }
    }

    pub fn find_intersection(&self, other: &Line2D) -> LineIntersection {
        // Check if the vector is parallel
        if self.v.dot(other.v.normalized()) == self.v.length() {
            // If the dot product of one v and the other one normalized is the magnitude of the first one then they are parallel because math

            let projected_tail = other.project_point(self.p0 + self.v);

            // If this is true then the projection did nothing, meaning that the lines are overlapping
            if projected_tail == self.p0 + self.v {
                return LineIntersection::Infinite; 
            }
            else {
                return LineIntersection::None;
            }
        }
        let displacement_line = Line2D::from_two_points(self.p0, other.p0);
        let t: f32;
        // First check which projection to use
        if displacement_line.v.dot(self.v) == 0.0 {
            // If the lines are parallel project onto self.v
            let proj_other_to_self = self.project_point(other.v);
            t = proj_other_to_self.length() / self.v.length();
        }
        else {
            // It's theoretically possible we go past the point, so pick two points equidistant away and choose the direction that gets closer
            let point1 = self.p0 + self.v;
            let point2 = self.p0 - self.v;

            let proj_p1 = other.project_point(point1);
            let proj_p2 = other.project_point(point2);

            let p1_magnitude = (proj_p1 - point1).length();
            let p2_magnitude = (proj_p2 - point2).length();
            t = (p2_magnitude - p1_magnitude) / f32_min(p1_magnitude, p2_magnitude);
        }

        return LineIntersection::Point(self.p0 + self.v * t);
    }

    fn project_point(&self, point: Vector2) -> Vector2 {
        let u = point - self.p0;
        // u_initial + projv(u) = projected point
        return self.p0 + self.v * self.v.dot(u) / self.v.dot(self.v);
    }
}

// p0 is the bottom-leftmost point that can contain the polygon, this is usually a combination of the leftmost and bottommost points, but in the case that they are the same then it's that point 1 unit left in the x direction
// p1 is the root point of the polygon
#[derive(Debug)]
pub struct ComparisonAxis {
    pub p0: Vector2,
    pub p1: Vector2,
}

impl ComparisonAxis {
    pub fn new(points: &[Vector2]) -> ComparisonAxis {
        let mut min_x_index = 0;
        let mut min_y_index = 0;
        for i in 0..points.len() {
            if points[i].y < points[min_y_index].y {
                min_y_index = i;
            }
            // Otherwise compare the x
            else if points[i].y == points[min_y_index].y && points[i].x < points[min_y_index].x {
                min_y_index = i;
            }

            // We don't any more complicated logic for this because we just need the x-coordinate
            if points[i].x < points[min_x_index].x {
                min_x_index = i;
            }
        }

        let p1;
        let p0;
        if min_x_index == min_y_index {
            p1 = points[min_y_index];
            p0 = Vector2::new(p1.x - 1.0, p1.y);
        } else {
            p1 = points[min_y_index];
            p0 = Vector2::new(points[min_x_index].x, p1.y);
        }

        ComparisonAxis { p0, p1 }
    }
}

/**
An implementation of Graham's scan that defines the convex polygon from a set of points in O(n) time
The reason the overall algorithm is nlogn is because the sorting takes nlogn
Yes I understand I'm working with a limited number of points, and often the data is nearly sorted so heapsort
and this worrying about time complexity is wrong, but I am doing this project TO LEARN, so it doesn't matter
 */
fn graham_scan(projected_points: &[Vector2]) -> Vec<Vector2> {
    fn is_counter_clockwise(point1: &Vector2, point2: &Vector2, point3: &Vector2) -> bool {
        return (point2.x - point1.x) * (point3.y - point1.y)
            - (point2.y - point1.y) * (point3.x - point1.x)
            > 0.0;
    }

    let mut stack: Vec<Vector2> = Vec::new();

    for point in projected_points {
        // The check for the length greater than 1 is done so it's possible to verify the angle
        while stack.len() > 1
            && !is_counter_clockwise(&stack[stack.len() - 2], &stack[stack.len() - 1], &point)
        {
            stack.pop();
        }

        stack.push(*point);
    }

    return stack;
}

// Planes that can be described by a single basis vector
#[derive(Debug)]
pub enum BasePlane {
    XZ,
    YZ,
    XY,
}

pub struct TwoDimensionalPointProjector {
    original_plane: Plane,
    base_plane: BasePlane,
}

impl TwoDimensionalPointProjector {
    pub fn new(plane: Plane) -> Self {
        // We don't actually need j hat
        let i = Vector3::new(1.0, 0.0, 0.0);
        //let j = Vector3::new(0.0, 1.0, 0.0,);
        let k = Vector3::new(0.0, 0.0, 1.0);

        // If the dot of the normal and the unit vector that defines the base plane is not the normal, then it's not parallel
        // i -> yz plane, j -> xz, k -> xy

        let magnitude_normal = plane.n.length();
        let base_plane: BasePlane;
        if plane.n.dot(k).abs() == magnitude_normal {
            base_plane = BasePlane::XY;
        } else if plane.n.dot(i).abs() == magnitude_normal {
            base_plane = BasePlane::YZ;
        } else {
            base_plane = BasePlane::XZ;
        }

        Self {
            original_plane: plane,
            base_plane,
        }
    }

    pub fn project_into_2d(&self, points: &[Vector3]) -> Vec<Vector2> {
        match self.base_plane {
            BasePlane::XY => {
                // We map onto the x y plane, so remove the z coordinate
                points
                    .iter()
                    .map(|point| Vector2::new(point.x, point.y))
                    .collect()
            }
            BasePlane::XZ => {
                // We map onto the x z plane, so remove the y coordinate
                points
                    .iter()
                    .map(|point| Vector2::new(point.x, point.z))
                    .collect()
            }
            BasePlane::YZ => {
                // We map onto the y z plane, so remove the x coordinate
                points
                    .iter()
                    .map(|point| Vector2::new(point.y, point.z))
                    .collect()
            }
        }
    }

    /**
    Takes the three initial points that can be used to define the points, the actual points that were mapped to 2D, and the plane that all the points were mapped onto to calculate the positions of the points in 3-space
     */
    pub fn project_into_3d(&self, points: &[Vector2]) -> Vec<Vector3> {
        // In two of the cases they are vertical planes, so account for those
        // XZ is the only plane we calculate the points for because it's also the default
        match self.base_plane {
            BasePlane::XZ => {
                // Find the y coordinate
                points
                    .iter()
                    .map(|point2d| {
                        let y = (self.original_plane.d
                            - point2d.x * self.original_plane.n.x
                            - point2d.y * self.original_plane.n.z)
                            / self.original_plane.n.y;
                        Vector3::new(point2d.x, y, point2d.y)
                    })
                    .collect()
            }

            // In this case calculate the z and put it at the end using the fact that only the z coordinate matters in the equation of the plane
            BasePlane::XY => {
                // Find the z coordinate
                points
                    .iter()
                    .map(|point2d| {
                        let z = self.original_plane.d / self.original_plane.n.z;
                        Vector3::new(point2d.x, point2d.y, z)
                    })
                    .collect()
            }

            // In this case calculate the x and put it at the end using the fact that only the x coordinate matters in the equation of the plane
            BasePlane::YZ => {
                // Find the x coordinate
                points
                    .iter()
                    .map(|point2d| {
                        let x = self.original_plane.d / self.original_plane.n.x;
                        Vector3::new(x, point2d.x, point2d.y)
                    })
                    .collect()
            }
        }
    }
}

/**
The two line segments passed in should be colinear
 */
pub fn line_segments_overlap(
    line_segment1: (Vector3, Vector3),
    line_segment2: (Vector3, Vector3),
) -> bool {
    let a1 = line_segment1.0;
    let b1 = line_segment1.1;
    let a2 = line_segment2.0;
    let b2 = line_segment2.1;

    // Direction vector of segment 1
    let v1 = b1 - a1;

    // Check if Segment 2 is a degenerate point
    if a2 == b2 {
        // Segment 2 is a point: Check if this point lies on Segment 1
        let v_point = a2 - a1;
        let is_collinear = v1.cross(v_point).length() == 0.0;

        if is_collinear {
            // Check if the point lies within Segment 1's bounds
            let dot_product = v_point.dot(v1);
            let is_within_bounds = dot_product >= 0.0 && dot_product <= v1.dot(v1);
            return is_within_bounds;
        }
        return false;
    }

    // Check if Segment 1 is a degenerate point
    if a1 == b1 {
        // Segment 1 is a point: Check if this point lies on Segment 2
        let v_point = a1 - a2;
        let v2 = b2 - a2;
        let is_collinear = v2.cross(v_point).length() == 0.0;

        if is_collinear {
            // Check if the point lies within Segment 2's bounds
            let dot_product = v_point.dot(v2);
            let is_within_bounds = dot_product >= 0.0 && dot_product <= v2.dot(v2);
            return is_within_bounds;
        }
        return false;
    }

    // General case: Both are actual segments
    let v2 = a2 - a1;
    let v3 = b2 - a1;

    // Check collinearity
    if v1.cross(v2).length() != 0.0 || v1.cross(v3).length() != 0.0 {
        return false; // Not collinear
    }

    // Parametrize A2 and B2 on Segment 1's line
    let t1 = (a2 - a1).dot(v1) / v1.dot(v1);
    let t2 = (b2 - a1).dot(v1) / v1.dot(v1);

    // Check if intervals overlap
    let t_min = t1.min(t2);
    let t_max = t1.max(t2);

    t_max >= 0.0 && t_min <= 1.0
}

// Calculates the difference in between the angles of the vectors using the base_vec as the base
// **in degrees**
pub fn calculate_difference_in_angle(
    base_vec: &Vector3,
    compare_to_vec: &Vector3,
) -> SphericalAngle {
    // To make things easier, we'll transform these vectors so that base_vec is on the z-axis
    // We can do this by calculating the angle the base_vec makes with the origin
    // We know x, y, and z, so we can use conversions between ro, theta, and phi to find these changes in angle
    let base_vec_magnitude = base_vec.length();
    let compare_to_vec_magnitude = compare_to_vec.length();

    let base_vec_theta = (base_vec.y / (base_vec.x.powi(2) + base_vec.y.powi(2)).sqrt()).asin();
    let base_vec_phi = base_vec.z / base_vec_magnitude;

    let compare_to_vec_theta =
        (compare_to_vec.y / (compare_to_vec.x.powi(2) + compare_to_vec.y.powi(2)).sqrt()).asin();
    let compare_to_vec_phi = compare_to_vec.z / compare_to_vec_magnitude;

    return SphericalAngle {
        theta: compare_to_vec_theta - base_vec_theta,
        phi: compare_to_vec_phi - base_vec_phi,
    };
}

// This function checks that given some angles, they surround the relative origin they were calculated from
pub fn check_angles_surround_relative_origin(angles: Vec<SphericalAngle>) -> bool {
    // TODO For now, brute force
    for i in 0..(angles.len() - 1) {
        for j in (i + 1)..angles.len() {
            // We want to check if the sign flips across the pair in both the theta and phi
            let angle1 = &angles[i];
            let angle2 = &angles[j];

            // Sign should be flipped, so check for that (the product should be less than 0)
            if angle1.theta * angle2.theta <= 0.0 && angle1.phi * angle2.phi <= 0.0 {
                return true;
            }
        }
    }

    return false;
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::strategy::{Strategy, ValueTree};
    use proptest::test_runner::TestRunner;

    fn point_strategy() -> impl Strategy<Value = Vector3> {
        let x = (-1000.0..1000.0_f32).prop_map(|x| x);
        let y = (-1000.0..1000.0_f32).prop_map(|y| y);
        let z = (-1000.0..1000.0_f32).prop_map(|z| z);
        (x, y, z).prop_map(|(x, y, z)| Vector3::new(x, y, z))
    }

    fn plane_strategy() -> impl Strategy<Value = Plane> {
        // Generate a random value for the point, yeah I could just generate a random value manually, but I'm using the features of proptest!
        let mut runner = TestRunner::default();
        let point_strategy = point_strategy();
        let point: Vector3 = point_strategy.new_tree(&mut runner).unwrap().current();
        (
            (-1000.0..1000.0_f32),
            (-1000.0..1000.0_f32),
            (-1000.0..1000.0_f32),
        )
            .prop_filter(
                "A plane can't have the zero vector as the normal",
                |(a, b, c)| *a != 0.0 || *b != 0.0 || *c != 0.0,
            )
            .prop_map(move |(a, b, c)| Plane::from_point_and_normal(point, Vector3::new(a, b, c)))
    }

    /*
    Proptested functions:
        - project_point_onto_plane
    */
}
