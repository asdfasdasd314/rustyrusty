use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};

use crate::heap::custom_heap::*;
use crate::{convert_polygon_to_plane, render::Polygon};

use raylib::prelude::*;

pub const fn cross_product(v1: &Vector3, v2: &Vector3) -> Vector3 {
    let x_component = v1.y * v2.z - v1.z * v2.y;
    let y_component = -1.0 * (v1.x * v2.z - v1.z * v2.x);
    let z_component = v1.x * v2.y - v1.y * v2.x;
    return Vector3::new(x_component, y_component, z_component);
}

pub const fn dot_product_3d(v1: &Vector3, v2: &Vector3) -> f32 {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

pub const fn dot_product_2d(v1: &Vector2, v2: &Vector2) -> f32 {
    return v1.x * v2.x + v1.y * v2.y;
}

pub fn calculate_magnitude_3d(v: &Vector3) -> f32 {
    return (v.x.powi(2) + v.y.powi(2) + v.z.powi(2)).sqrt();
}

pub fn calculate_magnitude_2d(v: &Vector2) -> f32 {
    return (v.x.powi(2) + v.y.powi(2)).sqrt();
}

pub enum PlaneIntersection {
    Line(Line3D),
    Infinite, // In this case the planes contain each other
}

pub enum LineIntersection {
    Point(Vector2),
    Infinite, // In this case the lines are parallel and contain each other
}

pub enum GeometricLine {
    RegularLine(Line3D),
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

#[derive(Debug)]
pub struct Plane {
    // These are the variables necessary to define a plane in 3 space
    // ax + by + cz = d
    pub a: f32,
    pub b: f32,
    pub c: f32,
    pub d: f32,

    // The point (x0, y0, z0) is an absolute point
    pub p0: Vector3,
}

impl Plane {
    // A plane can also be uniquely defined by three points
    pub fn from_three_points(points: (&Vector3, &Vector3, &Vector3)) -> Self {
        let u = *points.2 - *points.0;
        let v = *points.1 - *points.0;

        let normal = cross_product(&u, &v);

        Self {
            a: normal.x,
            b: normal.y,
            c: normal.z,
            d: normal.x * points.0.x + normal.y * points.0.y + normal.z * points.0.z,
            p0: *points.0,
        }
    }

    // A plane can be uniquely defined with a normal vector and a point
    pub fn from_point_and_normal(p0: Vector3, v: Vector3) -> Self {
        let d = p0.x * v.x + p0.y * v.y + p0.z * v.z;
        Self {
            a: v.x,
            b: v.y,
            c: v.z,
            d,
            p0,
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
}

pub fn convert_line_segment_to_line(line_segment: (Vector3, Vector3)) -> Line3D {
    let vec = line_segment.1 - line_segment.0;
    return Line3D::from_point_and_parallel_vec(line_segment.1, vec);
}

pub fn project_point_onto_plane(point: &Vector3, plane: &Plane) -> Vector3 {
    let normal_vec = Vector3::new(plane.a, plane.b, plane.c);
    let d_vec = *point - plane.p0;

    let normal_magnitude = calculate_magnitude_3d(&normal_vec);
    let d_magnitude = calculate_magnitude_3d(&d_vec);

    if d_magnitude == 0.0 {
        return *point;
    }

    // theta is in radians
    let theta = (dot_product_3d(&normal_vec, &d_vec) / (normal_magnitude * d_magnitude)).acos();

    let height = d_magnitude * (90.0 - theta.to_degrees()).to_radians().sin();
    let height_vec = normal_vec * (height / normal_magnitude);

    return *point - height_vec;
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

fn is_counter_clockwise(point1: &Vector2, point2: &Vector2, point3: &Vector2) -> bool {
    return (point2.x - point1.x) * (point3.y - point1.y)
        - (point2.y - point1.y) * (point3.x - point1.x)
        > 0.0;
}

/**
An implementation of Graham's scan that defines the convex polygon from a set of points in O(n) time
The reason the overall algorithm is nlogn is because the sorting takes nlogn
Yes I understand I'm working with a limited number of points, and often the data is nearly sorted so heapsort
and this worrying about time complexity is wrong, but I am doing this project TO LEARN, so it doesn't matter
 */
fn graham_scan(projected_points: &[Vector2]) -> Vec<Vector2> {
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

pub fn calculate_base_plane(plane: &Plane) -> BasePlane {
    // We don't actually need j hat
    let i = Vector3::new(1.0, 0.0, 0.0);
    //let j = Vector3::new(0.0, 1.0, 0.0,);
    let k = Vector3::new(0.0, 0.0, 1.0);

    let plane_normal = Vector3::new(plane.a, plane.b, plane.c);
    // If the dot of the normal and the unit vector that defines the base plane is not the normal, then it's not parallel
    // i -> yz plane, j -> xz, k -> xy

    let magnitude_normal = calculate_magnitude_3d(&plane_normal);
    if dot_product_3d(&plane_normal, &k).abs() == magnitude_normal {
        return BasePlane::XY;
    } else if dot_product_3d(&plane_normal, &i).abs() == magnitude_normal {
        return BasePlane::YZ;
    } else {
        return BasePlane::XZ;
    }
}

pub fn project_into_2d(points: &[Vector3], base_plane: &BasePlane) -> Vec<Vector2> {
    match base_plane {
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
pub fn project_into_3d(plane: &Plane, points: &[Vector2], base_plane: &BasePlane) -> Vec<Vector3> {
    // In two of the cases they are vertical planes, so account for those
    // XZ is the only plane we calculate the points for because it's also the default
    match base_plane {
        BasePlane::XZ => {
            // Find the y coordinate
            points
                .iter()
                .map(|point2d| {
                    let y = (plane.d - point2d.x * plane.a - point2d.y * plane.c) / plane.b;
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
                    let z = plane.d / plane.c;
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
                    let x = plane.d / plane.a;
                    Vector3::new(x, point2d.x, point2d.y)
                })
                .collect()
        }
    }
}

pub fn project_mesh_onto_plane(plane: &Plane, mesh_polygons: &[Polygon]) -> Polygon {
    // Project each individual point onto the plane
    let mut projected_points: Vec<Vector3> = Vec::new();

    for polygon in mesh_polygons {
        for point in &polygon.points {
            projected_points.push(project_point_onto_plane(point, plane));
        }
    }

    // Convert the points to 2D
    let base_plane = calculate_base_plane(plane);
    let mut two_dimensional_points = project_into_2d(&projected_points, &base_plane);

    // Really quickly remove duplicates because there's a very high chance of those existing
    // We have a small number of points, so just O(n^2) search

    //for i in (0..two_dimensional_points.len()).rev() {
    //    let point = two_dimensional_points[i];
    //    for j in 0..i {
    //        if point.x == two_dimensional_points[j].x && point.y == two_dimensional_points[j].y {
    //            two_dimensional_points.remove(i);
    //            break;
    //        }
    //    }
    //}

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

    return Polygon::new(project_into_3d(plane, &bounding_points, &base_plane));
}

pub fn project_polygon_onto_line(line: &Line3D, polygon: &Polygon) -> (Vector3, Vector3) {
    let line_dir = line.v.normalized();
    
    let mut t_min = f32::INFINITY;
    let mut t_max = f32::NEG_INFINITY;
    
    for point in &polygon.points {
        let vec_to_point = *point - line.p0;

        let t = vec_to_point.dot(line_dir);

        if t < t_min {
            t_min = t;
        }
        if t > t_max {
            t_max = t;
        }
    }

    let start_point = line.p0 + line_dir * t_min;
    let end_point = line.p0+  line_dir * t_max;

    (start_point, end_point)
}

/**
Helper function for finding a projection
 */
fn project_point_onto_line(line: &Line3D, point: Vector3) -> Vector3 {
    let v = line.v;
    let u = point - line.p0;
    // u_initial + projv(u) = projected point
    return line.p0 + v * (dot_product_3d(&v, &u) / dot_product_3d(&v, &v));
}

/**
Calculates the orthogonal planes necessary for the separating axis theorem
Some of these planes could be the same as ones already checked, so this isn't quite optimal
 */
pub fn calculate_orthogonal_planes(polygon: &Polygon) -> Vec<Plane> {
    let polygon_as_plane = convert_polygon_to_plane(polygon);
    let normal = Vector3::new(polygon_as_plane.a, polygon_as_plane.b, polygon_as_plane.c);
    let mut planes: Vec<Plane> = Vec::with_capacity(polygon.points.len());
    for edge in polygon.get_edges() {
        let edge_vec = edge.1 - edge.0;
        let new_normal = cross_product(&normal, &edge_vec);
        planes.push(Plane::from_point_and_normal(edge.0, new_normal));
    }
    return planes;
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
    return plane1.d / plane1.p0.x == plane2.d / plane2.p0.x;
}

// This function assumes it has been checked that the planes intersect
pub fn calculate_line_intersection_between_planes(
    plane1: &Plane,
    plane2: &Plane,
) -> PlaneIntersection {
    if !check_planes_intersect(plane1, plane2) {
        panic!("The planes passed tot calculate_line_intersection_between_planes should intersect at some point, and this should be checked");
    }

    let parallel = cross_product(
        &Vector3::new(plane1.a, plane1.b, plane1.c),
        &Vector3::new(plane2.a, plane2.b, plane2.c),
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
            // This is the most basic case where we have an actual point
            // Remember z0 = 0
            let p0 = Vector3::new(point.x, point.y, 0.0);
            return PlaneIntersection::Line(Line3D::from_point_and_parallel_vec(p0, parallel));
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
    let x0 = (d0 - b0 * y0) / a0;

    return LineIntersection::Point(Vector2::new(x0, y0));
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
    let base_vec_magnitude = calculate_magnitude_3d(base_vec);
    let compare_to_vec_magnitude = calculate_magnitude_3d(compare_to_vec);

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

pub fn sort_points_by_angle_from_centroid(points: &mut Vec<Vector3>) {
    let centroid = points
        .iter()
        .fold(Vector3::new(0.0, 0.0, 0.0), |acc, &p| acc + p)
        / points.len() as f32;

    // We need to project everything onto a 2D plane so we can calculate theta and sort by that

    // Calculate normal by crossing any two vectors in the polygon

    let normal: Vector3;

    let p0 = points[0];
    let p1 = points[1];
    let p2 = points[2];

    let v1 = p1 - p0;
    let v2 = p2 - p0;

    // Then we cross the other way
    if dot_product_3d(&v1, &v2) == 0.0 {
        let new_v1 = p2 - p1;
        let new_v2 = p0 - p1;
        normal = cross_product(&new_v1, &new_v2);
    } else {
        normal = cross_product(&v1, &v2);
    }

    let normal_magnitude = calculate_magnitude_3d(&normal);

    let base_vec = p0 - centroid;

    // Use the formula (normal dot v) / (magnitude(normal) * magnitude(v)) = cos(theta)
    points.sort_by(|a, b| {
        let a_vec = *a - centroid;
        let b_vec = *b - centroid;
        let a_dot = dot_product_3d(&a_vec, &base_vec);
        let b_dot = dot_product_3d(&b_vec, &base_vec);
        let a_magnitude = calculate_magnitude_3d(&a_vec);
        let b_magnitude = calculate_magnitude_3d(&b_vec);
        let a_angle = (a_dot / (a_magnitude * normal_magnitude)).acos();
        let b_angle = (b_dot / (b_magnitude * normal_magnitude)).acos();
        a_angle.partial_cmp(&b_angle).expect(&format!(
            "This should never fail? Angle A: {} Angle B: {}",
            a_angle, b_angle
        ))
    });
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
    use proptest::prelude::*;
    use proptest::test_runner::TestRunner;
    use proptest::strategy::{Strategy, ValueTree};
    use super::*;

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
        ((-1000.0..1000.0_f32), (-1000.0..1000.0_f32), (-1000.0..1000.0_f32)).prop_filter("A plane can't have the zero vector as the normal", |(a, b, c)| {
            *a != 0.0 || *b != 0.0 || *c != 0.0
        }).prop_map(move |(a, b, c)| {
            Plane::from_point_and_normal(point, Vector3::new(a, b, c))
        })
    }

    /*
    Proptested functions:
        - project_point_onto_plane
    */
}
