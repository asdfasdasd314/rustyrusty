use std::convert;

use crate::{convert_polygon_to_plane, render::Polygon};
use raylib::prelude::*;

pub const fn cross_product(v1: Vector3, v2: Vector3) -> Vector3 {
    let x_component = v1.y * v2.z - v1.z * v2.y;
    let y_component = -1.0 * (v1.x * v2.z - v1.z * v2.x);
    let z_component = v1.x * v2.y - v1.y * v2.x;
    return Vector3::new(x_component, y_component, z_component);
}

pub const fn dot_product(v1: &Vector3, v2: &Vector3) -> f32 {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

pub fn calculate_magnitude_3d(v: Vector3) -> f32 {
    return (v.x.powi(2) + v.y.powi(2) + v.z.powi(2)).sqrt();
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
    pub fn new(theta: f32, phi: f32) -> SphericalAngle {
        SphericalAngle { theta, phi }
    }
}

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

    let normal_magnitude = calculate_magnitude_3d(normal_vec);
    let d_magnitude = calculate_magnitude_3d(d_vec);

    // theta is in radians
    let theta = (dot_product(&normal_vec, &d_vec) / (normal_magnitude * d_magnitude)).acos();

    let height = d_magnitude * (90.0 - theta.to_degrees());
    let height_vec = normal_vec * (height / normal_magnitude);

    return *point - height_vec;
}

struct RootPoint {
    pub true_root: Vector3,
    // This "helper root" is necessary when calculating the angle from the true root to something else
    pub helper_root: Vector3,
}
fn calculate_root_point(projected_points: &[Vector3]) -> RootPoint {
    let mut min_point_index = 0;
    for i in 0..projected_points.len() {
        if projected_points[i].y <= projected_points[min_point_index].y {
            if projected_points[i].y < projected_points[min_point_index].y {
                min_point_index = i;
            }
            // Otherwise compare the x and z values
            else {
                if projected_points[i].z <= projected_points[min_point_index].z {
                    if projected_points[i].z < projected_points[min_point_index].z {
                        min_point_index = i;
                    } else {
                        if projected_points[i].x < projected_points[min_point_index].x {
                            min_point_index = i;
                        }
                    }
                }
            }
        }
    }

    return RootPoint {
        true_root: projected_points[min_point_index],
        helper_root: Vector3::new(
            projected_points[min_point_index].x - 1.0,
            projected_points[min_point_index].y,
            projected_points[min_point_index].z,
        ),
    };
}

fn graham_scan(projected_points: &[Vector3]) -> Polygon {
    todo!()
}

pub fn project_mesh_onto_plane(plane: &Plane, mesh_polygons: &[Polygon]) -> Polygon {
    // Project each individual point onto the plane
    let mut projected_points: Vec<Vector3> = Vec::new();

    for polygon in mesh_polygons {
        for point in &polygon.points {
            projected_points.push(project_point_onto_plane(point, plane));
        }
    }

    // Find the point with the lowest y-value when projected onto the x-y plane
    // It's not clear what the "bottom left point" is if the plane is three dimensional, and sometimes the minimum x-value is also necessary, so we project the points onto an actual 2D plane
    // Where there is no z coordinate
    let root_point = calculate_root_point(&projected_points);

    // Perform Graham's check to get the points of the polygon

    // Sort based on the root pointt
    sort_points_by_angle_from_root(&mut projected_points, &root_point);

    return graham_scan(&projected_points);
}

pub fn project_polygon_onto_line(line: &Line3D, polygon: &Polygon) -> (Vector3, Vector3) {
    let points = polygon.points.clone();
    // Two pointers  search for the furthest two points on the line when projected
    let mut left = 0;
    let mut right = points.len() - 1;

    let mut max_distance: f32 = calculate_magnitude_3d(
        project_point_onto_line(line, points[right]) - project_point_onto_line(line, points[left]),
    );

    // Move the left pointer to the right and go until the distance is less
    while left < right {
        let curr_distance = calculate_magnitude_3d(
            project_point_onto_line(line, points[right])
                - project_point_onto_line(line, points[left]),
        );
        if curr_distance > max_distance {
            max_distance = curr_distance;
            left += 1;
        } else {
            // This is a necessary step, because it's not getting any better on the left side
            right -= 1;
            break;
        }
    }

    // Move the right pointer to the left and go until the distance is less
    while left < right {
        let curr_distance = calculate_magnitude_3d(
            project_point_onto_line(line, points[right])
                - project_point_onto_line(line, points[left]),
        );
        if curr_distance > max_distance {
            max_distance = curr_distance;
            right -= 1;
        } else {
            break;
        }
    }

    return (points[left], points[right]);
}

/**
Helper function for finding a projection
 */
fn project_point_onto_line(line: &Line3D, point: Vector3) -> Vector3 {
    let v = line.v;
    let u = point - line.p0;
    // u_initial + projv(u) = projected point
    return point + v * (dot_product(&v, &u) / dot_product(&v, &v));
}

/**
Calculates the line orthogonal to the line segment relative to a coplanar point
I'll just note that technically this is completely unnecessary any of the infinite orthogonal lines will do, but I think geometrically this is just more accurate
 */
pub fn calculate_orthogonal_line(
    line_segment: &(Vector3, Vector3),
    third_point: &Vector3,
) -> Line3D {
    let parallel_line =
        Line3D::from_point_and_parallel_vec(*third_point, line_segment.1 - line_segment.0);
    let q = project_point_onto_line(&parallel_line, line_segment.1);
    return Line3D::from_point_and_parallel_vec(line_segment.1, q - line_segment.1);
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
        let new_normal = cross_product(normal, edge_vec);
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
    let line_segment1_magnitude = calculate_magnitude_3d(line_segment1.1 - line_segment1.0);
    let distance1 = calculate_magnitude_3d(line_segment1.1 - line_segment2.0);
    let distance2 = calculate_magnitude_3d(line_segment1.1 - line_segment2.1);
    return distance1 < line_segment1_magnitude || distance2 < line_segment1_magnitude;
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
    let base_vec_magnitude = calculate_magnitude_3d(*base_vec);
    let compare_to_vec_magnitude = calculate_magnitude_3d(*compare_to_vec);

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
    if dot_product(&v1, &v2) == 0.0 {
        let new_v1 = p2 - p1;
        let new_v2 = p0 - p1;
        normal = cross_product(new_v1, new_v2);
    } else {
        normal = cross_product(v1, v2);
    }

    let normal_magnitude = calculate_magnitude_3d(normal);

    let base_vec = p0 - centroid;

    // Use the formula (normal dot v) / (magnitude(normal) * magnitude(v)) = cos(theta)
    points.sort_by(|a, b| {
        let a_vec = *a - centroid;
        let b_vec = *b - centroid;
        let a_dot = dot_product(&a_vec, &base_vec);
        let b_dot = dot_product(&b_vec, &base_vec);
        let a_magnitude = calculate_magnitude_3d(a_vec);
        let b_magnitude = calculate_magnitude_3d(b_vec);
        let a_angle = (a_dot / (a_magnitude * normal_magnitude)).acos();
        let b_angle = (b_dot / (b_magnitude * normal_magnitude)).acos();
        a_angle.partial_cmp(&b_angle).expect(&format!(
            "This should never fail? Angle A: {} Angle B: {}",
            a_angle, b_angle
        ))
    });
}

fn heap_sort_helper(points: &mut Vec<Vector3>) {}

pub fn sort_points_by_angle_from_root(points: &mut Vec<Vector3>, root_point: &Vector3) {}

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
