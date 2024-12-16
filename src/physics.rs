use crate::math_util::*;
use crate::render::*;
use raylib::prelude::*;

pub trait Dynamic {
    // Moves the object and its components by change
    fn move_by(&mut self, change: Vector3);
}

// This object is a basic object that interacts with physics
#[derive(Debug)]
pub struct RigidBody {
    pub absolute_position: Vector3,
    pub velocity: Vector3,
    pub mesh: Box<dyn MeshShape>,
}

impl Dynamic for RigidBody {
    fn move_by(&mut self, change: Vector3) {
        self.absolute_position += change;
        self.mesh.move_by(change);
    }
}

impl RigidBody {
    pub fn new(absolute_position: Vector3, mesh: Box<dyn MeshShape>) -> Self {
        RigidBody {
            absolute_position,
            velocity: Vector3::new(0.0, 0.0, 0.0),
            mesh,
        }
    }
}

fn separating_axis_theorem(polygon1: &Polygon, polygon2: &Polygon) -> bool {
    for i in 0..polygon1.points.len() {
        let line_segment1 = polygon1.get_edges()[i];
        let third_point = polygon1.get_edges()[(i + 1) % polygon1.points.len()].1; // Borrow another point so we can calculate the normal vector to the line

        let projection_line = calculate_orthogonal_line(&line_segment1, &third_point);

        let line_segment2 = project_polygon_onto_line(&projection_line, polygon2);

        if !line_segments_overlap(line_segment1, line_segment2) {
            return false;
        }
    }

    for i in 0..polygon2.points.len() {
        let line_segment1 = polygon2.get_edges()[i];
        let third_point = polygon2.get_edges()[(i + 1) % polygon1.points.len()].1;

        let projection_line = calculate_orthogonal_line(&line_segment1, &third_point);

        let line_segment2 = project_polygon_onto_line(&projection_line, polygon2);

        if !line_segments_overlap(line_segment1, line_segment2) {
            return false;
        }
    }

    return true;
}
pub fn rigid_bodies_collide(rigid_body1: &RigidBody, rigid_body2: &RigidBody) -> bool {
    let self_polygons = rigid_body1.mesh.get_polygons();
    let other_polygons = rigid_body2.mesh.get_polygons();
    for polygon1 in &self_polygons {
        let projection_planes = calculate_orthogonal_planes(&polygon1);
        let projected_polygons: Vec<Polygon> = projection_planes
            .iter()
            .map(|plane| project_mesh_onto_plane(plane, &other_polygons))
            .collect();

        // Perform the separating axis theorem for each polygon
        for polygon2 in projected_polygons.iter() {
            if !separating_axis_theorem(&polygon1, polygon2) {
                return false;
            }
        }
    }

    // We have to the same procedure with the other shape according to the algorithm
    for polygon1 in other_polygons {
        let projection_planes = calculate_orthogonal_planes(&polygon1);
        let projected_polygons: Vec<Polygon> = projection_planes
            .iter()
            .map(|plane| project_mesh_onto_plane(plane, &self_polygons))
            .collect();

        // Perform the separating axis theorem for each polygon
        for polygon2 in projected_polygons.iter() {
            if !separating_axis_theorem(&polygon1, polygon2) {
                return false;
            }
        }
    }

    return true;
}

pub fn polygons_collide(polygon1: &Polygon, polygon2: &Polygon) -> bool {
    // This is why calc 3 was useful :)

    let plane1 = convert_polygon_to_plane(polygon1);
    let plane2 = convert_polygon_to_plane(polygon2);

    // If this is false, then the polygons are parallel and won't collide
    if !check_planes_intersect(&plane1, &plane2) {
        return false;
    }

    // Now we can calculate the line in 3 dimensions where these planes intersect
    let plane_intersection = calculate_line_intersection_between_planes(&plane1, &plane2);

    match plane_intersection {
        // If the polygons are to intersect, then this intersection line should be contained by both polygons so just use polygon 1
        PlaneIntersection::Line(line) => {
            // Use the line's parallel vector to get the bounding angles
            let bounding_angles = calculate_bounding_angles(&line.v, polygon1);

            return check_angles_surround_relative_origin(bounding_angles);
        }
        PlaneIntersection::Infinite => {
            // Immediately return true
            return true;
        }
    }
}

fn calculate_bounding_angles(base_vec: &Vector3, polygon: &Polygon) -> Vec<SphericalAngle> {
    let mut bounding_vecs: Vec<Vector3> = Vec::with_capacity(polygon.points.len());
    for vertex in polygon.points.iter() {
        let bounding_vec = Vector3::new(
            vertex.x - base_vec.x,
            vertex.y - base_vec.y,
            vertex.z - base_vec.z,
        );

        bounding_vecs.push(bounding_vec);
    }

    // Find the difference in angles
    let bounding_angles: Vec<SphericalAngle> = bounding_vecs
        .iter()
        .map(|vec| calculate_difference_in_angle(base_vec, vec))
        .collect();

    return bounding_angles;
}
