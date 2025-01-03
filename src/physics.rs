use crate::hashable_vector::HashableVector3;
use crate::math_util::{Line3D, Plane};
use crate::render::*;
use raylib::prelude::*;
use std::collections::HashSet;

pub trait Physical {
    fn get_center(&self) -> Vector3;
    fn get_bounding_circle_radius(&self) -> f32;
    fn get_mesh(&self) -> &Box<dyn MeshShape>;
}

pub trait Dynamic {
    // Moves the object and its components by change
    fn move_by(&mut self, change: Vector3);
}

pub trait Static {}

#[derive(Debug)]
pub struct StaticBody {
    pub absolute_position: Vector3,
    pub mesh: Box<dyn MeshShape>,
}

impl Static for StaticBody {}

impl Physical for StaticBody {
    fn get_mesh(&self) -> &Box<dyn MeshShape> {
        return &self.mesh;
    }
    fn get_bounding_circle_radius(&self) -> f32 {
        return self.mesh.get_bounding_circle_radius();
    }
    fn get_center(&self) -> Vector3 {
        return self.mesh.get_center();
    }
}

impl StaticBody {
    pub fn new(absolute_position: Vector3, mesh: Box<dyn MeshShape>) -> Self {
        Self {
            absolute_position,
            mesh,
        }
    }
}

// This object is a basic object that interacts with physics
#[derive(Debug)]
pub struct DynamicBody {
    pub absolute_position: Vector3,
    pub velocity: Vector3,
    pub mesh: Box<dyn MeshShape>,
}

impl Dynamic for DynamicBody {
    fn move_by(&mut self, change: Vector3) {
        self.absolute_position += change;
        self.mesh.move_by(change);
    }
}

impl Physical for DynamicBody {
    fn get_center(&self) -> Vector3 {
        return self.mesh.get_center();
    }

    fn get_bounding_circle_radius(&self) -> f32 {
        return self.mesh.get_bounding_circle_radius();
    }

    fn get_mesh(&self) -> &Box<dyn MeshShape> {
        return &self.mesh;
    }
}

impl DynamicBody {
    pub fn new(absolute_position: Vector3, mesh: Box<dyn MeshShape>) -> Self {
        DynamicBody {
            absolute_position,
            velocity: Vector3::new(0.0, 0.0, 0.0),
            mesh,
        }
    }

    /**
    Returns the mtv to snap the objects onto the edges of each other while not intersecting
    If the objects aren't colliding, then returns `None`
    The MTV should be applied to the dynamic body (root body) of the collision because it's designed to move it out of the other body
     */
    pub fn collides_with(&self, other: Box<&dyn Physical>) -> Option<Vector3> {
        // We need to find the angle when we look at there is a minimum overlap
        enum SATResult {
            MTV(f32, Vector3),
            AllAxesChecked,
            NoCollision,
        }

        fn separating_axis_theorem_3d(
            mesh1: &Box<dyn MeshShape>,
            mesh2: &Box<dyn MeshShape>,
            checked_planes: &mut HashSet<HashableVector3>,
        ) -> SATResult {
            let mut overlap: f32 = f32::MAX;
            let mut axis: Option<Vector3> = None;
            let first_mesh_polygons = mesh1.get_polygons();
            for polygon1 in &first_mesh_polygons {
                let projection_plane = polygon1.calculate_orthogonal_plane();
                let projected_polygon1 = projection_plane.project_mesh(mesh1);
                let projected_polygon2 = projection_plane.project_mesh(mesh2);

                // Determine if we've already seen this
                let mut direction = projection_plane.n.normalized();
                if direction.y < 0.0 {
                    direction *= -1.0;
                }
                if checked_planes.contains(&HashableVector3::from_vector3(direction)) {
                    continue;
                }

                let collision = projected_polygon1.collides_with(&projected_polygon2);
                match collision {
                    Some(collision) => {
                        if collision.0 < overlap {
                            overlap = collision.0;
                            axis = Some(collision.1);
                        }
                    }
                    None => {
                        return SATResult::NoCollision;
                    }
                }

                // We haven't seen it, so make sure it's clear that we've checked it
                checked_planes.insert(HashableVector3::from_vector3(direction));
            }

            match axis {
                Some(axis) => {
                    return SATResult::MTV(overlap, axis);
                }
                None => {
                   return SATResult::AllAxesChecked;
                }
            }
        }

        fn align_direction_vec(direction_vec: &mut Vector3, dynamic_body_mesh: &Box<dyn MeshShape>, other_mesh: &Box<dyn MeshShape>) {
            // The direction vector is going to either be correct, or flipped, so walk in each direction, and the one that gets further from the other mesh's center is the right direction
            let dynamic_body_center = dynamic_body_mesh.get_center();
            let other_body_center = other_mesh.get_center();
            
            let center1 = dynamic_body_center + *direction_vec;
            let center2 = dynamic_body_center - *direction_vec;

            // If adding the direction vector walks in the wrong direction
            if (other_body_center - center1).length() < (other_body_center - center2).length() {
                *direction_vec *= -1.0;
            }
        }

        let mut min_overlap = 0.0;
        let mut direction: Vector3 = Vector3::new(0.0, 0.0, 0.0,);
        let mut checked_planes: HashSet<HashableVector3> = HashSet::new();
        let collision1 =
            separating_axis_theorem_3d(&self.mesh, other.get_mesh(), &mut checked_planes);
        match collision1 {
            SATResult::MTV(overlap, direction_vec) => {
                min_overlap = overlap;

                // This is going to have to be reversed based on the direction it's facing, but that's a later problem, right now I just want to put the code down
                direction = direction_vec;
                // TODO
            }
            SATResult::NoCollision => {
                return None;
            }
            SATResult::AllAxesChecked => {}
        }

        let collision2 =
            separating_axis_theorem_3d(other.get_mesh(), &self.mesh, &mut checked_planes);
        match collision2 {
            SATResult::MTV(overlap, direction_vec) => {
                if overlap < min_overlap {
                    min_overlap = overlap;

                    // This is going to have to be reversed based on the direction it's facing, but that's a later problem, right now I just want to put the code down
                    direction = direction_vec;
                    // TODO
                }
            }
            SATResult::NoCollision => {
                return None;
            }
            SATResult::AllAxesChecked => {}
        }

        align_direction_vec(&mut direction, self.get_mesh(), other.get_mesh());

        return Some(direction.normalized() * min_overlap);
    }
}
