use crate::hashable_vector::HashableVector3;
use crate::render::*;
use raylib::prelude::*;
use std::collections::HashSet;

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

    pub fn collides_with(&self, other: &RigidBody) -> bool {
        let mut checked_planes: HashSet<HashableVector3> = HashSet::new();
        let self_polygons = self.mesh.get_polygons();
        let other_polygons = other.mesh.get_polygons();
        for polygon1 in &self_polygons {
            let projection_planes = polygon1.calculate_orthogonal_planes();
            let projected_polygons1: Vec<Polygon> = projection_planes
                .iter()
                .map(|plane| plane.project_mesh(&self.mesh))
                .collect();
            let projected_polygons2: Vec<Polygon> = projection_planes
                .iter()
                .map(|plane| plane.project_mesh(&other.mesh))
                .collect();

            // Perform the separating axis theorem for each polygon
            for i in 0..projection_planes.len() {
                // Determine if we've already seen this
                let projection_plane = &projection_planes[i];
                let mut direction = projection_plane.n.normalized();
                if direction.y < 0.0 {
                    direction *= -1.0;
                }
                if checked_planes.contains(&HashableVector3::from_vector3(direction)) {
                    continue;
                }

                if !projected_polygons1[i].collides_with(&projected_polygons2[i]) {
                    return false;
                }

                // We haven't seen it, so make sure it's clear that we've checked it
                checked_planes.insert(HashableVector3::from_vector3(direction));
            }
        }

        // We have to the same procedure with the other shape according to the algorithm
        for polygon1 in &other_polygons {
            let projection_planes = polygon1.calculate_orthogonal_planes();
            let projected_polygons1: Vec<Polygon> = projection_planes
                .iter()
                .map(|plane| plane.project_mesh(&self.mesh))
                .collect();
            let projected_polygons2: Vec<Polygon> = projection_planes
                .iter()
                .map(|plane| plane.project_mesh(&other.mesh))
                .collect();
            // Perform the separating axis theorem for each polygon
            for i in 0..projection_planes.len() {
                // Determine if we've already seen this
                let projection_plane = &projection_planes[i];
                let mut direction = projection_plane.n.normalized();
                if direction.y < 0.0 {
                    direction *= -1.0;
                }
                if checked_planes.contains(&HashableVector3::from_vector3(direction)) {
                    continue;
                }

                if !projected_polygons1[i].collides_with(&projected_polygons2[i]) {
                    return false;
                }

                // We haven't seen it, so make sure it's clear that we've checked it
                checked_planes.insert(HashableVector3::from_vector3(direction));
            }
        }

        return true;
    }

    pub fn get_center(&self) -> Vector3 {
        return self.mesh.get_center();
    }

    pub fn get_bounding_circle_radius(&self) -> f32 {
        return self.mesh.get_bounding_circle_radius();
    }
}
