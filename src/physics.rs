use raylib::prelude::*;

pub trait DetectCollision {
    fn get_previous_position(&self) -> Vector3;
    fn get_current_position(&self) -> Vector3;
    fn get_collision_vertices(&self) -> Vec<Vector3>;

    fn colliding_with<T: DetectCollision>(&self, other: &T) -> bool {
        todo!()
    }
}

pub struct PhysicsBody {
    // Both of these positions are measured from the center of the object (whatever that means)
    // The previous position is necessary for determining collisions
    previous_position: Vector3,
    current_position: Vector3,

    // These are the vertices that define the shape, and they should be relative to the center
    pub collision_vertices: Vec<Vector3>,
}

impl PhysicsBody {
    pub fn new(
        collision_vertices: Vec<Vector3>,
        previous_position: Vector3,
        current_position: Vector3,
    ) -> Self {
        PhysicsBody {
            collision_vertices,
            previous_position,
            current_position,
        }
    }
}

impl DetectCollision for PhysicsBody {
    fn get_previous_position(&self) -> Vector3 {
        self.previous_position
    }

    fn get_current_position(&self) -> Vector3 {
        self.current_position
    }

    fn get_collision_vertices(&self) -> Vec<Vector3> {
        // I know this isn't the best, but for now to get things running this should be alright
        self.collision_vertices.clone()
    }
}
