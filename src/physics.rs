use crate::math_util::*;
use raylib::prelude::*;

pub struct Polygon {
    pub points: Vec<Vector3>,
}

impl Polygon {
    pub fn new(points: Vec<Vector3>) -> Self {
        Polygon { points }
    }

    pub fn as_plane(&self) -> Plane {
        // Just make it the first point because why not, we just need *a* point
        let p0 = self.points.first().unwrap();

        // Calculate cross product

        // Every polygon should have three points, so should be able to safely retrieve the second
        // and third points in the vector
        let p1 = self.points.get(1).unwrap();
        let p2 = self.points.get(2).unwrap();

        let v1 = Vector3::new(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);

        let v2 = Vector3::new(p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);

        let cross = cross_product(v1, v2);

        Plane::from_point_and_normal(*p0, cross)
    }
}

pub trait MeshShape {
    fn get_vertices(&self) -> Vec<Vector3>;
    fn get_polygons(&self) -> Vec<Polygon>;
}

pub struct RectangularPrism {
    // The root is the bottom left point of the front face
    pub root: Vector3,

    pub length: f32,
    pub width: f32,
    pub height: f32,
}

impl MeshShape for RectangularPrism {
    fn get_vertices(&self) -> Vec<Vector3> {
        vec![
            Vector3::new(self.root.x, self.root.y, self.root.z),
            Vector3::new(self.root.x + self.width, self.root.y, self.root.z),
            Vector3::new(self.root.x, self.root.y + self.height, self.root.z),
            Vector3::new(
                self.root.x + self.width,
                self.root.y + self.height,
                self.root.z,
            ),
            Vector3::new(self.root.x, self.root.y, self.root.z + self.length),
            Vector3::new(
                self.root.x + self.width,
                self.root.y,
                self.root.z + self.length,
            ),
            Vector3::new(
                self.root.x,
                self.root.y + self.height,
                self.root.z + self.length,
            ),
            Vector3::new(
                self.root.x + self.width,
                self.root.y + self.height,
                self.root.z + self.length,
            ),
        ]
    }

    fn get_polygons(&self) -> Vec<Polygon> {
        let vertices = self.get_vertices();

        // Define all of the polygons
        vec![
            Polygon::new(vec![vertices[0], vertices[1], vertices[2], vertices[3]]),
            Polygon::new(vec![vertices[4], vertices[5], vertices[6], vertices[7]]),
            Polygon::new(vec![vertices[0], vertices[1], vertices[4], vertices[5]]),
            Polygon::new(vec![vertices[2], vertices[3], vertices[6], vertices[7]]),
            Polygon::new(vec![vertices[0], vertices[2], vertices[4], vertices[6]]),
            Polygon::new(vec![vertices[1], vertices[3], vertices[5], vertices[7]]),
        ]
    }
}

pub struct PhysicsBody<'a> {
    // Both of these positions are measured from the root of the object (whatever that means)
    // The previous position is necessary for determining collisions
    pub previous_position: Vector3,
    pub current_position: Vector3,

    pub simple_shape: &'a dyn MeshShape,
}

impl<'a> PhysicsBody<'a> {
    pub fn new(
        previous_position: Vector3,
        current_position: Vector3,
        simple_shape: &'a dyn MeshShape,
    ) -> Self {
        PhysicsBody {
            previous_position,
            current_position,
            simple_shape,
        }
    }
}

fn check_polygons_collide(polygon1: &Polygon, polygon2: &Polygon) -> bool {
    // This is why calc 3 was useful :)

    let plane1 = polygon1.as_plane();
    let plane2 = polygon2.as_plane();

    // If this is false, then the polygons are parallel and won't collide
    if !check_planes_intersect(&plane1, &plane2) {
        return false;
    }

    // Now we can calculate the line in 3 dimensions where these planes intersect
    let plane_intersection = calculate_line_intersection_between_planes(&plane1, &plane2);
    
    match plane_intersection {
        PlaneIntersection::Line(line) => {
            // If the polygons are to intersect, then this intersection line should be contained by both polygons

            // We need to check if the line and a plane intersect because if it intersects with one plane, then it intersects with both
        }
        PlaneIntersection::VerticalLine(vertical_line) => {

        }
        PlaneIntersection::Infinite => {

        }
    }

    todo!()
}
